
'''
The ``cornish-view`` daemon (SPEC-09A §5).

Stdlib only (D2): ``http.server.ThreadingHTTPServer`` with one thread per
connection. An SSE connection holds its thread for its lifetime, which is fine
at local-viewer scale — the documented ceiling is "a handful of browser tabs",
and aiohttp remains the upgrade path if that ever binds.

All state lives in one :class:`~cornish.view.registry.PlotRegistry`; the
handlers here are thin shims over it, and the Flask blueprint drives the same
registry through the same helpers.

**Threat model (§5.4), stated rather than discovered.** The socket binds
127.0.0.1 only, which is NOT the same as single-user on a shared host: any local
account can connect to the loopback interface. Writes are protected — the token
file is mode 0600, so only the owner can read it — but **reads are open by
design** (the browser's ``EventSource`` API cannot send an ``Authorization``
header, so gating SSE on the token would break the viewer). On a multi-user
machine any local user can therefore view any plot: ``GET /view/{id}``,
``.json``, ``.svg`` and ``/events/``. That is acceptable for the target use —
plots of your own analysis on your own machine — but do not route sensitive plot
content through the daemon on a shared host. Opt-in read auth is the documented
upgrade path if that constraint ever binds. (Loopback traffic is not exposed to
the network; the open-reads *port* is the concern, not sniffing.)
'''

import argparse
import json
import logging
import os
import queue
import re
import secrets
import select
import socket as socket_module
import stat
import sys
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Optional

from ..exc import EnvelopeError
from ..plot.envelope import SCHEMA_VERSION, canonical_json, envelope_hash
from .registry import PlotRegistry

__all__ = ['DEFAULT_PORT', 'ViewServer', 'serve', 'main', 'ensure_daemon',
           'token_path', 'read_or_create_token']

logger = logging.getLogger("cornish")

DEFAULT_PORT = 8765

#: v1 hard-codes loopback. Exposing plots beyond localhost is a deliberate
#: non-feature until someone needs it and thinks about it.
BIND_HOST = "127.0.0.1"

#: SVG width bounds for the query parameter (§5.2).
MINIMUM_SVG_WIDTH = 64
MAXIMUM_SVG_WIDTH = 4096

#: Keep-alive interval; detects dead connections through proxies.
PING_SECONDS = 15.0

SSE_RETRY_MS = 2000

STATIC_DIRECTORY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")

# The id group is LAZY on purpose. The legal plot-id charset (§3.1) contains
# '.', so a greedy group swallows the extension and '/plot/demo.json' resolves
# to the id 'demo.json' with no suffix — which silently 404s the document route.
# Lazy matching resolves the inherent ambiguity in favour of the suffix, which
# is what `GET /plot/{plot_id}.json` means.
_PLOT_ROUTE = re.compile(r"^/plot/([A-Za-z0-9._-]{1,64}?)(\.json|\.svg)?$")
_VIEW_ROUTE = re.compile(r"^/view/([A-Za-z0-9._-]{1,64})$")
_EVENTS_ROUTE = re.compile(r"^/events/([A-Za-z0-9._-]{1,64})$")
_STATIC_ROUTE = re.compile(r"^/static/([A-Za-z0-9._-]{1,64})$")


# ---------------------------------------------------------------------------
# token handling (§5.4)
# ---------------------------------------------------------------------------

def token_path() -> str:
	''' ``~/.cornish/view_token``. '''
	return os.path.join(os.path.expanduser("~"), ".cornish", "view_token")


def read_or_create_token(path: Optional[str] = None) -> str:
	'''
	Read the bearer token, creating one if absent.

	The directory is created 0700 and the file 0600: the token is the only thing
	standing between a local account and write access to the daemon.
	'''
	path = path or token_path()
	directory = os.path.dirname(path)
	if directory:
		os.makedirs(directory, mode=0o700, exist_ok=True)
	if os.path.exists(path):
		with open(path, encoding="utf-8") as stream:
			existing = stream.read().strip()
		if existing:
			return existing

	token = secrets.token_urlsafe(32)
	descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
	with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
		stream.write(token + "\n")
	os.chmod(path, stat.S_IRUSR | stat.S_IWUSR)
	return token


# ---------------------------------------------------------------------------
# the handler
# ---------------------------------------------------------------------------

class _Handler(BaseHTTPRequestHandler):
	''' Thin shims over the registry; the blueprint mirrors these route by route. '''

	server_version = "cornish-view"
	protocol_version = "HTTP/1.1"

	# -- plumbing --------------------------------------------------------

	@property
	def registry(self) -> PlotRegistry:
		return self.server.registry

	def log_message(self, format, *args):
		logger.debug("cornish-view %s - %s", self.address_string(), format % args)

	def _send(self, status, body: bytes = b"", content_type: str = "text/plain; charset=utf-8",
	          extra_headers=None) -> None:
		self.send_response(int(status))
		self.send_header("Content-Type", content_type)
		self.send_header("Content-Length", str(len(body)))
		for name, value in (extra_headers or {}).items():
			self.send_header(name, value)
		self.end_headers()
		if self.command != "HEAD" and body:
			self.wfile.write(body)

	def _send_json(self, status, payload: dict, extra_headers=None) -> None:
		self._send(status, json.dumps(payload).encode("utf-8"),
		           "application/json; charset=utf-8", extra_headers)

	def _send_error_json(self, status, message: str) -> None:
		self._send_json(status, {"error": message})

	def _authorize(self) -> bool:
		'''
		Bearer-token check for mutating routes. Auth runs FIRST in the frozen
		precedence so an unauthenticated caller learns nothing about how bodies
		are handled.
		'''
		if self.server.token is None:
			return True
		header = self.headers.get("Authorization", "")
		if not header.startswith("Bearer "):
			self._send_error_json(HTTPStatus.UNAUTHORIZED,
			                      "A bearer token is required: Authorization: Bearer <token>.")
			return False
		presented = header[len("Bearer "):].strip()
		if not secrets.compare_digest(presented, self.server.token):
			self._send_error_json(HTTPStatus.FORBIDDEN, "The bearer token is not valid.")
			return False
		return True

	def _read_json_body(self):
		'''
		Read and parse the request body.

		:returns: ``(payload, True)`` on success, ``(None, False)`` when a response was already sent
		'''
		content_type = (self.headers.get("Content-Type") or "").split(";")[0].strip().lower()
		if content_type != "application/json":
			self._send_error_json(HTTPStatus.UNSUPPORTED_MEDIA_TYPE,
			                      "Content-Type must be application/json; "
			                      f"got '{content_type or 'nothing'}'.")
			return None, False
		try:
			length = int(self.headers.get("Content-Length", "0"))
		except ValueError:
			self._send_error_json(HTTPStatus.BAD_REQUEST, "Content-Length is not an integer.")
			return None, False
		raw = self.rfile.read(length) if length > 0 else b""
		try:
			return json.loads(raw.decode("utf-8")), True
		except (UnicodeDecodeError, json.JSONDecodeError) as error:
			self._send_error_json(HTTPStatus.BAD_REQUEST, f"The body is not valid JSON: {error}")
			return None, False

	# -- routes ----------------------------------------------------------

	def do_GET(self):
		path, _, query = self.path.partition("?")

		if path == "/healthz":
			self._send_json(HTTPStatus.OK,
			                {"cornish_view": SCHEMA_VERSION, "plots": len(self.registry.plot_ids())})
			return

		if path == "/":
			self._send(HTTPStatus.OK, self._index_html().encode("utf-8"),
			           "text/html; charset=utf-8")
			return

		match = _VIEW_ROUTE.match(path)
		if match:
			# an unknown id still returns 200: the page subscribes and waits,
			# because opening a viewer before the first POST is supported
			self._send_static("viewer.html", "text/html; charset=utf-8")
			return

		match = _STATIC_ROUTE.match(path)
		if match:
			name = match.group(1)
			content_type = {"js": "text/javascript; charset=utf-8",
			                "mjs": "text/javascript; charset=utf-8",
			                "html": "text/html; charset=utf-8",
			                "json": "application/json; charset=utf-8"}.get(
				name.rsplit(".", 1)[-1], "application/octet-stream")
			self._send_static(name, content_type)
			return

		match = _EVENTS_ROUTE.match(path)
		if match:
			self._stream_events(match.group(1))
			return

		match = _PLOT_ROUTE.match(path)
		if match:
			plot_id, suffix = match.group(1), match.group(2)
			if suffix == ".json":
				self._send_document(plot_id)
				return
			if suffix == ".svg":
				self._send_svg(plot_id, query)
				return

		self._send_error_json(HTTPStatus.NOT_FOUND, f"No route for {path}.")

	def do_POST(self):
		path = self.path.partition("?")[0]
		match = _PLOT_ROUTE.match(path)
		if not match or match.group(2):
			self._send_error_json(HTTPStatus.NOT_FOUND, f"No route for POST {path}.")
			return
		plot_id = match.group(1)

		# frozen precedence: auth -> Content-Type -> parse/validate -> mismatch -> state
		if not self._authorize():
			return
		payload, ok = self._read_json_body()
		if not ok:
			return

		try:
			revision = self.registry.apply(plot_id, payload)
		except EnvelopeError as error:
			self._send_error_json(HTTPStatus.BAD_REQUEST, str(error))
		except ValueError as error:
			self._send_error_json(HTTPStatus.CONFLICT, str(error))
		except KeyError as error:
			self._send_error_json(HTTPStatus.NOT_FOUND, str(error).strip("'"))
		else:
			self._send_json(HTTPStatus.OK, {"revision": revision})

	def do_DELETE(self):
		path = self.path.partition("?")[0]
		match = _PLOT_ROUTE.match(path)
		if not match or match.group(2):
			self._send_error_json(HTTPStatus.NOT_FOUND, f"No route for DELETE {path}.")
			return
		if not self._authorize():
			return
		try:
			self.registry.delete_plot(match.group(1))
		except KeyError as error:
			self._send_error_json(HTTPStatus.NOT_FOUND, str(error).strip("'"))
		else:
			self._send_json(HTTPStatus.OK, {"deleted": match.group(1)})

	# -- route bodies ----------------------------------------------------

	def _send_static(self, name: str, content_type: str) -> None:
		path = os.path.normpath(os.path.join(STATIC_DIRECTORY, name))
		if not path.startswith(STATIC_DIRECTORY) or not os.path.isfile(path):
			self._send_error_json(HTTPStatus.NOT_FOUND, f"No static asset '{name}'.")
			return
		with open(path, "rb") as stream:
			body = stream.read()
		self._send(HTTPStatus.OK, body, content_type, {"Cache-Control": "no-cache"})

	def _send_document(self, plot_id: str) -> None:
		try:
			document = self.registry.document(plot_id)
		except KeyError:
			self._send_error_json(HTTPStatus.NOT_FOUND, f"No plot '{plot_id}'.")
			return
		body = canonical_json(document).encode("utf-8")
		self._send(HTTPStatus.OK, body, "application/json; charset=utf-8",
		           {"ETag": f'"{envelope_hash(document)}"'})

	def _send_svg(self, plot_id: str, query: str) -> None:
		width = 640
		if query:
			from urllib.parse import parse_qs
			values = parse_qs(query).get("width")
			if values:
				try:
					width = int(values[0])
				except ValueError:
					self._send_error_json(HTTPStatus.BAD_REQUEST,
					                      f"width must be an integer; got '{values[0]}'.")
					return
				if not MINIMUM_SVG_WIDTH <= width <= MAXIMUM_SVG_WIDTH:
					# never silently clamped: the caller asked for a size the
					# server will not produce, and should be told so
					self._send_error_json(
						HTTPStatus.BAD_REQUEST,
						f"width must lie in [{MINIMUM_SVG_WIDTH}, {MAXIMUM_SVG_WIDTH}]; got {width}.")
					return
		try:
			document = self.registry.document(plot_id)
		except KeyError:
			self._send_error_json(HTTPStatus.NOT_FOUND, f"No plot '{plot_id}'.")
			return
		from ..plot.svg import render_svg
		body = render_svg(document, width=width).encode("utf-8")
		# width changes the bytes, so it must change the ETag
		self._send(HTTPStatus.OK, body, "image/svg+xml; charset=utf-8",
		           {"ETag": f'"{envelope_hash(document)}-w{width}"'})

	def _index_html(self) -> str:
		items = "".join(f'<li><a href="/view/{plot_id}">{plot_id}</a></li>'
		                for plot_id in self.registry.plot_ids())
		return ("<!doctype html><html><head><meta charset=\"utf-8\">"
		        "<title>cornish-view</title></head><body><h1>cornish-view</h1>"
		        f"<ul>{items or '<li><em>no plots yet</em></li>'}</ul></body></html>")

	def _stream_events(self, plot_id: str) -> None:
		'''
		The SSE stream (§5.3): ``retry:``, one snapshot, then deltas.

		Snapshot-then-deltas makes the protocol stateless per connection, so an
		EventSource auto-reconnect mid-stream is always consistent.
		'''
		subscription = self.registry.subscribe(plot_id)
		try:
			self.send_response(int(HTTPStatus.OK))
			self.send_header("Content-Type", "text/event-stream; charset=utf-8")
			self.send_header("Cache-Control", "no-cache")
			self.send_header("Connection", "keep-alive")
			self.end_headers()
			self.wfile.write(f"retry: {SSE_RETRY_MS}\n\n".encode("utf-8"))
			self.wfile.flush()

			last_ping = time.monotonic()
			while not self.server.stopping:
				if self._peer_has_gone():
					logger.debug("cornish-view: viewer for '%s' closed the connection.", plot_id)
					return
				if subscription.overflowed:
					# drop the connection so the viewer reconnects and takes a
					# clean snapshot, rather than silently missing events
					logger.warning("cornish-view: subscriber for '%s' fell more than %d events "
					               "behind; dropping the connection so it reconnects cleanly.",
					               plot_id, subscription.queue.maxsize)
					return
				try:
					event, payload = subscription.queue.get(timeout=1.0)
				except queue.Empty:
					if time.monotonic() - last_ping >= PING_SECONDS:
						self.wfile.write(b": ping\n\n")
						self.wfile.flush()
						last_ping = time.monotonic()
					continue
				self._write_event(event, payload)
				last_ping = time.monotonic()
		except (BrokenPipeError, ConnectionResetError):
			logger.debug("cornish-view: viewer for '%s' disconnected.", plot_id)
		finally:
			self.registry.unsubscribe(subscription)

	def _peer_has_gone(self) -> bool:
		'''
		Whether the viewer has closed the connection.

		A failed *write* is not a reliable signal: writing to a half-closed
		socket routinely succeeds (the bytes go to the kernel buffer and the RST
		arrives later), so a departed viewer can otherwise hold its subscription
		until the next keep-alive ping 15 seconds later — verified on macOS.
		Polling the read side is definite: an SSE client sends nothing after its
		request, so a readable socket means EOF.

		If the peer did send something, this consumes it — acceptable because
		the connection is dedicated to the event stream and nothing else is
		read from it.
		'''
		try:
			readable, _, _ = select.select([self.connection], [], [], 0)
			if not readable:
				return False
			return self.connection.recv(1, socket_module.MSG_PEEK) == b""
		except OSError:
			return True

	def _write_event(self, event: str, payload: dict) -> None:
		'''
		Write one SSE event. Canonical JSON is single-line, so this is always
		one ``data:`` line — but a conformant writer must still split on any
		embedded newline, so it does. It simply never triggers.
		'''
		serialized = canonical_json(payload)
		lines = "".join(f"data: {line}\n" for line in serialized.split("\n"))
		self.wfile.write(f"event: {event}\n{lines}\n".encode("utf-8"))
		self.wfile.flush()


class ViewServer(ThreadingHTTPServer):
	''' The daemon. Started by :func:`serve`; stop it with :meth:`shutdown_now`. '''

	daemon_threads = True
	allow_reuse_address = True

	def __init__(self, address, handler, registry: PlotRegistry, token: Optional[str]):
		super().__init__(address, handler)
		self.registry = registry
		self.token = token
		self.stopping = False
		self._thread = None

	def handle_error(self, request, client_address):
		'''
		Route handler exceptions to the logger instead of printing a traceback
		to stderr.

		A browser closing a tab mid-stream raises ``ConnectionResetError`` /
		``BrokenPipeError`` here; that is normal operation, not a fault, and
		socketserver's default of dumping a traceback makes routine use look
		broken. **Everything else is still reported** — at error level, with the
		traceback — because a genuine handler bug must not be quieted by the
		same change.
		'''
		error = sys.exc_info()[1]
		if isinstance(error, (BrokenPipeError, ConnectionResetError)):
			logger.debug("cornish-view: client %s disconnected (%s).",
			             client_address, type(error).__name__)
			return
		logger.exception("cornish-view: unhandled error serving %s", client_address)

	@property
	def port(self) -> int:
		return self.server_address[1]

	@property
	def base_url(self) -> str:
		return f"http://{BIND_HOST}:{self.port}"

	def shutdown_now(self) -> None:
		''' Stop serving and release the socket. '''
		self.stopping = True
		self.shutdown()
		self.server_close()
		if self._thread is not None:
			self._thread.join(timeout=5)

	#: alias matching the spec's ``.shutdown()`` naming; the base class already
	#: defines shutdown(), so the explicit name above is the one to call
	def stop(self) -> None:
		self.shutdown_now()


def serve(registry: Optional[PlotRegistry] = None, *, port: int = DEFAULT_PORT,
          token: Optional[str] = None) -> ViewServer:
	'''
	Start the daemon on a background thread and return it.

	:param registry: the registry to serve; a fresh one is created when omitted
	:param port: TCP port; 0 picks a free one (useful in tests)
	:param token: bearer token for writes; ``None`` disables write auth entirely
	'''
	server = ViewServer((BIND_HOST, port), _Handler, registry or PlotRegistry(), token)
	thread = threading.Thread(target=server.serve_forever, name="cornish-view", daemon=True)
	server._thread = thread
	thread.start()
	return server


def main(argv: Optional[list] = None) -> None:
	''' ``python -m cornish.view`` — run the daemon in the foreground. '''
	parser = argparse.ArgumentParser(prog="cornish.view",
	                                 description="Serve cornish plots to a browser (localhost only).")
	parser.add_argument("-p", "--port", type=int, default=DEFAULT_PORT)
	parser.add_argument("--token-file", default=None,
	                    help="path to the bearer-token file (default: ~/.cornish/view_token)")
	arguments = parser.parse_args(argv)

	token = read_or_create_token(arguments.token_file)
	server = serve(PlotRegistry(), port=arguments.port, token=token)
	print(f"cornish-view listening on {server.base_url} (localhost only)", file=sys.stderr)
	print(f"token file: {arguments.token_file or token_path()}", file=sys.stderr)
	try:
		while True:
			time.sleep(3600)
	except KeyboardInterrupt:
		server.shutdown_now()


def ensure_daemon(*, port: int = DEFAULT_PORT, launch: bool = True) -> tuple:
	'''
	Return ``(base_url, token)`` for a running daemon, starting one if needed.

	Probes ``GET /healthz``; if nothing answers and *launch* is set, spawns
	``sys.executable -m cornish.view`` detached and polls until it is ready.

	:raises RuntimeError: if no daemon is reachable within the timeout
	'''
	import subprocess
	import urllib.error
	import urllib.request

	base_url = f"http://{BIND_HOST}:{port}"

	def healthy() -> bool:
		try:
			with urllib.request.urlopen(f"{base_url}/healthz", timeout=1.0) as response:
				return response.status == 200
		except (urllib.error.URLError, OSError):
			return False

	token = read_or_create_token()
	if healthy():
		return base_url, token
	if not launch:
		raise RuntimeError(f"No cornish-view daemon is listening on {base_url} "
		                   "and launch=False was requested.")

	subprocess.Popen([sys.executable, "-m", "cornish.view", "--port", str(port)],
	                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
	                 start_new_session=True)

	deadline = time.monotonic() + 5.0
	while time.monotonic() < deadline:
		if healthy():
			return base_url, token
		time.sleep(0.1)
	raise RuntimeError(f"A cornish-view daemon was launched but did not answer on {base_url} "
	                   "within 5 seconds.")
