
'''
The Flask blueprint (SPEC-09A §5.2, §5.6).

This is the **only** module in cornish that imports Flask, and it does so behind
a guarded import with an actionable message. It exposes the same route set as
the standalone daemon over the same :class:`~cornish.view.registry.PlotRegistry`,
so the two share their behaviour rather than merely resembling each other.
'''

import json
import os
from typing import Optional

from ..exc import EnvelopeError
from ..plot.envelope import SCHEMA_VERSION, canonical_json, envelope_hash
from ..plot.svg import render_svg
from .registry import PlotRegistry
from .server import (MAXIMUM_SVG_WIDTH, MINIMUM_SVG_WIDTH, STATIC_DIRECTORY,
                     read_or_create_token)

__all__ = ['flask_blueprint']


def flask_blueprint(registry: Optional[PlotRegistry] = None, *,
                    require_token: bool = True):
	'''
	Build a Flask blueprint exposing the §5.2 routes::

	    app.register_blueprint(flask_blueprint(), url_prefix='/plots')

	:param registry: the registry to serve; a fresh one is created when omitted
	:param require_token: gate the *write* routes on the ``~/.cornish/view_token`` bearer token
	:raises ImportError: if Flask is not installed

	**Worker-model caveats — read these before deploying.**

	1. **SSE pins a worker.** Under WSGI, ``GET /events/…`` holds its worker for
	   the connection's life. The Flask dev server needs ``threaded=True`` (its
	   default since 1.0). Under gunicorn use ``gthread`` (with ``threads``
	   comfortably above the expected viewer count) or ``gevent``/``eventlet``;
	   under plain ``sync`` workers the events route **will** exhaust the worker
	   pool. If that is a problem, use pattern 1, or run the standalone daemon
	   beside the app and have the app POST to it — cross-process works because
	   the protocol is already HTTP.

	2. **The registry is per-process.** With ``workers > 1`` a POST lands in one
	   worker's registry while a viewer may be subscribed in another. The
	   blueprint is therefore only correct with a single web process (or process
	   affinity), or with the standalone daemon owning all plot state. This is an
	   architectural property, not a bug queued for later.

	3. **Blueprint auth default.** ``require_token=True`` uses the same
	   ``~/.cornish/view_token``. An app with its own auth layer passes
	   ``require_token=False`` *and takes responsibility*.

	4. **Pattern 1 needs none of this.** ``render_svg`` is a pure call and is the
	   recommended production shape.

	5. **Read routes are unauthenticated BY DESIGN.** The browser's
	   ``EventSource`` API cannot send an ``Authorization`` header, so gating SSE
	   on the bearer token would break the viewer. The *daemon* is safe because
	   it hard-binds 127.0.0.1. **The blueprint inherits the open reads but NOT
	   the binding guarantee: mounting it on a publicly reachable app exposes
	   every plot's full data — documents, SVGs, live SSE — unauthenticated.**
	   The blueprint is for internal/localhost apps. A public app must wrap it in
	   its own auth — one documented line,
	   ``bp.before_request(your_auth_check)``, which can use cookies/sessions
	   (which EventSource *does* send). ``require_token=True`` protects only the
	   write routes, never the reads.
	'''
	try:
		import flask
	except ImportError as error:
		raise ImportError(
			"cornish.view.blueprint requires Flask, which is not installed and is NOT a "
			"cornish dependency. Install it (pip install flask) to embed plots in an "
			"existing app, or use the standalone daemon (python -m cornish.view), or "
			"render server-side with cornish.plot.render_svg — which needs nothing extra."
		) from error

	from flask import Response, request

	registry = registry if registry is not None else PlotRegistry()
	token = read_or_create_token() if require_token else None
	blueprint = flask.Blueprint("cornish_view", __name__)

	def authorize():
		''' ``None`` when the caller may write, else the error response. '''
		if token is None:
			return None
		header = request.headers.get("Authorization", "")
		if not header.startswith("Bearer "):
			return _json_error(flask, 401,
			                   "A bearer token is required: Authorization: Bearer <token>.")
		import secrets
		if not secrets.compare_digest(header[len("Bearer "):].strip(), token):
			return _json_error(flask, 403, "The bearer token is not valid.")
		return None

	@blueprint.route("/healthz", methods=["GET"])
	def healthz():
		return flask.jsonify({"cornish_view": SCHEMA_VERSION, "plots": len(registry.plot_ids())})

	@blueprint.route("/", methods=["GET"])
	def index():
		items = "".join(f'<li><a href="view/{plot_id}">{plot_id}</a></li>'
		                for plot_id in registry.plot_ids())
		return Response("<!doctype html><html><head><meta charset=\"utf-8\">"
		                "<title>cornish-view</title></head><body><h1>cornish-view</h1>"
		                f"<ul>{items or '<li><em>no plots yet</em></li>'}</ul></body></html>",
		                mimetype="text/html")

	@blueprint.route("/view/<plot_id>", methods=["GET"])
	def view(plot_id):
		# an unknown id still returns 200: the page subscribes and waits
		return _static_response(flask, "viewer.html", "text/html")

	@blueprint.route("/static/<path:name>", methods=["GET"])
	def static_asset(name):
		return _static_response(flask, name, None)

	@blueprint.route("/plot/<plot_id>", methods=["POST"])
	def post_plot(plot_id):
		# the frozen precedence, identical to the daemon's
		denied = authorize()
		if denied is not None:
			return denied
		if (request.mimetype or "").lower() != "application/json":
			return _json_error(flask, 415, "Content-Type must be application/json; "
			                               f"got '{request.mimetype or 'nothing'}'.")
		try:
			payload = json.loads(request.get_data(as_text=True) or "")
		except json.JSONDecodeError as error:
			return _json_error(flask, 400, f"The body is not valid JSON: {error}")

		try:
			revision = registry.apply(plot_id, payload)
		except EnvelopeError as error:
			return _json_error(flask, 400, str(error))
		except ValueError as error:
			return _json_error(flask, 409, str(error))
		except KeyError as error:
			return _json_error(flask, 404, str(error).strip("'"))
		return flask.jsonify({"revision": revision})

	@blueprint.route("/plot/<plot_id>", methods=["DELETE"])
	def delete_plot(plot_id):
		denied = authorize()
		if denied is not None:
			return denied
		try:
			registry.delete_plot(plot_id)
		except KeyError as error:
			return _json_error(flask, 404, str(error).strip("'"))
		return flask.jsonify({"deleted": plot_id})

	@blueprint.route("/plot/<plot_id>.json", methods=["GET"])
	def plot_document(plot_id):
		try:
			document = registry.document(plot_id)
		except KeyError:
			return _json_error(flask, 404, f"No plot '{plot_id}'.")
		response = Response(canonical_json(document), mimetype="application/json")
		response.headers["ETag"] = f'"{envelope_hash(document)}"'
		return response

	@blueprint.route("/plot/<plot_id>.svg", methods=["GET"])
	def plot_svg(plot_id):
		width = 640
		raw = request.args.get("width")
		if raw is not None:
			try:
				width = int(raw)
			except ValueError:
				return _json_error(flask, 400, f"width must be an integer; got '{raw}'.")
			if not MINIMUM_SVG_WIDTH <= width <= MAXIMUM_SVG_WIDTH:
				return _json_error(flask, 400, f"width must lie in [{MINIMUM_SVG_WIDTH}, "
				                               f"{MAXIMUM_SVG_WIDTH}]; got {width}.")
		try:
			document = registry.document(plot_id)
		except KeyError:
			return _json_error(flask, 404, f"No plot '{plot_id}'.")
		response = Response(render_svg(document, width=width), mimetype="image/svg+xml")
		response.headers["ETag"] = f'"{envelope_hash(document)}-w{width}"'
		return response

	@blueprint.route("/events/<plot_id>", methods=["GET"])
	def events(plot_id):
		''' SSE. See caveat 1: this holds a worker for the connection's life. '''
		import queue as queue_module
		import time

		subscription = registry.subscribe(plot_id)

		def stream():
			try:
				yield f"retry: 2000\n\n"
				last_ping = time.monotonic()
				while True:
					if subscription.overflowed:
						return          # drop; the viewer reconnects and re-snapshots
					try:
						event, payload = subscription.queue.get(timeout=1.0)
					except queue_module.Empty:
						if time.monotonic() - last_ping >= 15.0:
							yield ": ping\n\n"
							last_ping = time.monotonic()
						continue
					serialized = canonical_json(payload)
					lines = "".join(f"data: {line}\n" for line in serialized.split("\n"))
					yield f"event: {event}\n{lines}\n"
					last_ping = time.monotonic()
			finally:
				registry.unsubscribe(subscription)

		return Response(stream(), mimetype="text/event-stream",
		                headers={"Cache-Control": "no-cache", "Connection": "keep-alive"})

	blueprint.registry = registry
	return blueprint


def _json_error(flask, status: int, message: str):
	response = flask.jsonify({"error": message})
	response.status_code = status
	return response


def _static_response(flask, name: str, mimetype):
	path = os.path.normpath(os.path.join(STATIC_DIRECTORY, name))
	if not path.startswith(STATIC_DIRECTORY) or not os.path.isfile(path):
		return _json_error(flask, 404, f"No static asset '{name}'.")
	if mimetype is None:
		mimetype = {"js": "text/javascript", "mjs": "text/javascript",
		            "html": "text/html", "json": "application/json"}.get(
			name.rsplit(".", 1)[-1], "application/octet-stream")
	with open(path, "rb") as stream:
		body = stream.read()
	return flask.Response(body, mimetype=mimetype, headers={"Cache-Control": "no-cache"})
