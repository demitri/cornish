
'''
The registry, the daemon and its SSE stream (SPEC-09A §7.5, plus the §7.2
registry state machine).

The daemon is exercised over real HTTP with the stdlib client on an ephemeral
port — not through a mock — because the frozen failure-mode precedence
(auth → Content-Type → parse/validate → plot_id mismatch → state) is a property
of the wire behaviour, not of the handler's internals.
'''

import http.client
import json
import os
import socket
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ElementTree

import pytest

from cornish.exc import EnvelopeError
from cornish.plot.envelope import (SCHEMA_VERSION, canonical_json, envelope_hash,
                                   make_envelope, LayerAction)
from cornish.view.registry import SUBSCRIBER_QUEUE_LIMIT, PlotRegistry
from cornish.view.server import serve

TOKEN = "test-token-not-a-secret"
CANVAS = {"gbox": [0, 0, 1, 1], "aspect": 1.0, "margin": 0.15, "char_height": 0.02, "width": 640}
LINE = {"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}


def envelope(**overrides) -> dict:
	base = {
		"cornish_view": SCHEMA_VERSION, "type": "draw_stream",
		"canvas": dict(CANVAS),
		"layer": {"id": "grid", "action": "replace", "kind": "grid", "z": -10},
		"commands": [LINE],
	}
	base.update(overrides)
	return base


def document(**overrides) -> dict:
	base = {
		"cornish_view": SCHEMA_VERSION, "type": "draw_stream",
		"canvas": dict(CANVAS),
		"layers": [{"id": "grid", "kind": "grid", "z": -10, "commands": [LINE]}],
	}
	base.update(overrides)
	return base


# ---------------------------------------------------------------------------
# registry state machine (§3.4, §7.2)
# ---------------------------------------------------------------------------

def test_create_by_replace():
	registry = PlotRegistry()
	assert registry.apply("p", envelope()) == 1
	assert registry.plot_ids() == ["p"]
	assert [layer["id"] for layer in registry.document("p")["layers"]] == ["grid"]


def test_a_plot_exists_only_once_it_has_a_canvas():
	'''
	§3.1: the creating envelope must carry a canvas, so no code path can ever
	observe a canvas-less document.
	'''
	registry = PlotRegistry()
	body = envelope()
	del body["canvas"]
	with pytest.raises(KeyError):
		registry.apply("p", body)
	assert registry.plot_ids() == []


def test_append_extends():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	registry.apply("p", envelope(layer={"id": "grid", "action": "append"},
	                             commands=[LINE], canvas=None) | {"canvas": dict(CANVAS)})
	assert len(registry.document("p")["layers"][0]["commands"]) == 2


def test_append_to_a_missing_layer_raises():
	'''
	Never a silent auto-create: an append is a statement about existing
	content, and creating the layer would mask a lost earlier envelope.
	'''
	registry = PlotRegistry()
	registry.apply("p", envelope())
	with pytest.raises(KeyError):
		registry.apply("p", envelope(layer={"id": "ghost", "action": "append"}, commands=[LINE]))


def test_delete_of_a_missing_layer_raises():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	body = envelope(layer={"id": "ghost", "action": "delete"})
	del body["commands"]
	with pytest.raises(KeyError):
		registry.apply("p", body)


def test_delete_removes_the_layer():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	body = envelope(layer={"id": "grid", "action": "delete"})
	del body["commands"]
	registry.apply("p", body)
	assert registry.document("p")["layers"] == []


def test_append_with_mismatched_z_or_kind_is_invalid():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	with pytest.raises(EnvelopeError, match="stored z"):
		registry.apply("p", envelope(layer={"id": "grid", "action": "append", "z": 5},
		                             commands=[LINE]))
	with pytest.raises(EnvelopeError, match="stored kind"):
		registry.apply("p", envelope(layer={"id": "grid", "action": "append", "kind": "points"},
		                             commands=[LINE]))


def test_append_with_matching_z_and_kind_is_fine():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	registry.apply("p", envelope(layer={"id": "grid", "action": "append", "z": -10, "kind": "grid"},
	                             commands=[LINE]))
	assert len(registry.document("p")["layers"][0]["commands"]) == 2


def test_replace_preserves_creation_order_among_z_ties():
	''' §3.4: a restyled layer keeps its slot. '''
	registry = PlotRegistry()
	registry.apply("p", envelope(layer={"id": "a", "action": "replace"}))
	registry.apply("p", envelope(layer={"id": "b", "action": "replace"}))
	registry.apply("p", envelope(layer={"id": "a", "action": "replace"}, commands=[]))
	assert [layer["id"] for layer in registry.document("p")["layers"]] == ["a", "b"]


def test_layers_are_ordered_by_z_then_creation():
	registry = PlotRegistry()
	registry.apply("p", envelope(layer={"id": "mid", "action": "replace"}))
	registry.apply("p", envelope(layer={"id": "back", "action": "replace", "z": -20},
	                             commands=[]))
	registry.apply("p", envelope(layer={"id": "front", "action": "replace", "z": 5}, commands=[]))
	assert [layer["id"] for layer in registry.document("p")["layers"]] == ["back", "mid", "front"]


def test_canvas_only_envelope_creates_and_rethemes():
	registry = PlotRegistry()
	registry.apply("p", make_envelope(canvas=dict(CANVAS)))
	assert registry.plot_ids() == ["p"]
	retheme = dict(CANVAS)
	retheme["background"] = "#ffffff"
	registry.apply("p", make_envelope(canvas=retheme))
	assert registry.document("p")["canvas"]["background"] == "#ffffff"


def test_document_apply_replaces_atomically():
	registry = PlotRegistry()
	registry.apply("p", envelope(layer={"id": "old", "action": "replace"}))
	registry.apply("p", document())
	assert [layer["id"] for layer in registry.document("p")["layers"]] == ["grid"]


def test_revision_is_monotonic():
	registry = PlotRegistry()
	revisions = [registry.apply("p", envelope()) for _ in range(4)]
	assert revisions == [1, 2, 3, 4]


def test_unknown_plot_reads_raise_key_error():
	registry = PlotRegistry()
	with pytest.raises(KeyError):
		registry.document("nope")
	with pytest.raises(KeyError):
		registry.svg("nope")
	with pytest.raises(KeyError):
		registry.delete_plot("nope")


def test_body_plot_id_mismatch_raises_value_error():
	registry = PlotRegistry()
	with pytest.raises(ValueError, match="does not match"):
		registry.apply("p", envelope(plot_id="other"))


def test_illegal_plot_id_raises():
	registry = PlotRegistry()
	with pytest.raises(ValueError, match="legal plot id"):
		registry.apply("not a plot id", envelope())


# ---------------------------------------------------------------------------
# subscriptions
# ---------------------------------------------------------------------------

def test_subscribing_to_an_existing_plot_gets_an_immediate_snapshot():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	subscription = registry.subscribe("p")
	event, payload = subscription.queue.get_nowait()
	assert event == "snapshot"
	assert payload["layers"][0]["id"] == "grid"


def test_subscribing_to_a_missing_plot_creates_a_slot_not_a_plot():
	'''
	§ registry: the id does not appear in plot_ids(), document() still raises,
	and the snapshot is deferred — a viewer opened before the first POST is a
	supported workflow.
	'''
	registry = PlotRegistry()
	subscription = registry.subscribe("later")
	assert registry.plot_ids() == []
	with pytest.raises(KeyError):
		registry.document("later")
	assert subscription.queue.empty()
	assert subscription.snapshot_pending


def test_deferred_snapshot_arrives_as_a_snapshot_not_an_envelope():
	'''
	The frozen §5.3 delivery exception: a subscriber still awaiting its deferred
	snapshot receives the plot-creating envelope AS a snapshot built from the
	just-committed document. The viewer's snapshot-then-deltas contract must
	hold for every subscriber regardless of subscription timing.
	'''
	registry = PlotRegistry()
	subscription = registry.subscribe("later")
	registry.apply("later", envelope())

	event, payload = subscription.queue.get_nowait()
	assert event == "snapshot"
	assert "layers" in payload and "layer" not in payload

	registry.apply("later", envelope(layer={"id": "second", "action": "replace"}, commands=[]))
	event, payload = subscription.queue.get_nowait()
	assert event == "envelope"


def test_document_apply_fans_out_as_a_snapshot():
	''' Deltas cannot express "everything you knew is gone". '''
	registry = PlotRegistry()
	registry.apply("p", envelope())
	subscription = registry.subscribe("p")
	subscription.queue.get_nowait()                     # the subscribe snapshot
	registry.apply("p", document())
	event, _ = subscription.queue.get_nowait()
	assert event == "snapshot"


def test_delete_plot_notifies_subscribers():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	subscription = registry.subscribe("p")
	subscription.queue.get_nowait()
	registry.delete_plot("p")
	event, _ = subscription.queue.get_nowait()
	assert event == "plot_deleted"


def test_unsubscribe_stops_delivery():
	registry = PlotRegistry()
	registry.apply("p", envelope())
	subscription = registry.subscribe("p")
	registry.unsubscribe(subscription)
	registry.unsubscribe(subscription)                  # idempotent
	assert registry.subscriber_count("p") == 0


def test_overflowing_subscriber_is_marked_rather_than_losing_events():
	'''
	A lagging viewer's connection is dropped so it reconnects and takes a clean
	snapshot — strictly better than silently discarding events and leaving it
	subtly wrong forever.
	'''
	registry = PlotRegistry()
	registry.apply("p", envelope())
	subscription = registry.subscribe("p")
	for _ in range(SUBSCRIBER_QUEUE_LIMIT + 10):
		registry.apply("p", envelope(layer={"id": "grid", "action": "replace"}, commands=[]))
	assert subscription.overflowed


def test_fan_out_happens_after_the_commit():
	''' A slow viewer can lag, never corrupt. '''
	registry = PlotRegistry()
	subscription = registry.subscribe("p")
	registry.apply("p", envelope())
	_, payload = subscription.queue.get_nowait()
	assert payload["layers"][0]["commands"] == [LINE]


# ---------------------------------------------------------------------------
# the daemon over real HTTP
# ---------------------------------------------------------------------------

@pytest.fixture
def daemon():
	server = serve(PlotRegistry(), port=0, token=TOKEN)
	try:
		yield server
	finally:
		server.shutdown_now()


def request(server, method, path, body=None, *, token=TOKEN,
            content_type="application/json", headers=None):
	''' One HTTP round trip against the daemon. '''
	connection = http.client.HTTPConnection("127.0.0.1", server.port, timeout=5)
	all_headers = dict(headers or {})
	if content_type is not None:
		all_headers["Content-Type"] = content_type
	if token is not None:
		all_headers["Authorization"] = f"Bearer {token}"
	payload = None if body is None else json.dumps(body).encode("utf-8")
	connection.request(method, path, body=payload, headers=all_headers)
	response = connection.getresponse()
	data = response.read()
	connection.close()
	return response.status, dict(response.getheaders()), data


def test_healthz(daemon):
	status, _, body = request(daemon, "GET", "/healthz", token=None, content_type=None)
	assert status == 200
	assert json.loads(body) == {"cornish_view": SCHEMA_VERSION, "plots": 0}


def test_post_envelope_returns_revision(daemon):
	status, _, body = request(daemon, "POST", "/plot/demo", envelope())
	assert status == 200
	assert json.loads(body) == {"revision": 1}


def test_post_document_replaces_atomically(daemon):
	request(daemon, "POST", "/plot/demo", envelope(layer={"id": "old", "action": "replace"}))
	status, _, _ = request(daemon, "POST", "/plot/demo", document())
	assert status == 200
	_, _, body = request(daemon, "GET", "/plot/demo.json", token=None, content_type=None)
	assert [layer["id"] for layer in json.loads(body)["layers"]] == ["grid"]


def test_get_json_matches_the_document_and_carries_the_hash_etag(daemon):
	request(daemon, "POST", "/plot/demo", envelope())
	status, headers, body = request(daemon, "GET", "/plot/demo.json", token=None, content_type=None)
	assert status == 200
	stored = json.loads(body)
	assert body.decode("utf-8") == canonical_json(stored)
	assert headers["ETag"] == f'"{envelope_hash(stored)}"'


def test_get_svg_parses_and_etag_tracks_width(daemon):
	request(daemon, "POST", "/plot/demo", envelope())
	status, headers, body = request(daemon, "GET", "/plot/demo.svg", token=None, content_type=None)
	assert status == 200
	ElementTree.fromstring(body.decode("utf-8"))
	assert headers["Content-Type"].startswith("image/svg+xml")
	assert headers["ETag"].endswith('-w640"')

	# width changes the bytes, so it must change the ETag
	_, wide_headers, wide_body = request(daemon, "GET", "/plot/demo.svg?width=1280",
	                                     token=None, content_type=None)
	assert wide_headers["ETag"] != headers["ETag"]
	assert wide_body != body


@pytest.mark.parametrize("query", ["width=10", "width=9999", "width=abc", "width=640.5"])
def test_svg_width_out_of_range_is_400_never_silently_clamped(daemon, query):
	request(daemon, "POST", "/plot/demo", envelope())
	status, _, _ = request(daemon, "GET", f"/plot/demo.svg?{query}", token=None, content_type=None)
	assert status == 400


def test_post_invalid_envelope_is_400_with_the_codec_message(daemon):
	broken = envelope(commands=[{"op": "line", "x": [0.0], "y": [0.0]}])
	status, _, body = request(daemon, "POST", "/plot/demo", broken)
	assert status == 400
	assert "at least 2 point" in json.loads(body)["error"]


@pytest.mark.parametrize("factory, fragment", [
	(lambda: {key: value for key, value in document().items() if key != "canvas"}, "canvas"),
	(lambda: document(layers=[{"id": "a", "commands": []}, {"id": "a", "commands": []}]), "duplicate"),
	(lambda: dict(document(), layer={"id": "grid", "action": "replace"}), "'layer' key"),
])
def test_post_invalid_document_is_400(daemon, factory, fragment):
	status, _, body = request(daemon, "POST", "/plot/demo", factory())
	assert status == 400
	assert fragment in json.loads(body)["error"]


def test_missing_token_is_401_and_wrong_token_is_403(daemon):
	status, _, _ = request(daemon, "POST", "/plot/demo", envelope(), token=None)
	assert status == 401
	status, _, _ = request(daemon, "POST", "/plot/demo", envelope(), token="wrong")
	assert status == 403


def test_auth_precedes_body_handling(daemon):
	'''
	Frozen precedence: an unauthenticated caller must learn nothing about how
	bodies are handled, so auth answers even for a body that is also invalid.
	'''
	status, _, _ = request(daemon, "POST", "/plot/demo", {"garbage": True},
	                       token=None, content_type="text/plain")
	assert status == 401


def test_wrong_content_type_is_415(daemon):
	status, _, _ = request(daemon, "POST", "/plot/demo", envelope(), content_type="text/plain")
	assert status == 415


def test_content_type_precedes_body_validation(daemon):
	status, _, _ = request(daemon, "POST", "/plot/demo", {"nonsense": 1}, content_type="text/plain")
	assert status == 415


def test_unparseable_body_is_400(daemon):
	connection = http.client.HTTPConnection("127.0.0.1", daemon.port, timeout=5)
	connection.request("POST", "/plot/demo", body=b"{not json",
	                   headers={"Content-Type": "application/json",
	                            "Authorization": f"Bearer {TOKEN}"})
	response = connection.getresponse()
	status = response.status
	response.read()
	connection.close()
	assert status == 400


def test_plot_id_mismatch_is_409(daemon):
	status, _, _ = request(daemon, "POST", "/plot/demo", envelope(plot_id="somethingelse"))
	assert status == 409


def test_validation_precedes_the_mismatch_check(daemon):
	''' The frozen order is parse/validate (400) THEN mismatch (409). '''
	broken = envelope(plot_id="somethingelse",
	                  commands=[{"op": "line", "x": [0.0], "y": [0.0]}])
	status, _, _ = request(daemon, "POST", "/plot/demo", broken)
	assert status == 400


def test_append_to_a_missing_layer_is_404(daemon):
	request(daemon, "POST", "/plot/demo", envelope())
	status, _, _ = request(daemon, "POST", "/plot/demo",
	                       envelope(layer={"id": "ghost", "action": "append"}, commands=[LINE]))
	assert status == 404


def test_envelope_for_an_unknown_plot_is_404(daemon):
	body = envelope()
	del body["canvas"]
	status, _, _ = request(daemon, "POST", "/plot/nothere", body)
	assert status == 404


def test_delete_then_get_is_404(daemon):
	request(daemon, "POST", "/plot/demo", envelope())
	status, _, _ = request(daemon, "DELETE", "/plot/demo", content_type=None)
	assert status == 200
	status, _, _ = request(daemon, "GET", "/plot/demo.json", token=None, content_type=None)
	assert status == 404
	status, _, _ = request(daemon, "DELETE", "/plot/demo", content_type=None)
	assert status == 404


def test_get_unknown_plot_is_404(daemon):
	for path in ("/plot/nope.json", "/plot/nope.svg"):
		status, _, _ = request(daemon, "GET", path, token=None, content_type=None)
		assert status == 404


def test_index_and_viewer_and_static_are_served(daemon):
	status, headers, body = request(daemon, "GET", "/", token=None, content_type=None)
	assert status == 200 and b"cornish-view" in body

	# an unknown id still returns 200: the page subscribes and waits
	status, _, body = request(daemon, "GET", "/view/not-yet", token=None, content_type=None)
	assert status == 200 and b"CornishView" in body

	status, headers, body = request(daemon, "GET", "/static/cornish-view.js",
	                                token=None, content_type=None)
	assert status == 200
	assert headers["Content-Type"].startswith("text/javascript")
	assert headers["Cache-Control"] == "no-cache"
	assert b"export class CornishView" in body


def test_static_path_traversal_is_refused(daemon):
	status, _, _ = request(daemon, "GET", "/static/..%2f..%2fsetup.py",
	                       token=None, content_type=None)
	assert status == 404


# ---------------------------------------------------------------------------
# SSE framing (§5.3)
# ---------------------------------------------------------------------------

def read_sse(server, plot_id, *, expected_events, timeout=10.0):
	'''
	Open ``/events/{plot_id}`` on a raw socket and return the text read until
	*expected_events* ``event:`` lines have arrived (or the timeout expires).
	'''
	connection = socket.create_connection(("127.0.0.1", server.port), timeout=timeout)
	connection.sendall(f"GET /events/{plot_id} HTTP/1.1\r\nHost: localhost\r\n\r\n".encode("utf-8"))
	connection.settimeout(0.5)
	text = ""
	deadline = time.monotonic() + timeout
	while time.monotonic() < deadline:
		if text.count("\nevent: ") >= expected_events:
			break
		try:
			chunk = connection.recv(65536)
		except socket.timeout:
			continue
		if not chunk:
			break
		text += chunk.decode("utf-8", "replace")
	connection.close()
	return text


def test_sse_opens_with_retry_then_a_snapshot(daemon):
	request(daemon, "POST", "/plot/demo", envelope())
	text = read_sse(daemon, "demo", expected_events=1)
	assert "text/event-stream" in text
	assert "retry: 2000" in text
	assert "\nevent: snapshot\n" in text
	payload = json.loads(text.split("event: snapshot\ndata: ", 1)[1].split("\n", 1)[0])
	assert payload["type"] == "draw_stream"
	assert payload["layers"][0]["id"] == "grid"


def test_sse_delivers_a_later_post_as_an_envelope(daemon):
	request(daemon, "POST", "/plot/demo", envelope())

	def post_later():
		time.sleep(0.4)
		request(daemon, "POST", "/plot/demo",
		        envelope(layer={"id": "second", "action": "replace"}, commands=[LINE]))

	thread = threading.Thread(target=post_later)
	thread.start()
	text = read_sse(daemon, "demo", expected_events=2)
	thread.join()

	assert "\nevent: snapshot\n" in text
	assert "\nevent: envelope\n" in text
	assert text.index("event: snapshot") < text.index("event: envelope")


def test_sse_deferred_snapshot_for_a_viewer_opened_before_the_first_post(daemon):
	'''
	The §5.3 pin: subscribe to a not-yet-POSTed id, then create the plot — the
	FIRST event must be a snapshot carrying the full document, never an
	envelope. A further POST then arrives as an envelope.
	'''
	def post_later():
		time.sleep(0.4)
		request(daemon, "POST", "/plot/deferred", envelope())
		time.sleep(0.3)
		request(daemon, "POST", "/plot/deferred",
		        envelope(layer={"id": "second", "action": "replace"}, commands=[LINE]))

	thread = threading.Thread(target=post_later)
	thread.start()
	text = read_sse(daemon, "deferred", expected_events=2)
	thread.join()

	first = text.index("event: ")
	assert text[first:].startswith("event: snapshot")
	snapshot = json.loads(text.split("event: snapshot\ndata: ", 1)[1].split("\n", 1)[0])
	assert "layers" in snapshot and "layer" not in snapshot
	assert "\nevent: envelope\n" in text


def test_sse_data_is_a_single_line(daemon):
	'''
	Canonical JSON is single-line, so each event is exactly one data: line. The
	writer still splits defensively on newlines — it simply never triggers.
	'''
	request(daemon, "POST", "/plot/demo", envelope())
	text = read_sse(daemon, "demo", expected_events=1)
	block = text.split("event: snapshot\n", 1)[1].split("\n\n", 1)[0]
	assert block.count("data: ") == 1


def test_sse_unsubscribes_when_the_viewer_disconnects(daemon):
	'''
	The stream notices a dead viewer on its next write and drops the
	subscription. A POST is what forces that write here — otherwise the next
	write would be the 15-second keep-alive ping, and the test would just be
	measuring the ping interval.
	'''
	request(daemon, "POST", "/plot/demo", envelope())
	read_sse(daemon, "demo", expected_events=1)          # opens and closes the socket
	assert daemon.registry.subscriber_count("demo") == 1

	request(daemon, "POST", "/plot/demo",
	        envelope(layer={"id": "second", "action": "replace"}, commands=[LINE]))
	deadline = time.monotonic() + 5
	while time.monotonic() < deadline and daemon.registry.subscriber_count("demo") != 0:
		time.sleep(0.05)
	assert daemon.registry.subscriber_count("demo") == 0


# ---------------------------------------------------------------------------
# token file and ensure_daemon
# ---------------------------------------------------------------------------

def test_token_file_is_created_private(tmp_path):
	from cornish.view.server import read_or_create_token
	path = tmp_path / "nested" / "view_token"
	token = read_or_create_token(str(path))
	assert len(token) > 20
	assert read_or_create_token(str(path)) == token          # stable across reads
	assert oct(os.stat(path).st_mode)[-3:] == "600"
	assert oct(os.stat(path.parent).st_mode)[-3:] == "700"


@pytest.mark.slow
def test_ensure_daemon_spawns_and_answers():
	'''
	Spawn a real daemon in a subprocess, probe healthz, and shut it down. Marked
	slow: it starts an interpreter.
	'''
	from cornish.view.server import read_or_create_token

	port = _free_port()
	source_directory = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
	process = subprocess.Popen([sys.executable, "-m", "cornish.view", "--port", str(port)],
	                           cwd=source_directory,
	                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	try:
		deadline = time.monotonic() + 20
		healthy = False
		while time.monotonic() < deadline and not healthy:
			try:
				connection = http.client.HTTPConnection("127.0.0.1", port, timeout=1)
				connection.request("GET", "/healthz")
				healthy = connection.getresponse().status == 200
				connection.close()
			except OSError:
				time.sleep(0.2)
		assert healthy, "the spawned daemon never answered /healthz"

		from cornish.view.server import ensure_daemon
		base_url, token = ensure_daemon(port=port, launch=False)
		assert base_url == f"http://127.0.0.1:{port}"
		assert token == read_or_create_token()
	finally:
		process.terminate()
		process.wait(timeout=10)


def test_ensure_daemon_without_launch_raises_when_nothing_is_listening():
	from cornish.view.server import ensure_daemon
	with pytest.raises(RuntimeError, match="No cornish-view daemon"):
		ensure_daemon(port=_free_port(), launch=False)


def _free_port() -> int:
	with socket.socket() as probe:
		probe.bind(("127.0.0.1", 0))
		return probe.getsockname()[1]


# ---------------------------------------------------------------------------
# SkyPlot live mode (§2.4 show)
# ---------------------------------------------------------------------------

def test_skyplot_live_mode_posts_each_layer(daemon, monkeypatch):
	'''
	``show(view='browser')`` bootstraps with the full document, then every
	subsequent draw method POSTs its own single-layer envelope.
	'''
	import astropy.units as u

	from cornish import ASTCircle, ASTICRSFrame
	from cornish.plot.sky import SkyPlot
	import cornish.view.server as server_module

	monkeypatch.setattr(server_module, "ensure_daemon",
	                    lambda **kwargs: (daemon.base_url, TOKEN))
	opened = []
	monkeypatch.setattr("webbrowser.open", lambda url: opened.append(url))

	plot = SkyPlot(extent=ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg),
	               plot_id="live")
	plot.show(view="browser")

	assert plot.live
	assert opened == [f"{daemon.base_url}/view/live"]
	assert daemon.registry.plot_ids() == ["live"]

	before = daemon.registry.revision("live")
	plot.addPoints(ra=[30.0], dec=[45.0])
	assert daemon.registry.revision("live") > before
	# the grid layer has a FIXED id, so it does not consume a counter slot
	assert "points-1" in [layer["id"] for layer in daemon.registry.document("live")["layers"]]


def test_skyplot_live_mode_removal_propagates(daemon, monkeypatch):
	''' A viewer must not keep showing a layer the plot no longer has. '''
	import astropy.units as u

	from cornish import ASTCircle, ASTICRSFrame
	from cornish.plot.sky import SkyPlot
	import cornish.view.server as server_module

	monkeypatch.setattr(server_module, "ensure_daemon",
	                    lambda **kwargs: (daemon.base_url, TOKEN))
	monkeypatch.setattr("webbrowser.open", lambda url: None)

	plot = SkyPlot(extent=ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg),
	               plot_id="live2")
	plot.show(view="browser")
	plot.removeLayer("grid")
	assert daemon.registry.document("live2")["layers"] == []


def test_skyplot_live_post_failure_surfaces(monkeypatch):
	'''
	Live posts are synchronous precisely so a failure raises. A fire-and-forget
	design would lose the LAST draw call of a session silently — the one plot
	state the user most cares about.
	'''
	import astropy.units as u

	from cornish import ASTCircle, ASTICRSFrame
	from cornish.plot.sky import SkyPlot
	import cornish.view.server as server_module

	dead_url = f"http://127.0.0.1:{_free_port()}"
	monkeypatch.setattr(server_module, "ensure_daemon", lambda **kwargs: (dead_url, TOKEN))
	monkeypatch.setattr("webbrowser.open", lambda url: None)

	plot = SkyPlot(extent=ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg),
	               plot_id="dead")
	with pytest.raises(RuntimeError, match="Could not reach"):
		plot.show(view="browser")
	assert not plot.live, "a failed bootstrap must not leave live mode half-armed"


def test_show_rejects_an_unknown_view():
	import astropy.units as u

	from cornish import ASTCircle, ASTICRSFrame
	from cornish.plot.sky import SkyPlot

	plot = SkyPlot(extent=ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg))
	with pytest.raises(ValueError, match="not a known view"):
		plot.show(view="hologram")
