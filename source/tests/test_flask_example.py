
'''
The Flask pattern-1 example, driven with Flask's test client (SPEC-09A §7.5).

The example file *is* the documentation, so it is imported and executed rather
than merely referenced — documentation that is never run is documentation that
silently rots. Also covers the blueprint's own routes, since it shares its
handler behaviour with the daemon.

Skipped when Flask is not installed; Flask is not a cornish dependency.
'''

import importlib.util
import os
import xml.etree.ElementTree as ElementTree

import pytest

pytest.importorskip("flask")

EXAMPLE_PATH = os.path.normpath(os.path.join(
	os.path.dirname(__file__), "..", "..", "docs", "examples", "flask_footprint.py"))


@pytest.fixture(scope="module")
def example():
	''' Import the example file by path — it lives in docs/, not on sys.path. '''
	assert os.path.exists(EXAMPLE_PATH), f"the documented example is missing: {EXAMPLE_PATH}"
	spec = importlib.util.spec_from_file_location("flask_footprint", EXAMPLE_PATH)
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	return module


@pytest.fixture
def client(example):
	example.app.config["TESTING"] = True
	return example.app.test_client()


def test_footprint_renders_an_svg(client):
	response = client.get("/footprint.svg?ra=30&dec=45&radius=2")
	assert response.status_code == 200
	assert response.mimetype == "image/svg+xml"
	ElementTree.fromstring(response.get_data(as_text=True))


def test_etag_is_stable_and_supports_conditional_requests(client):
	'''
	render_svg is deterministic, so equal documents give equal bytes and equal
	ETags — which is what makes the 304 correct rather than merely plausible.
	'''
	first = client.get("/footprint.svg?ra=30&dec=45&radius=2")
	second = client.get("/footprint.svg?ra=30&dec=45&radius=2")
	assert first.headers["ETag"] == second.headers["ETag"]
	assert first.get_data() == second.get_data()

	conditional = client.get("/footprint.svg?ra=30&dec=45&radius=2",
	                         headers={"If-None-Match": first.headers["ETag"]})
	assert conditional.status_code == 304


def test_width_changes_the_etag(client):
	''' Width changes the bytes, so it must change the ETag. '''
	narrow = client.get("/footprint.svg?ra=30&dec=45&radius=2&width=640")
	wide = client.get("/footprint.svg?ra=30&dec=45&radius=2&width=1280")
	assert narrow.headers["ETag"] != wide.headers["ETag"]
	assert narrow.get_data() != wide.get_data()


@pytest.mark.parametrize("query", [
	"ra=abc&dec=45&radius=2", "ra=30&dec=45&radius=0", "ra=30&dec=45&radius=1000",
	"ra=30&dec=45&radius=2&width=1", "ra=30&dec=45&radius=2&width=99999",
])
def test_bad_parameters_are_400_never_silently_defaulted(client, query):
	assert client.get(f"/footprint.svg?{query}").status_code == 400


def test_index_embeds_the_svg(client):
	response = client.get("/")
	assert response.status_code == 200
	assert b"/footprint.svg" in response.get_data()


def test_build_document_is_a_valid_draw_stream(example):
	from cornish.plot.envelope import validate_document
	document = example.build_document(30.0, 45.0, 2.0)
	validate_document(document)
	assert "footprint" in [layer["id"] for layer in document["layers"]]


# ---------------------------------------------------------------------------
# the blueprint (§5.2, §5.6)
# ---------------------------------------------------------------------------

@pytest.fixture
def blueprint_client():
	import flask

	from cornish.view import PlotRegistry, flask_blueprint

	registry = PlotRegistry()
	app = flask.Flask(__name__)
	app.config["TESTING"] = True
	app.register_blueprint(flask_blueprint(registry, require_token=False), url_prefix="/plots")
	return app.test_client(), registry


CANVAS = {"gbox": [0, 0, 1, 1], "aspect": 1.0, "margin": 0.15, "char_height": 0.02}
ENVELOPE = {
	"cornish_view": "1.0", "type": "draw_stream", "canvas": CANVAS,
	"layer": {"id": "grid", "action": "replace", "kind": "grid", "z": -10},
	"commands": [{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}],
}


def test_blueprint_exposes_the_same_routes(blueprint_client):
	client, registry = blueprint_client

	assert client.get("/plots/healthz").get_json()["plots"] == 0

	response = client.post("/plots/plot/demo", json=ENVELOPE)
	assert response.status_code == 200
	assert response.get_json() == {"revision": 1}

	document = client.get("/plots/plot/demo.json")
	assert document.status_code == 200
	assert document.headers["ETag"]

	svg = client.get("/plots/plot/demo.svg")
	assert svg.status_code == 200
	ElementTree.fromstring(svg.get_data(as_text=True))

	assert client.delete("/plots/plot/demo").status_code == 200
	assert client.get("/plots/plot/demo.json").status_code == 404


def test_blueprint_shares_the_daemon_failure_precedence(blueprint_client):
	''' The blueprint and the daemon share handler behaviour, not just routes. '''
	client, _ = blueprint_client

	# 415 before body validation
	assert client.post("/plots/plot/demo", data="{}", content_type="text/plain").status_code == 415

	# 400 for an invalid body
	broken = dict(ENVELOPE, commands=[{"op": "line", "x": [0.0], "y": [0.0]}])
	assert client.post("/plots/plot/demo", json=broken).status_code == 400

	# 409 for a body/URL plot_id mismatch
	mismatched = dict(ENVELOPE, plot_id="elsewhere")
	assert client.post("/plots/plot/demo", json=mismatched).status_code == 409

	# 404 for an append to a layer that does not exist
	client.post("/plots/plot/demo", json=ENVELOPE)
	appended = dict(ENVELOPE, layer={"id": "ghost", "action": "append"})
	assert client.post("/plots/plot/demo", json=appended).status_code == 404


def test_blueprint_svg_width_is_validated(blueprint_client):
	client, _ = blueprint_client
	client.post("/plots/plot/demo", json=ENVELOPE)
	assert client.get("/plots/plot/demo.svg?width=10").status_code == 400
	assert client.get("/plots/plot/demo.svg?width=abc").status_code == 400
	assert client.get("/plots/plot/demo.svg?width=1280").status_code == 200


def test_blueprint_requires_a_token_when_asked(blueprint_client, monkeypatch):
	import flask

	from cornish.view import PlotRegistry
	import cornish.view.blueprint as blueprint_module

	monkeypatch.setattr(blueprint_module, "read_or_create_token", lambda: "a-token")
	app = flask.Flask(__name__)
	app.config["TESTING"] = True
	app.register_blueprint(blueprint_module.flask_blueprint(PlotRegistry(), require_token=True),
	                       url_prefix="/plots")
	client = app.test_client()

	assert client.post("/plots/plot/demo", json=ENVELOPE).status_code == 401
	assert client.post("/plots/plot/demo", json=ENVELOPE,
	                   headers={"Authorization": "Bearer wrong"}).status_code == 403
	assert client.post("/plots/plot/demo", json=ENVELOPE,
	                   headers={"Authorization": "Bearer a-token"}).status_code == 200
	# reads stay open by design — EventSource cannot send an Authorization header
	assert client.get("/plots/plot/demo.json").status_code == 200


def test_blueprint_documents_its_worker_model_caveats():
	'''
	The §5.6 caveats are load-bearing operational facts (SSE pins a worker; the
	registry is per-process; public mounting exposes every plot's data). They
	must live where a deployer will meet them.
	'''
	from cornish.view.blueprint import flask_blueprint

	docstring = flask_blueprint.__doc__
	for fragment in ["SSE pins a worker", "per-process", "unauthenticated",
	                 "gthread", "before_request"]:
		assert fragment in docstring, f"the blueprint docstring omits '{fragment}'"


def test_blueprint_import_error_is_actionable(monkeypatch):
	'''
	A missing Flask must produce guidance, not a bare ImportError traceback:
	Flask is not a cornish dependency and there are two ways to proceed without it.
	'''
	import builtins

	import cornish.view.blueprint as blueprint_module

	real_import = builtins.__import__

	def fake_import(name, *args, **kwargs):
		if name == "flask":
			raise ImportError("No module named 'flask'")
		return real_import(name, *args, **kwargs)

	monkeypatch.setattr(builtins, "__import__", fake_import)
	with pytest.raises(ImportError) as raised:
		blueprint_module.flask_blueprint()
	message = str(raised.value)
	assert "pip install flask" in message
	assert "render_svg" in message and "python -m cornish.view" in message
