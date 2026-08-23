
'''
The machine-readable JSON Schema, held in sync with the codec (SPEC-09A §1, §3.7).

``cornish/view/static/envelope.schema.json`` is **hand-maintained**, not
generated. That is a deliberate choice — but a hand-maintained mirror decays
unless something checks it, so this module is that something: it validates every
conformance fixture against the schema AND asserts schema/codec agreement on the
invalid-document battery, so the two cannot drift apart silently.

``jsonschema`` is a test-only dependency.
'''

import json
import os

import pytest

from cornish.exc import EnvelopeError
from cornish.plot.envelope import validate_document, validate_envelope

jsonschema = pytest.importorskip("jsonschema")

SCHEMA_PATH = os.path.join(os.path.dirname(__file__), "..", "cornish", "view", "static",
                           "envelope.schema.json")
FIXTURE_DIRECTORY = os.path.join(os.path.dirname(__file__), "fixtures", "draw_streams")

CANVAS = {"gbox": [0, 0, 1, 1], "aspect": 1.0, "margin": 0.15, "char_height": 0.02}


@pytest.fixture(scope="module")
def schema():
	with open(os.path.normpath(SCHEMA_PATH), encoding="utf-8") as stream:
		return json.load(stream)


@pytest.fixture(scope="module")
def validator(schema):
	cls = jsonschema.validators.validator_for(schema)
	cls.check_schema(schema)
	return cls(schema)


def schema_accepts(validator, document) -> bool:
	return validator.is_valid(document)


def codec_accepts(document) -> bool:
	''' Whichever validator owns the body's shape (§3.1 dispatch). '''
	try:
		if "layers" in document:
			validate_document(document)
		else:
			validate_envelope(document)
	except EnvelopeError:
		return False
	return True


def envelope(**overrides) -> dict:
	base = {
		"cornish_view": "1.0", "type": "draw_stream", "canvas": dict(CANVAS),
		"layer": {"id": "grid", "action": "replace"}, "commands": [],
	}
	base.update(overrides)
	return base


def document(**overrides) -> dict:
	base = {
		"cornish_view": "1.0", "type": "draw_stream", "canvas": dict(CANVAS),
		"layers": [{"id": "grid", "kind": "grid", "z": -10, "commands": []}],
	}
	base.update(overrides)
	return base


def test_the_schema_is_itself_valid(schema):
	cls = jsonschema.validators.validator_for(schema)
	cls.check_schema(schema)


def test_every_golden_fixture_validates(validator):
	names = sorted(name for name in os.listdir(FIXTURE_DIRECTORY) if name.endswith(".json"))
	assert names, "no golden fixtures found"
	for name in names:
		with open(os.path.join(FIXTURE_DIRECTORY, name), encoding="utf-8") as stream:
			fixture = json.load(stream)
		errors = sorted(validator.iter_errors(fixture), key=lambda error: error.path)
		assert not errors, f"{name} failed the schema: {[error.message for error in errors]}"


#: Bodies both the schema and the codec must ACCEPT.
VALID_BODIES = [
	("envelope", envelope()),
	("envelope with every op", envelope(commands=[
		{"op": "attr", "prim": "line", "colour": 11, "width": 0.5, "style": 2},
		{"op": "line", "x": [0.0, 1.0], "y": [0.0, 1.0]},
		{"op": "mark", "x": [0.5], "y": [0.5], "marker": 4},
		{"op": "text", "text": "M51", "x": 0.5, "y": 0.5, "just": "BC", "up": [0, 1]},
		{"op": "image", "x1": 0, "y1": 0, "x2": 1, "y2": 1, "href": "data:image/png;base64,AA=="},
	])),
	("canvas-only envelope", {"cornish_view": "1.0", "type": "draw_stream", "canvas": dict(CANVAS)}),
	("delete envelope", {"cornish_view": "1.0", "type": "draw_stream",
	                     "layer": {"id": "grid", "action": "delete"}}),
	("document", document()),
	("unknown op tolerated", envelope(commands=[{"op": "hologram", "intensity": 3}])),
	("unknown minor tolerated", envelope(cornish_view="1.7")),
]

#: Bodies both the schema and the codec must REJECT.
INVALID_BODIES = [
	("bad major", envelope(cornish_view="2.0")),
	("unknown type", envelope(type="plot_spec")),
	("commands on delete", {"cornish_view": "1.0", "type": "draw_stream",
	                        "layer": {"id": "grid", "action": "delete"}, "commands": []}),
	("commands on a canvas-only envelope", {"cornish_view": "1.0", "type": "draw_stream",
	                                        "canvas": dict(CANVAS), "commands": []}),
	("missing commands", {"cornish_view": "1.0", "type": "draw_stream", "canvas": dict(CANVAS),
	                      "layer": {"id": "grid", "action": "replace"}}),
	("bad layer id", envelope(layer={"id": "a b", "action": "replace"})),
	("bad action", envelope(layer={"id": "grid", "action": "paint"})),
	("bad plot_id", envelope(plot_id="not a plot id")),
	("envelope carrying layers", envelope(layers=[])),
	("document carrying layer", dict(document(), layer={"id": "grid", "action": "replace"})),
	("document without canvas", {"cornish_view": "1.0", "type": "draw_stream",
	                             "layers": [{"id": "a", "commands": []}]}),
	("layer without commands", document(layers=[{"id": "a"}])),
	("layer with an action", document(layers=[{"id": "a", "action": "replace", "commands": []}])),
	("line with one point", envelope(commands=[{"op": "line", "x": [0.0], "y": [0.0]}])),
	("mark without a marker", envelope(commands=[{"op": "mark", "x": [0.0], "y": [0.0]}])),
	("attr without prim", envelope(commands=[{"op": "attr", "colour": 1}])),
	("bad justification", envelope(commands=[{"op": "text", "text": "x", "x": 0, "y": 0,
	                                          "just": "ZZ", "up": [0, 1]}])),
	("remote image href", envelope(commands=[{"op": "image", "x1": 0, "y1": 0, "x2": 1, "y2": 1,
	                                          "href": "https://example.org/x.png"}])),
	("margin out of domain", envelope(canvas={"margin": 0.9})),
	("palette too short", envelope(canvas={"palette": ["#000000"]})),
	("neither layer nor layers nor canvas", {"cornish_view": "1.0", "type": "draw_stream"}),
]


@pytest.mark.parametrize("label, body", VALID_BODIES, ids=[case[0] for case in VALID_BODIES])
def test_schema_and_codec_agree_on_valid_bodies(validator, label, body):
	assert codec_accepts(body), f"the codec rejected {label}"
	errors = [error.message for error in validator.iter_errors(body)]
	assert not errors, f"the schema rejected {label}: {errors}"


@pytest.mark.parametrize("label, body", INVALID_BODIES, ids=[case[0] for case in INVALID_BODIES])
def test_schema_and_codec_agree_on_invalid_bodies(validator, label, body):
	'''
	The sync check that keeps a hand-maintained schema honest: for every case in
	the battery, the schema and the codec must reach the SAME verdict.
	'''
	assert not codec_accepts(body), f"the codec accepted {label}, which should be invalid"
	assert not schema_accepts(validator, body), f"the schema accepted {label}, which should be invalid"


def test_hybrid_body_is_rejected_by_both(validator):
	''' §1.1's ownership rule, expressed in the schema by `layer: false` / `layers: false`. '''
	hybrid = dict(document(), layer={"id": "grid", "action": "replace"})
	assert not codec_accepts(hybrid)
	assert not schema_accepts(validator, hybrid)


def test_schema_ships_next_to_the_renderer():
	''' It is part of the conformance fixture set astplot-js vendors (§7.6). '''
	assert os.path.exists(os.path.normpath(SCHEMA_PATH))
