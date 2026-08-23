
'''
Envelope schema v1.0 codec tests (SPEC-09A §7.2).

Three things are pinned here:

1. The two validators, including the §3.1 "exactly one of layer/layers" rule —
   which is enforced by each validator treating the *other* shape's marker key
   as a hard error, so a hybrid body fails whichever validator sees it first.
2. RFC 8785 (JCS) canonicalization, especially ``_jcs_number``. The number
   vectors are given as IEEE 754 bit patterns (the RFC's own appendix style)
   with the expected output taken from a real JavaScript engine, because the
   whole point of JCS here is that a hash pinned in Python reproduces in JS.
3. ``make_envelope``'s legal shapes.

The registry state machine (also listed under §7.2) is exercised in
``test_view_daemon.py``, where the registry lives.
'''

import json
import struct

import pytest

from cornish.exc import CornishError, EnvelopeError
from cornish.plot.envelope import (SCHEMA_VERSION, GrfOp, LayerAction, PrimType,
                                   canonical_json, envelope_hash, make_envelope,
                                   validate_document, validate_envelope)
from cornish.plot.envelope import _jcs_number


CANVAS = {"gbox": [0, 0, 1, 1], "aspect": 1.0, "margin": 0.15, "char_height": 0.02}


def envelope(**overrides) -> dict:
	''' A minimal valid envelope, with fields overridden/removed for the negative cases. '''
	base = {
		"cornish_view": SCHEMA_VERSION,
		"type": "draw_stream",
		"canvas": dict(CANVAS),
		"layer": {"id": "grid", "action": "replace"},
		"commands": [],
	}
	for key, value in overrides.items():
		if value is REMOVE:
			base.pop(key, None)
		else:
			base[key] = value
	return base


class _Remove:
	def __repr__(self):
		return "REMOVE"


REMOVE = _Remove()


# ---------------------------------------------------------------------------
# validate_envelope
# ---------------------------------------------------------------------------

def test_envelope_error_is_a_cornish_error():
	''' SPEC-10 policy: every cornish-specific exception descends from CornishError. '''
	assert issubclass(EnvelopeError, CornishError)


def test_validate_envelope_accepts_every_command_op():
	'''
	Acceptance of every §3.3 op in one envelope — the positive half of the
	codec contract.
	'''
	commands = [
		{"op": GrfOp.ATTR.value, "prim": PrimType.LINE.value, "colour": 11, "width": 0.5, "style": 2},
		{"op": GrfOp.ATTR.value, "prim": PrimType.TEXT.value, "colour": "#4a7f7b", "size": 1.2, "font": 3},
		{"op": GrfOp.LINE.value, "x": [0.0, 0.5, 1.0], "y": [0.0, 0.5, 1.0]},
		{"op": GrfOp.MARK.value, "x": [0.25], "y": [0.75], "marker": 4},
		{"op": GrfOp.TEXT.value, "text": "M51", "x": 0.5, "y": 0.5, "just": "BC", "up": [0.0, 1.0]},
		{"op": GrfOp.IMAGE.value, "x1": 0.0, "y1": 0.0, "x2": 1.0, "y2": 1.0,
		 "href": "data:image/png;base64,iVBORw0KGgo="},
	]
	validate_envelope(envelope(commands=commands))


def test_validate_envelope_accepts_canvas_only_envelope():
	''' §3.1†: the plot-creating / re-theming form carries a canvas and nothing else. '''
	validate_envelope(envelope(layer=REMOVE, commands=REMOVE))


def test_validate_envelope_accepts_delete_without_commands():
	validate_envelope(envelope(layer={"id": "grid", "action": "delete"}, commands=REMOVE))


def test_validate_envelope_accepts_layer_without_canvas():
	'''
	``canvas`` is required only on the envelope that CREATES a plot — a rule
	about registry state, not about shape. A later delta legitimately omits it.
	'''
	validate_envelope(envelope(canvas=REMOVE))


def test_validate_envelope_accepts_empty_commands():
	''' §3.1: an empty replace clears a layer but keeps its slot. '''
	validate_envelope(envelope(commands=[]))


@pytest.mark.parametrize("doc, fragment", [
	(envelope(cornish_view="2.0"), "cornish-view 2.x"),
	(envelope(cornish_view="one point oh"), "not a '<major>.<minor>'"),
	(envelope(type="plot_spec"), "draw_stream"),
	# no layer AND no canvas: there is no legal layer-less shape without a canvas
	(envelope(layer=REMOVE, canvas=REMOVE), "'layer' is missing"),
	(envelope(layer=REMOVE, canvas=REMOVE, commands=REMOVE), "'layer' is missing"),
	(envelope(commands=REMOVE), "'commands' is missing"),
	(envelope(layer={"id": "grid", "action": "delete"}), "must not carry 'commands'"),
	(envelope(layer=REMOVE, commands=[]), "canvas-only envelope must not carry 'commands'"),
	(envelope(layer={"id": "a b", "action": "replace"}), "not a legal id"),
	(envelope(layer={"id": "grid", "action": "paint"}), "expected one of"),
	(envelope(layer={"id": "grid"}), "'action' is missing"),
	(envelope(layer={"action": "replace"}), "'id' is missing"),
	(envelope(plot_id="not a plot id"), "not a legal id"),
	(envelope(layers=[]), "must not carry a 'layers' key"),
	(envelope(canvas={"gbox": [0, 0, 1, 2], "aspect": 1.0}), "contradicts the gbox"),
	(envelope(canvas={"margin": 0.9}), "[0, 0.5]"),
	(envelope(canvas={"char_height": 0}), "must be > 0"),
	(envelope(canvas={"width": 0}), "positive integer"),
	(envelope(canvas={"palette": ["#000000"]}), "at least 2 entries"),
	(envelope(canvas={"palette": ["#000000", "not-a-colour"]}), "not a valid palette entry"),
	(envelope(commands=[{"op": "line", "x": [0.0, float("nan")], "y": [0.0, 1.0]}]), "finite"),
	(envelope(commands=[{"op": "line", "x": [0.0, float("inf")], "y": [0.0, 1.0]}]), "finite"),
	(envelope(commands=[{"op": "line", "x": [0.0], "y": [0.0]}]), "at least 2 point"),
	(envelope(commands=[{"op": "line", "x": [0.0, 1.0], "y": [0.0]}]), "same length"),
	(envelope(commands=[{"op": "mark", "x": [0.0], "y": [0.0]}]), "needs a 'marker'"),
	(envelope(commands=[{"op": "attr", "colour": 1}]), "needs a 'prim'"),
	(envelope(commands=[{"op": "attr", "prim": "border", "colour": 1}]), "expected one of"),
	(envelope(commands=[{"op": "text", "text": "x", "x": 0, "y": 0, "just": "ZZ", "up": [0, 1]}]), "not a legal justification"),
	(envelope(commands=[{"op": "text", "text": "x", "x": 0, "y": 0, "just": "BC", "up": [0, 0]}]), "must be non-zero"),
	(envelope(commands=[{"op": "image", "x1": 0, "y1": 0, "x2": 1, "y2": 1,
	                     "href": "https://example.org/x.png"}]), "'data:' URI"),
	(envelope(commands=[{"x": [0, 1], "y": [0, 1]}]), "'op' is missing"),
	(envelope(commands="not-an-array"), "expected an array of commands"),
])
def test_validate_envelope_rejections(doc, fragment):
	''' Every §7.2 rejection, each with a pointed message. '''
	with pytest.raises(EnvelopeError) as raised:
		validate_envelope(doc)
	assert fragment in str(raised.value)


def test_validate_envelope_message_carries_a_path():
	''' Messages locate the offending field, JSON-pointer style. '''
	with pytest.raises(EnvelopeError) as raised:
		validate_envelope(envelope(commands=[{"op": "line", "x": [0.0, float("nan")], "y": [0.0, 1.0]}]))
	assert "/commands/0/x/1" in str(raised.value)


def test_validate_envelope_tolerates_unknown_minor_and_unknown_op():
	'''
	§3.7: a 1.0 consumer must accept an unknown *minor* and must not hard-fail
	on an unknown op (the renderer skips it with a diagnostic instead) — that is
	what keeps old viewers working against slightly-newer streams.
	'''
	validate_envelope(envelope(cornish_view="1.7",
	                           commands=[{"op": "hologram", "intensity": 3}]))


def test_validate_envelope_ignores_unknown_fields():
	''' §3.7: unknown fields are ignored, not rejected — that is how v1.x grows. '''
	doc = envelope()
	doc["some_future_field"] = {"anything": 1}
	doc["canvas"]["future_canvas_key"] = "ok"
	validate_envelope(doc)


def test_validate_envelope_rejects_boolean_where_a_number_belongs():
	'''
	``True`` is an ``int`` in Python; accepting it would silently mean 1 and
	corrupt a canonical hash.
	'''
	with pytest.raises(EnvelopeError):
		validate_envelope(envelope(commands=[{"op": "line", "x": [0.0, True], "y": [0.0, 1.0]}]))


# ---------------------------------------------------------------------------
# validate_document
# ---------------------------------------------------------------------------

def document(**overrides) -> dict:
	base = {
		"cornish_view": SCHEMA_VERSION,
		"type": "draw_stream",
		"canvas": dict(CANVAS),
		"layers": [{"id": "grid", "kind": "grid", "z": -10, "commands": []}],
	}
	for key, value in overrides.items():
		if value is REMOVE:
			base.pop(key, None)
		else:
			base[key] = value
	return base


def test_validate_document_accepts_well_formed():
	validate_document(document())


@pytest.mark.parametrize("doc, fragment", [
	(document(canvas=REMOVE), "requires a 'canvas'"),
	(document(layers=REMOVE), "requires a 'layers' array"),
	(document(layers=[{"kind": "grid", "commands": []}]), "needs an 'id'"),
	(document(layers=[{"id": "grid"}]), "needs a 'commands' array"),
	(document(layers=[{"id": "grid", "commands": []}, {"id": "grid", "commands": []}]), "duplicate layer id"),
	(document(layer={"id": "grid", "action": "replace"}), "must not carry a 'layer' key"),
	(document(layers=[{"id": "grid", "action": "replace", "commands": []}]), "'action' belongs to an envelope"),
	(document(layers="nope"), "expected an array"),
])
def test_validate_document_rejections(doc, fragment):
	with pytest.raises(EnvelopeError) as raised:
		validate_document(doc)
	assert fragment in str(raised.value)


def test_hybrid_body_fails_through_either_validator():
	'''
	§1.1 ownership rule: a body carrying BOTH ``layer`` and ``layers`` must be
	rejected no matter which validator is reached first — that is what makes
	the exactly-one-of rule dispatch-order independent.
	'''
	hybrid = document()
	hybrid["layer"] = {"id": "grid", "action": "replace"}
	with pytest.raises(EnvelopeError):
		validate_document(hybrid)
	with pytest.raises(EnvelopeError):
		validate_envelope(hybrid)


def test_body_with_neither_layer_nor_layers_and_no_canvas_is_invalid():
	naked = {"cornish_view": SCHEMA_VERSION, "type": "draw_stream"}
	with pytest.raises(EnvelopeError):
		validate_envelope(naked)
	with pytest.raises(EnvelopeError):
		validate_document(naked)


# ---------------------------------------------------------------------------
# canonical JSON (RFC 8785 / JCS)
# ---------------------------------------------------------------------------

def double_from_hex(bits: str) -> float:
	''' Decode an IEEE 754 bit pattern, the way RFC 8785's appendix states its vectors. '''
	return struct.unpack(">d", struct.pack(">Q", int(bits, 16)))[0]


#: (bit pattern, expected ES Number::toString output). Values verified against a
#: real JavaScript engine — the point of JCS here is cross-language agreement,
#: so the reference is the other language, not our own reading of the RFC.
JCS_NUMBER_VECTORS = [
	("0000000000000000", "0"),
	("8000000000000000", "0"),                       # -0 prints as 0 per JCS
	("0000000000000001", "5e-324"),                  # smallest subnormal
	("8000000000000001", "-5e-324"),
	("7fefffffffffffff", "1.7976931348623157e+308"), # largest finite double
	("4340000000000000", "9007199254740992"),        # 2**53
	("4430000000000000", "295147905179352830000"),
	("44b52d02c7e14af6", "1e+23"),
	("444b1ae4d6e2ef50", "1e+21"),                   # the decimal/exponential threshold
	("3eb0c6f7a0b5ed8d", "0.000001"),                # 1e-6 stays decimal
	("3ff0000000000000", "1"),                       # 1.0 -> "1", NOT "1.0"
	("3fb999999999999a", "0.1"),
	("412e848000000000", "1000000"),
	("3f847ae147ae147b", "0.01"),
]


@pytest.mark.parametrize("bits, expected", JCS_NUMBER_VECTORS)
def test_jcs_number_vectors(bits, expected):
	assert _jcs_number(double_from_hex(bits)) == expected


@pytest.mark.parametrize("value, expected", [
	(1.0, "1"),               # Python repr says '1.0'; JS says '1'
	(1e-07, "1e-7"),          # Python repr pads the exponent to '1e-07'
	(0.00001, "0.00001"),     # stays decimal (|v| >= 1e-6)
	(1e-06, "0.000001"),
	(1e-7, "1e-7"),
	(-0.0, "0"),
	(1e21, "1e+21"),
	(1e20, "100000000000000000000"),
	(1, "1"),                 # a Python int is a double to the schema
	(-1.5, "-1.5"),
	(123.456, "123.456"),
])
def test_jcs_number_divergence_cases(value, expected):
	''' The exact cases where Python's native formatting disagrees with JS (§3.6). '''
	assert _jcs_number(value) == expected


def test_jcs_number_rejects_non_finite():
	for value in (float("nan"), float("inf"), float("-inf")):
		with pytest.raises(EnvelopeError):
			_jcs_number(value)


def test_canonical_json_sorts_keys_and_omits_whitespace():
	assert canonical_json({"b": 1, "a": 2, "A": 3}) == '{"A":3,"a":2,"b":1}'


def test_canonical_json_passes_non_ascii_through():
	''' JCS keeps non-ASCII literal (UTF-8); only the C0 controls and " \\ are escaped. '''
	assert canonical_json({"t": "α δ ⊙"}) == '{"t":"α δ ⊙"}'
	assert canonical_json({"t": "a\nb\"c\\d"}) == '{"t":"a\\nb\\"c\\\\d"}'


def test_canonical_json_rejects_non_finite_and_alien_types():
	with pytest.raises(EnvelopeError):
		canonical_json({"x": float("nan")})
	with pytest.raises(EnvelopeError):
		canonical_json({"x": {1, 2}})


def test_canonical_json_is_parseable_json():
	doc = document()
	assert json.loads(canonical_json(doc)) == json.loads(json.dumps(doc))


#: Pinned canonical hashes, cross-checked into the JS conformance fixtures.
PINNED_CANVAS_ONLY = {
	"cornish_view": "1.0", "type": "draw_stream",
	"canvas": {"gbox": [0, 0, 1, 1], "aspect": 1.0, "margin": 0.15, "char_height": 0.02, "width": 640},
}
PINNED_CANVAS_ONLY_HASH = "58b9ed579dc52f5929b8f2f3292e08e0d2403af5aa5d58650f8f5cc254f57c58"

PINNED_DOCUMENT = {
	"cornish_view": "1.0", "type": "draw_stream", "plot_id": "demo",
	"canvas": {"gbox": [0, 0, 1, 0.5], "aspect": 0.5},
	"layers": [{"id": "grid", "kind": "grid", "z": -10, "commands": [
		{"op": "attr", "prim": "line", "colour": 11, "width": 0.5},
		{"op": "line", "x": [0.0, 1.0], "y": [0.25, 0.25]},
		{"op": "text", "text": "2:00", "x": 0.5, "y": -0.05, "just": "TC", "up": [0, 1]},
	]}],
}
PINNED_DOCUMENT_HASH = "dfd69b5a2e4f6cb41437fed730df2403367ebcf812561c45ba8aa43a1e7cd017"


def test_envelope_hash_is_pinned():
	'''
	These two hashes are the cross-language contract: the JS conformance suite
	reproduces them with its own ``canonicalJSON``. A change here is a schema
	break, not a refactor.
	'''
	validate_envelope(PINNED_CANVAS_ONLY)
	validate_document(PINNED_DOCUMENT)
	assert envelope_hash(PINNED_CANVAS_ONLY) == PINNED_CANVAS_ONLY_HASH
	assert envelope_hash(PINNED_DOCUMENT) == PINNED_DOCUMENT_HASH


def test_envelope_hash_is_insensitive_to_key_insertion_order():
	first = {"cornish_view": "1.0", "type": "draw_stream", "canvas": dict(CANVAS)}
	second = {"canvas": dict(CANVAS), "type": "draw_stream", "cornish_view": "1.0"}
	assert envelope_hash(first) == envelope_hash(second)


def test_canonical_json_collapses_integral_floats():
	''' 1.0 and 1 hash identically — the int/float distinction is not semantic here. '''
	assert envelope_hash({"a": 1.0}) == envelope_hash({"a": 1})


# ---------------------------------------------------------------------------
# make_envelope
# ---------------------------------------------------------------------------

def test_make_envelope_layer_form():
	built = make_envelope(layer_id="region-1", action=LayerAction.REPLACE, kind="region", z=0,
	                      commands=[{"op": "line", "x": [0, 1], "y": [0, 1]}], canvas=dict(CANVAS))
	validate_envelope(built)
	assert built["layer"] == {"id": "region-1", "action": "replace", "kind": "region", "z": 0}
	assert built["cornish_view"] == SCHEMA_VERSION


def test_make_envelope_canvas_only_form():
	built = make_envelope(canvas=dict(CANVAS))
	validate_envelope(built)
	assert "layer" not in built and "commands" not in built


def test_make_envelope_delete_form():
	built = make_envelope(layer_id="grid", action=LayerAction.DELETE)
	validate_envelope(built)
	assert "commands" not in built


def test_make_envelope_accepts_empty_command_list():
	built = make_envelope(layer_id="grid", action=LayerAction.REPLACE, commands=[])
	assert built["commands"] == []


@pytest.mark.parametrize("kwargs", [
	{"layer_id": "grid"},                                              # action missing
	{"action": LayerAction.REPLACE},                                   # layer_id missing
	{},                                                                 # no legal shape at all
	{"canvas": dict(CANVAS), "commands": [{"op": "line", "x": [0, 1], "y": [0, 1]}]},  # canvas-only + commands
	{"canvas": dict(CANVAS), "kind": "grid"},                          # kind without a layer
	{"layer_id": "grid", "action": LayerAction.DELETE,
	 "commands": [{"op": "line", "x": [0, 1], "y": [0, 1]}]},          # commands on delete
	{"layer_id": "grid", "action": LayerAction.REPLACE},                # replace without commands
])
def test_make_envelope_illegal_shapes_raise(kwargs):
	with pytest.raises(ValueError):
		make_envelope(**kwargs)


def test_make_envelope_plot_id_is_optional():
	''' §3.1: pattern-3 embedders omit plot_id. '''
	assert "plot_id" not in make_envelope(canvas=dict(CANVAS))
	assert make_envelope(canvas=dict(CANVAS), plot_id="demo")["plot_id"] == "demo"
