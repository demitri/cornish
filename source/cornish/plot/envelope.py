
'''
Envelope schema v1.0 — the codec (SPEC-09A §3).

An *envelope* is one self-contained message about one layer of one plot; a
*document* is the whole plot (canvas + every layer in render order). The two
shapes are distinguished by exactly one key — ``layer`` vs ``layers`` — and
each validator here treats the *other* shape's marker key as a hard error, so
a hybrid body fails no matter which validator sees it first (§3.1).

This module is **pure**: no I/O, no transport, no AST. It is the contract
shared by cornish's own producers and consumers and by any other
implementation of the schema, so its rules are written to be mechanical.

Canonicalization is RFC 8785 (JCS), not ``json.dumps(sort_keys=True)``:
Python's native float formatting disagrees with JavaScript's on essentially
every envelope (``1.0`` vs ``1``, ``1e-07`` vs ``1e-7``), which would make
pinned hashes irreproducible in the JS consumers (§3.6).
'''

import hashlib
import json
import math
import re
from enum import StrEnum
from typing import Optional

from ..exc import EnvelopeError
from .palette import normalize_palette_entry

__all__ = ['SCHEMA_VERSION', 'EnvelopeError', 'LayerAction', 'GrfOp', 'PrimType',
           'validate_envelope', 'validate_document', 'canonical_json', 'envelope_hash',
           'make_envelope', 'ID_PATTERN']

#: The schema version this codec produces. Consumers MUST reject an unknown
#: *major* and MUST accept an unknown *minor* while ignoring unknown fields.
SCHEMA_VERSION = "1.0"

SCHEMA_MAJOR = 1

DOCUMENT_TYPE = "draw_stream"

#: Legal plot ids and layer ids (§3.1, §3.4).
ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]{1,64}$")

_VERSION_PATTERN = re.compile(r"^(\d+)\.(\d+)$")

_JUST_VERTICAL = "TCBM"
_JUST_HORIZONTAL = "LCR"

#: attr-state defaults at the start of every layer (§3.3).
ATTR_DEFAULTS = {"colour": 1, "width": 1.0, "style": 1, "size": 1.0, "font": 1}


class LayerAction(StrEnum):
	''' What an envelope does to the layer it addresses (§3.4). '''
	REPLACE = "replace"
	APPEND = "append"
	DELETE = "delete"


class GrfOp(StrEnum):
	''' Command discriminants in a draw stream (§3.3). '''
	ATTR = "attr"
	LINE = "line"
	MARK = "mark"
	TEXT = "text"
	IMAGE = "image"


class PrimType(StrEnum):
	''' The primitive classes an ``attr`` command can address (§3.3). '''
	LINE = "line"
	MARK = "mark"
	TEXT = "text"


# ---------------------------------------------------------------------------
# RFC 8785 (JCS) canonicalization
# ---------------------------------------------------------------------------

def _jcs_number(value) -> str:
	'''
	Serialize *value* exactly as ECMAScript ``Number::toString`` would — the
	number grammar RFC 8785 mandates (§3.6).

	Python's ``repr`` already yields the shortest round-tripping digit string;
	this function only re-places the decimal point according to the ES rules,
	so ``1.0`` becomes ``1``, ``1e-07`` becomes ``1e-7``, ``0.00001`` stays
	decimal, and ``-0.0`` becomes ``0``.
	'''
	if isinstance(value, bool):
		# guarded by the callers too, but bool is an int subclass and a
		# silently-numeric True would corrupt a hash
		raise EnvelopeError("A boolean cannot be serialized as a number.")
	if isinstance(value, int):
		if abs(value) > 2 ** 53:
			raise EnvelopeError(f"Integer {value} exceeds the exact range of an IEEE 754 double; "
			                    "envelope numbers are doubles (RFC 8785).")
		value = float(value)
	number = float(value)

	if not math.isfinite(number):
		raise EnvelopeError(f"Non-finite number {number!r} cannot appear in an envelope "
		                    "(NaN/Infinity are illegal in JSON and in the schema).")
	if number == 0:
		return "0"                              # also collapses -0.0, per JCS
	if number < 0:
		return "-" + _jcs_number(-number)

	representation = repr(number)
	if "e" in representation or "E" in representation:
		mantissa, _, exponent_text = representation.partition("e")
		if not exponent_text:
			mantissa, _, exponent_text = representation.partition("E")
		exponent = int(exponent_text)
	else:
		mantissa, exponent = representation, 0

	integer_part, _, fraction_part = mantissa.partition(".")
	all_digits = integer_part + fraction_part
	stripped = all_digits.lstrip("0")
	leading_zeros = len(all_digits) - len(stripped)
	# value == 0.<digits> x 10**point, with `digits` carrying no leading zero
	point = len(integer_part) + exponent - leading_zeros
	digits = stripped.rstrip("0") or "0"
	count = len(digits)

	if count <= point <= 21:
		return digits + "0" * (point - count)
	if 0 < point <= 21:
		return digits[:point] + "." + digits[point:]
	if -6 < point <= 0:
		return "0." + "0" * (-point) + digits
	# exponential notation; the ES exponent is point - 1
	exponent_value = point - 1
	sign = "+" if exponent_value >= 0 else "-"
	significand = digits if count == 1 else digits[0] + "." + digits[1:]
	return f"{significand}e{sign}{abs(exponent_value)}"


def _jcs_string(text: str) -> str:
	'''
	JCS §3.2.2.2 string escaping. Python's ``json.dumps`` with
	``ensure_ascii=False`` implements exactly this set: escape ``"`` and ``\\``,
	use the short forms for the five named control characters, ``\\u00XX`` for
	the remaining C0 controls, and leave every other code point literal.
	'''
	return json.dumps(text, ensure_ascii=False)


def _canonicalize(value, path: str) -> str:
	'''
	Recursively canonicalize *value*; *path* is a JSON-pointer-ish location used
	in error messages.
	'''
	if value is None:
		return "null"
	if isinstance(value, bool):
		return "true" if value else "false"
	if isinstance(value, str):
		return _jcs_string(value)
	if isinstance(value, (int, float)):
		try:
			return _jcs_number(value)
		except EnvelopeError as error:
			raise EnvelopeError(f"{path}: {error}") from error
	if isinstance(value, (list, tuple)):
		return "[" + ",".join(_canonicalize(item, f"{path}/{index}")
		                      for index, item in enumerate(value)) + "]"
	if isinstance(value, dict):
		parts = []
		for key in _sorted_keys(value, path):
			parts.append(_jcs_string(key) + ":" + _canonicalize(value[key], f"{path}/{key}"))
		return "{" + ",".join(parts) + "}"
	raise EnvelopeError(f"{path}: values of type {type(value).__name__} cannot appear in an "
	                    "envelope (allowed: object, array, string, number, boolean, null).")


def _sorted_keys(mapping: dict, path: str):
	'''
	Sort object keys by their UTF-16 code units, as RFC 8785 requires. (Python's
	native str ordering is by code point, which differs from UTF-16 ordering for
	astral characters — envelope keys are ASCII, but the codec is a contract
	other implementations are checked against, so it follows the RFC exactly.)
	'''
	for key in mapping:
		if not isinstance(key, str):
			raise EnvelopeError(f"{path}: object keys must be strings; got {type(key).__name__}.")
	return sorted(mapping, key=lambda key: key.encode("utf-16-be"))


def canonical_json(doc: dict) -> str:
	'''
	The RFC 8785 (JCS) canonical serialization of *doc* (§3.6).

	:raises EnvelopeError: on a non-finite number or a value of an unserializable type
	'''
	return _canonicalize(doc, "")


def envelope_hash(doc: dict) -> str:
	''' SHA-256 hexdigest of the canonical UTF-8 bytes of *doc* (§3.6). '''
	return hashlib.sha256(canonical_json(doc).encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# validation helpers
# ---------------------------------------------------------------------------

def _fail(path: str, message: str):
	raise EnvelopeError(f"{path or '/'}: {message}")


def _require_mapping(value, path: str) -> dict:
	if not isinstance(value, dict):
		_fail(path, f"expected an object; got {type(value).__name__}.")
	return value


def _require_number(value, path: str) -> float:
	'''
	A finite JSON number. Booleans are rejected outright: ``True`` is an ``int``
	in Python, and letting it through would silently mean 1.
	'''
	if isinstance(value, bool) or not isinstance(value, (int, float)):
		_fail(path, f"expected a number; got {type(value).__name__}.")
	if not math.isfinite(value):
		_fail(path, f"expected a finite number; got {value!r} "
		            "(NaN/Infinity are illegal — an Ast.BAD-contaminated point must never be emitted).")
	return float(value)


def _require_int(value, path: str) -> int:
	if isinstance(value, bool) or not isinstance(value, int):
		_fail(path, f"expected an integer; got {type(value).__name__}.")
	return value


def _require_string(value, path: str) -> str:
	if not isinstance(value, str):
		_fail(path, f"expected a string; got {type(value).__name__}.")
	return value


def _require_id(value, path: str) -> str:
	text = _require_string(value, path)
	if not ID_PATTERN.match(text):
		_fail(path, f"'{text}' is not a legal id (must match {ID_PATTERN.pattern}).")
	return text


def _check_version(doc: dict, path: str) -> None:
	if "cornish_view" not in doc:
		_fail(path, "the required field 'cornish_view' (schema version) is missing.")
	version = _require_string(doc["cornish_view"], f"{path}/cornish_view")
	match = _VERSION_PATTERN.match(version)
	if match is None:
		_fail(f"{path}/cornish_view", f"'{version}' is not a '<major>.<minor>' version string.")
	major = int(match.group(1))
	if major != SCHEMA_MAJOR:
		# a MAJOR bump changes existing semantics, so refusing is the whole point;
		# an unknown MINOR is accepted (§3.7) and its unknown fields ignored
		_fail(f"{path}/cornish_view", f"stream requires cornish-view {major}.x; "
		                              f"this codec implements {SCHEMA_MAJOR}.x.")


def _check_type(doc: dict, path: str) -> None:
	if "type" not in doc:
		_fail(path, "the required field 'type' is missing.")
	document_type = _require_string(doc["type"], f"{path}/type")
	if document_type != DOCUMENT_TYPE:
		_fail(f"{path}/type", f"expected '{DOCUMENT_TYPE}'; got '{document_type}' "
		                      "(a plot_spec is a different document and is never mixed in).")


def _check_canvas(canvas, path: str) -> None:
	''' §3.2. Unknown keys are ignored — that is how v1.x minor growth works. '''
	canvas = _require_mapping(canvas, path)

	gbox = None
	if "gbox" in canvas:
		raw = canvas["gbox"]
		if not isinstance(raw, (list, tuple)) or len(raw) != 4:
			_fail(f"{path}/gbox", "expected [x1, y1, x2, y2].")
		gbox = [_require_number(value, f"{path}/gbox/{index}") for index, value in enumerate(raw)]
		if gbox[2] == gbox[0] or gbox[3] == gbox[1]:
			_fail(f"{path}/gbox", f"the graphics box must have non-zero extent in both axes; got {gbox}.")

	aspect = None
	if "aspect" in canvas:
		aspect = _require_number(canvas["aspect"], f"{path}/aspect")
		if aspect <= 0:
			_fail(f"{path}/aspect", f"must be > 0; got {aspect}.")

	if gbox is not None and aspect is not None:
		implied = (gbox[3] - gbox[1]) / (gbox[2] - gbox[0])
		if not math.isclose(implied, aspect, rel_tol=1e-9, abs_tol=1e-12):
			_fail(f"{path}/aspect", f"{aspect} contradicts the gbox, which implies {implied} "
			                        "(they must agree when both are present).")

	if "margin" in canvas:
		margin = _require_number(canvas["margin"], f"{path}/margin")
		if not 0 <= margin <= 0.5:
			_fail(f"{path}/margin", f"must lie in [0, 0.5]; got {margin}.")

	if "char_height" in canvas:
		char_height = _require_number(canvas["char_height"], f"{path}/char_height")
		if char_height <= 0:
			_fail(f"{path}/char_height", f"must be > 0; got {char_height}.")

	if "font_family" in canvas:
		_require_string(canvas["font_family"], f"{path}/font_family")

	if "palette" in canvas:
		palette = canvas["palette"]
		if not isinstance(palette, (list, tuple)):
			_fail(f"{path}/palette", f"expected an array; got {type(palette).__name__}.")
		if len(palette) < 2:
			_fail(f"{path}/palette", f"needs at least 2 entries (0 = background, 1 = foreground); got {len(palette)}.")
		for index, entry in enumerate(palette):
			entry_path = f"{path}/palette/{index}"
			_require_string(entry, entry_path)
			try:
				normalize_palette_entry(entry)
			except (TypeError, ValueError) as error:
				_fail(entry_path, str(error))

	if "background" in canvas and canvas["background"] is not None:
		background = _require_string(canvas["background"], f"{path}/background")
		try:
			normalize_palette_entry(background)
		except (TypeError, ValueError) as error:
			_fail(f"{path}/background", str(error))

	if "width" in canvas:
		width = _require_int(canvas["width"], f"{path}/width")
		if width <= 0:
			_fail(f"{path}/width", f"must be a positive integer; got {width}.")


def _check_point_arrays(command: dict, path: str, minimum: int) -> None:
	for axis in ("x", "y"):
		if axis not in command:
			_fail(path, f"the required field '{axis}' is missing.")
		if not isinstance(command[axis], (list, tuple)):
			_fail(f"{path}/{axis}", f"expected an array; got {type(command[axis]).__name__}.")
	length_x, length_y = len(command["x"]), len(command["y"])
	if length_x != length_y:
		_fail(path, f"'x' and 'y' must be the same length; got {length_x} and {length_y}.")
	if length_x < minimum:
		_fail(path, f"needs at least {minimum} point(s); got {length_x}.")
	for axis in ("x", "y"):
		for index, value in enumerate(command[axis]):
			_require_number(value, f"{path}/{axis}/{index}")


def _check_command(command, path: str) -> None:
	'''
	§3.3. Known ops are validated strictly; an **unknown op is accepted** here
	because §3.7 requires consumers to skip it with a diagnostic rather than
	reject the stream — that is how a 1.0 consumer survives a 1.1 producer.
	Its numbers are still checked for finiteness by the caller's sweep.
	'''
	command = _require_mapping(command, path)
	if "op" not in command:
		_fail(path, "the required field 'op' is missing.")
	op = _require_string(command["op"], f"{path}/op")

	if op == GrfOp.ATTR:
		if "prim" not in command:
			_fail(path, "an 'attr' command needs a 'prim' field.")
		prim = _require_string(command["prim"], f"{path}/prim")
		if prim not in tuple(PrimType):
			_fail(f"{path}/prim", f"expected one of {[str(value) for value in PrimType]}; got '{prim}'.")
		# Out-of-RANGE values are NOT rejected: §3.3 freezes visible render-time
		# fallbacks for them (out-of-range colour -> foreground, style -> solid,
		# font -> normal), which also covers a canvas re-theme stranding indices.
		if "colour" in command:
			colour = command["colour"]
			if isinstance(colour, bool) or not isinstance(colour, (int, str)):
				_fail(f"{path}/colour", "expected a base-palette index (int) or a CSS colour string; "
				                        f"got {type(colour).__name__}.")
		for key in ("width", "size"):
			if key in command:
				_require_number(command[key], f"{path}/{key}")
		for key in ("style", "font"):
			if key in command:
				_require_int(command[key], f"{path}/{key}")

	elif op == GrfOp.LINE:
		_check_point_arrays(command, path, minimum=2)

	elif op == GrfOp.MARK:
		_check_point_arrays(command, path, minimum=1)
		if "marker" not in command:
			_fail(path, "a 'mark' command needs a 'marker' field.")
		_require_int(command["marker"], f"{path}/marker")

	elif op == GrfOp.TEXT:
		for key in ("text", "just"):
			if key not in command:
				_fail(path, f"a 'text' command needs a '{key}' field.")
		_require_string(command["text"], f"{path}/text")
		just = _require_string(command["just"], f"{path}/just")
		if len(just) != 2 or just[0] not in _JUST_VERTICAL or just[1] not in _JUST_HORIZONTAL:
			_fail(f"{path}/just", f"'{just}' is not a legal justification "
			                      f"(vertical {tuple(_JUST_VERTICAL)} then horizontal {tuple(_JUST_HORIZONTAL)}).")
		for key in ("x", "y"):
			if key not in command:
				_fail(path, f"a 'text' command needs a '{key}' field.")
			_require_number(command[key], f"{path}/{key}")
		if "up" not in command:
			_fail(path, "a 'text' command needs an 'up' field.")
		up = command["up"]
		if not isinstance(up, (list, tuple)) or len(up) != 2:
			_fail(f"{path}/up", "expected a 2-element up-vector [ux, uy].")
		up_x = _require_number(up[0], f"{path}/up/0")
		up_y = _require_number(up[1], f"{path}/up/1")
		if up_x == 0 and up_y == 0:
			_fail(f"{path}/up", "the up-vector must be non-zero.")

	elif op == GrfOp.IMAGE:
		for key in ("x1", "y1", "x2", "y2"):
			if key not in command:
				_fail(path, f"an 'image' command needs a '{key}' field.")
			_require_number(command[key], f"{path}/{key}")
		if "href" not in command:
			_fail(path, "an 'image' command needs an 'href' field.")
		href = _require_string(command["href"], f"{path}/href")
		if not href.startswith("data:"):
			# envelopes are self-contained documents; a remote fetch would make
			# replay impure and turn a plot into a tracking vector
			_fail(f"{path}/href", "must be a 'data:' URI — remote URLs are rejected so that an "
			                      "envelope stays a self-contained, replayable document.")


def _check_commands(commands, path: str) -> None:
	if not isinstance(commands, (list, tuple)):
		_fail(path, f"expected an array of commands; got {type(commands).__name__}.")
	for index, command in enumerate(commands):
		_check_command(command, f"{path}/{index}")


def _check_layer_object(layer, path: str) -> str:
	'''
	Validate the ``layer`` object of an envelope (§3.4); returns its action.
	'''
	layer = _require_mapping(layer, path)
	if "id" not in layer:
		_fail(path, "the required field 'id' is missing.")
	_require_id(layer["id"], f"{path}/id")
	if "action" not in layer:
		_fail(path, "the required field 'action' is missing.")
	action = _require_string(layer["action"], f"{path}/action")
	if action not in tuple(LayerAction):
		_fail(f"{path}/action", f"expected one of {[str(value) for value in LayerAction]}; got '{action}'.")
	if "kind" in layer:
		_require_string(layer["kind"], f"{path}/kind")
	if "z" in layer:
		_require_int(layer["z"], f"{path}/z")
	return action


# ---------------------------------------------------------------------------
# the two public validators
# ---------------------------------------------------------------------------

def validate_envelope(doc: dict) -> None:
	'''
	Validate a single **envelope** (§3.1).

	A body carrying a ``layers`` key is rejected outright: that is the document
	shape, and each validator owning the other's marker key as a hard error is
	what makes §3.1's exactly-one-of rule hold regardless of dispatch order.

	:raises EnvelopeError: with a JSON-pointer-ish path to the offending field
	'''
	doc = _require_mapping(doc, "")
	if "layers" in doc:
		_fail("/layers", "an envelope must not carry a 'layers' key — that is the document shape "
		                 "(exactly one of 'layer'/'layers', §3.1).")
	_check_version(doc, "")
	_check_type(doc, "")

	if "plot_id" in doc:
		_require_id(doc["plot_id"], "/plot_id")
	if "canvas" in doc:
		_check_canvas(doc["canvas"], "/canvas")

	has_layer = "layer" in doc
	has_commands = "commands" in doc

	if not has_layer:
		# the only legal layer-less body is the canvas-only envelope
		if "canvas" not in doc:
			_fail("", "the required field 'layer' is missing (only a canvas-only envelope — "
			          "'canvas' present, 'commands' absent — may omit it).")
		if has_commands:
			_fail("/commands", "a canvas-only envelope must not carry 'commands' "
			                   "(add a 'layer' object to draw something).")
		return

	action = _check_layer_object(doc["layer"], "/layer")

	if action == LayerAction.DELETE:
		if has_commands:
			_fail("/commands", "a 'delete' envelope must not carry 'commands'.")
		return

	if not has_commands:
		_fail("", f"the required field 'commands' is missing for a '{action}' envelope "
		          "(an empty array is legal and clears the layer).")
	_check_commands(doc["commands"], "/commands")


def validate_document(doc: dict) -> None:
	'''
	Validate a full **document** — the ``layers`` form (§3.1), which is what
	``SkyPlot.document()``, ``save('.json')`` and ``GET /plot/{id}.json``
	produce, and what a client may POST for an atomic whole-plot replace.

	A body carrying a ``layer`` key is rejected (mirror of
	:func:`validate_envelope`'s rule).

	:raises EnvelopeError: with a JSON-pointer-ish path to the offending field
	'''
	doc = _require_mapping(doc, "")
	if "layer" in doc:
		_fail("/layer", "a document must not carry a 'layer' key — that is the envelope shape "
		                "(exactly one of 'layer'/'layers', §3.1).")
	_check_version(doc, "")
	_check_type(doc, "")

	if "plot_id" in doc:
		_require_id(doc["plot_id"], "/plot_id")

	if "canvas" not in doc:
		_fail("", "a document requires a 'canvas' object.")
	_check_canvas(doc["canvas"], "/canvas")

	if "layers" not in doc:
		_fail("", "a document requires a 'layers' array (an envelope uses 'layer' instead).")
	layers = doc["layers"]
	if not isinstance(layers, (list, tuple)):
		_fail("/layers", f"expected an array; got {type(layers).__name__}.")

	seen = {}
	for index, layer in enumerate(layers):
		path = f"/layers/{index}"
		layer = _require_mapping(layer, path)
		if "id" not in layer:
			_fail(path, "every layer needs an 'id'.")
		layer_id = _require_id(layer["id"], f"{path}/id")
		if layer_id in seen:
			_fail(f"{path}/id", f"duplicate layer id '{layer_id}' (already used at /layers/{seen[layer_id]}).")
		seen[layer_id] = index
		if "kind" in layer:
			_require_string(layer["kind"], f"{path}/kind")
		if "z" in layer:
			_require_int(layer["z"], f"{path}/z")
		if "action" in layer:
			_fail(f"{path}/action", "a document's layers carry content, not actions — "
			                        "'action' belongs to an envelope's 'layer' object.")
		if "commands" not in layer:
			_fail(path, "every layer needs a 'commands' array (an empty array is legal).")
		_check_commands(layer["commands"], f"{path}/commands")


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------

def make_envelope(*, layer_id: Optional[str] = None,
                  action: Optional[LayerAction] = None,
                  commands: Optional[list] = None,
                  canvas: Optional[dict] = None,
                  plot_id: Optional[str] = None,
                  kind: Optional[str] = None,
                  z: Optional[int] = None) -> dict:
	'''
	Build one envelope (§3.1). The output always passes
	:func:`validate_envelope`.

	``layer_id`` and ``action`` go together; both ``None`` **with** a canvas
	produces the canvas-only envelope; both ``None`` **without** a canvas has no
	legal shape at all.

	:raises ValueError: if the arguments describe no legal envelope shape
	:raises EnvelopeError: if the assembled envelope fails validation (a bug in the caller's field values)
	'''
	if (layer_id is None) != (action is None):
		raise ValueError("'layer_id' and 'action' must be given together "
		                 f"(got layer_id={layer_id!r}, action={action!r}).")

	envelope = {"cornish_view": SCHEMA_VERSION, "type": DOCUMENT_TYPE}
	if plot_id is not None:
		envelope["plot_id"] = plot_id
	if canvas is not None:
		envelope["canvas"] = canvas

	if layer_id is None:
		if canvas is None:
			raise ValueError("An envelope with no layer must carry a canvas "
			                 "(the canvas-only, plot-creating form); there is no other legal shape.")
		if commands is not None:
			raise ValueError("A canvas-only envelope cannot carry commands; "
			                 "supply 'layer_id' and 'action' to draw something.")
		if kind is not None or z is not None:
			raise ValueError("'kind'/'z' describe a layer, so they need 'layer_id' and 'action'.")
		validate_envelope(envelope)
		return envelope

	action = LayerAction(action)
	layer = {"id": layer_id, "action": str(action)}
	if kind is not None:
		layer["kind"] = kind
	if z is not None:
		layer["z"] = z
	envelope["layer"] = layer

	if action == LayerAction.DELETE:
		if commands is not None:
			raise ValueError("A 'delete' envelope cannot carry commands.")
	else:
		if commands is None:
			raise ValueError(f"A '{action}' envelope requires a commands list "
			                 "(an empty list is legal and clears the layer).")
		envelope["commands"] = list(commands)

	validate_envelope(envelope)
	return envelope
