
'''
SVG writer purity, determinism and well-formedness (SPEC-09A §7.4).

Also home to the **render/serve purity gate** — the machine-enforced rule that
nothing under ``cornish/plot/`` may reach for a server or transport library.
That rule is what keeps ``render_svg`` usable as a pure call in Flask pattern 1
(the recommended production shape), and it joins the existing gate family
(``test_no_kwargs_into_pyast_constructors``, ``test_no_conversions_outside_bridge``).

The golden draw streams are exercised in ``test_golden_streams.py``.
'''

import ast as python_ast
import base64
import glob
import logging
import os
import struct
import subprocess
import sys
import xml.etree.ElementTree as ElementTree
import zlib

import pytest

from cornish.plot.envelope import envelope_hash, validate_document
from cornish.plot.svg import (SVGRenderer, array_to_png_data_uri, encode_greyscale_png,
                              render_svg)

SVG_NS = "{http://www.w3.org/2000/svg}"


def document(commands, *, canvas=None, layers=None) -> dict:
	doc = {
		"cornish_view": "1.0", "type": "draw_stream",
		"canvas": canvas or {"gbox": [0, 0, 1, 1], "aspect": 1.0, "margin": 0.15, "char_height": 0.02},
		"layers": layers if layers is not None else [{"id": "layer-1", "commands": commands}],
	}
	validate_document(doc)
	return doc


LINE_DOCUMENT = document([{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}])


# ---------------------------------------------------------------------------
# purity gate
# ---------------------------------------------------------------------------

#: Modules no part of the render path may import — at module level OR inside a
#: function body (SPEC-09A §1, §7.4).
BANNED_IMPORTS = ("http", "socketserver", "flask", "aiohttp", "asyncio", "socket",
                  "urllib.request", "cornish.view")

#: The ONE sanctioned exception (§1): browser display and its live-mode POST
#: helper use function-local transport imports, mirroring to_matplotlib()'s
#: "imports matplotlib here, and only here" pattern. The gate names them so a
#: THIRD occurrence fails.
TRANSPORT_ALLOWLIST = {("sky.py", "SkyPlot.show"), ("sky.py", "SkyPlot._live_post")}


def plot_package_directory() -> str:
	return os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "cornish", "plot"))


def _imported_modules(node, module_package: str):
	'''
	The absolute module names an Import/ImportFrom node brings in, resolving
	relative imports against *module_package* (so ``from ..view import x`` inside
	``cornish/plot/`` resolves to ``cornish.view`` and cannot hide from the gate).
	'''
	if isinstance(node, python_ast.Import):
		return [alias.name for alias in node.names]
	if isinstance(node, python_ast.ImportFrom):
		if node.level:
			parts = module_package.split(".")
			base = ".".join(parts[:len(parts) - node.level + 1])
			module = f"{base}.{node.module}" if node.module else base
		else:
			module = node.module or ""
		return [module] + [f"{module}.{alias.name}" for alias in node.names]
	return []


def _enclosing_function(path_stack) -> str:
	names = [node.name for node in path_stack
	         if isinstance(node, (python_ast.FunctionDef, python_ast.AsyncFunctionDef,
	                              python_ast.ClassDef))]
	return ".".join(names)


def test_plot_package_imports_no_transport():
	'''
	The render/serve separation, machine-enforced by AST walk over every module
	under ``cornish/plot/`` — module-level and function-local alike, since a
	function-local ``import flask`` would defeat a module-level-only check.

	Blind spots, stated so they are not assumed covered: a dynamic
	``importlib.import_module(name)`` with a computed name, and a transport
	object handed in by a caller. Neither shape exists today; both remain a
	human-review duty.
	'''
	offenders = []
	for path in glob.glob(os.path.join(plot_package_directory(), "**", "*.py"), recursive=True):
		relative = os.path.relpath(path, plot_package_directory())
		with open(path) as source_file:
			tree = python_ast.parse(source_file.read(), filename=path)

		stack = []

		def visit(node, stack=stack, relative=relative):
			for child in python_ast.iter_child_nodes(node):
				if isinstance(child, (python_ast.Import, python_ast.ImportFrom)):
					for module in _imported_modules(child, "cornish.plot"):
						for banned in BANNED_IMPORTS:
							if module == banned or module.startswith(banned + "."):
								if (relative, _enclosing_function(stack)) in TRANSPORT_ALLOWLIST:
									continue
								offenders.append(
									f"{relative}:{child.lineno}: imports '{module}' "
									f"(banned; enclosing function "
									f"'{_enclosing_function(stack) or '<module level>'}')")
				stack.append(child)
				visit(child)
				stack.pop()

		visit(tree)

	assert not offenders, ("transport imports found in the render path:\n" + "\n".join(offenders))


#: What the RUNTIME half of the gate can honestly assert.
#:
#: The spec asks that ``import cornish.plot`` pull in no banned module. Three of
#: them arrive whatever cornish does: ``astropy.coordinates`` — a hard
#: dependency of the SPEC-04 bridge — imports ``http``, ``http.client``,
#: ``socket`` and ``urllib.request`` at import time (verified directly: they are
#: absent after ``import numpy`` and after ``import astropy.units``, and present
#: after ``import astropy.coordinates``). Their presence in ``sys.modules``
#: therefore says nothing about cornish's own code.
#:
#: So the runtime check covers what it can actually attribute to cornish: the
#: view package and the web frameworks, none of which any dependency drags in.
#: Enforcement of cornish's own source is the AST-walk gate above, which is not
#: weakened by this — do not "fix" this list by deleting the gate.
RUNTIME_BANNED_IMPORTS = ("socketserver", "flask", "aiohttp", "cornish.view")


def test_importing_cornish_plot_pulls_in_no_transport():
	'''
	The runtime half of the gate: importing the render path must not drag the
	view package or a web framework into ``sys.modules``. Run in a subprocess so
	the check sees a clean interpreter rather than whatever this test session
	has already imported.
	'''
	source_directory = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
	program = (
		"import sys\n"
		"import cornish.plot\n"
		"banned = %r\n"
		"found = sorted(name for name in sys.modules\n"
		"               for b in banned if name == b or name.startswith(b + '.'))\n"
		"print(','.join(found))\n" % (RUNTIME_BANNED_IMPORTS,)
	)
	result = subprocess.run([sys.executable, "-c", program], cwd=source_directory,
	                        capture_output=True, text=True, check=True)
	assert result.stdout.strip() == "", f"cornish.plot imported transport modules: {result.stdout}"


def test_astropy_is_the_source_of_the_stdlib_transport_imports():
	'''
	Pins the reason ``RUNTIME_BANNED_IMPORTS`` is narrower than
	``BANNED_IMPORTS``, so the exemption stays justified rather than assumed. If
	astropy ever stops importing these, this test fails and the runtime check
	can be widened back.
	'''
	program = (
		"import sys\n"
		"import astropy.coordinates\n"
		"names = ('http', 'socket', 'urllib.request')\n"
		"print(','.join(sorted(n for n in names if n in sys.modules)))\n"
	)
	result = subprocess.run([sys.executable, "-c", program], capture_output=True,
	                        text=True, check=True)
	assert result.stdout.strip() == "http,socket,urllib.request"


# ---------------------------------------------------------------------------
# determinism / ETag
# ---------------------------------------------------------------------------

def test_render_is_byte_stable():
	''' The ETag pattern depends on this. '''
	assert render_svg(LINE_DOCUMENT) == render_svg(LINE_DOCUMENT)


def test_render_is_stable_across_equal_documents():
	first = document([{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}])
	second = document([{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}])
	assert render_svg(first) == render_svg(second)
	assert envelope_hash(first) == envelope_hash(second)


def test_document_key_order_does_not_change_the_hash_or_the_render():
	reordered = {key: LINE_DOCUMENT[key] for key in reversed(list(LINE_DOCUMENT))}
	assert envelope_hash(reordered) == envelope_hash(LINE_DOCUMENT)
	assert render_svg(reordered) == render_svg(LINE_DOCUMENT)


# ---------------------------------------------------------------------------
# well-formedness
# ---------------------------------------------------------------------------

def test_output_parses_as_xml_and_has_no_nan():
	svg = render_svg(LINE_DOCUMENT)
	root = ElementTree.fromstring(svg)
	assert root.tag == f"{SVG_NS}svg"
	assert "NaN" not in svg and "Infinity" not in svg


def test_svg_is_responsive():
	''' §4.3: viewBox + width="100%" so it scales inside its container. '''
	root = ElementTree.fromstring(render_svg(LINE_DOCUMENT))
	assert root.get("width") == "100%"
	assert root.get("viewBox") == "0 0 640 640"


def test_layer_groups_carry_data_attributes_and_no_id():
	'''
	§3.5: no ``id`` attribute at all. The legal layer-id charset contains ``.``
	and ``-``, which break ``#...`` CSS selectors — lookup is by attribute
	selector, where a quoted value needs no escaping.
	'''
	doc = document(None, layers=[{"id": "grid", "kind": "grid", "commands": []}])
	group = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}g")
	assert group.get("data-layer-id") == "grid"
	assert group.get("data-kind") == "grid"
	assert group.get("id") is None


def test_layer_without_a_kind_omits_the_attribute_entirely():
	''' §3.4: never ``data-kind=""``. '''
	doc = document(None, layers=[{"id": "bare", "commands": []}])
	group = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}g")
	assert "data-kind" not in group.attrib


def test_layers_render_in_z_order_then_creation_order():
	doc = document(None, layers=[
		{"id": "top", "z": 5, "commands": []},
		{"id": "back", "z": -20, "commands": []},
		{"id": "mid-a", "commands": []},
		{"id": "mid-b", "commands": []},
	])
	root = ElementTree.fromstring(render_svg(doc))
	order = [group.get("data-layer-id") for group in root.findall(f".//{SVG_NS}g")]
	assert order == ["back", "mid-a", "mid-b", "top"]


def test_polyline_and_point_counts_match_the_document():
	commands = [
		{"op": "line", "x": [0.1, 0.5, 0.9], "y": [0.1, 0.5, 0.9]},
		{"op": "line", "x": [0.2, 0.8], "y": [0.8, 0.2]},
	]
	root = ElementTree.fromstring(render_svg(document(commands)))
	polylines = root.findall(f".//{SVG_NS}polyline")
	assert len(polylines) == 2
	assert [len(element.get("points").split()) for element in polylines] == [3, 2]


def test_width_scaling_preserves_aspect():
	''' Doubling the width doubles the viewBox height too. '''
	first = ElementTree.fromstring(render_svg(LINE_DOCUMENT, width=640)).get("viewBox").split()
	second = ElementTree.fromstring(render_svg(LINE_DOCUMENT, width=1280)).get("viewBox").split()
	assert float(second[2]) == 2 * float(first[2])
	assert float(second[3]) == pytest.approx(2 * float(first[3]))


def test_non_square_aspect_is_honoured():
	canvas = {"gbox": [0, 0, 1, 0.5], "aspect": 0.5, "margin": 0.15}
	doc = document([{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.4]}], canvas=canvas)
	viewbox = ElementTree.fromstring(render_svg(doc)).get("viewBox").split()
	# (Gh + 2m) / (Gw + 2m) with Gw=1, Gh=0.5, m=0.15
	assert float(viewbox[3]) / float(viewbox[2]) == pytest.approx((0.5 + 0.3) / (1 + 0.3))


def test_render_svg_rejects_a_bare_envelope():
	'''
	§4.5: a DOCUMENT dict only. A caller holding envelopes assembles a document
	first — SkyPlot.document() exists precisely for that.
	'''
	from cornish.exc import EnvelopeError
	bare = {"cornish_view": "1.0", "type": "draw_stream",
	        "canvas": {"gbox": [0, 0, 1, 1]},
	        "layer": {"id": "grid", "action": "replace"}, "commands": []}
	with pytest.raises(EnvelopeError):
		render_svg(bare)


def test_render_svg_rejects_a_bad_width():
	with pytest.raises(ValueError):
		render_svg(LINE_DOCUMENT, width=0)
	with pytest.raises(TypeError):
		render_svg(LINE_DOCUMENT, width="640")


def test_background_is_painted_when_given_and_absent_when_null():
	with_background = document([], canvas={"gbox": [0, 0, 1, 1], "background": "#ffffff"})
	assert '<rect x="0" y="0"' in render_svg(with_background)
	assert '<rect x="0" y="0"' not in render_svg(document([], canvas={"gbox": [0, 0, 1, 1]}))


def test_text_is_escaped():
	doc = document([{"op": "text", "text": "a<b & c", "x": 0.5, "y": 0.5,
	                 "just": "CC", "up": [0, 1]}])
	svg = render_svg(doc)
	assert "&lt;" in svg and "&amp;" in svg
	assert ElementTree.fromstring(svg).find(f".//{SVG_NS}text").text == "a<b & c"


# ---------------------------------------------------------------------------
# the frozen mapping tables (§4.2)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("style, dashed", [(1, False), (2, True), (3, True), (4, True), (5, True)])
def test_line_styles_map_to_dasharrays(style, dashed):
	doc = document([{"op": "attr", "prim": "line", "style": style},
	                {"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}])
	polyline = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}polyline")
	assert (polyline.get("stroke-dasharray") is not None) == dashed


def test_dash_lengths_scale_with_stroke_width():
	''' §4.2: dash lengths are multiples of the stroke width. '''
	def dasharray(width):
		doc = document([{"op": "attr", "prim": "line", "style": 2, "width": width},
		                {"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]}])
		return ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}polyline").get("stroke-dasharray")
	thin = [float(value) for value in dasharray(1).split(",")]
	thick = [float(value) for value in dasharray(4).split(",")]
	assert thick[0] == pytest.approx(4 * thin[0])


@pytest.mark.parametrize("font, weight, style", [
	(1, "normal", "normal"), (2, "bold", "normal"), (3, "normal", "italic"),
	(4, "bold", "italic"), (5, "500", "normal"), (6, "600", "normal"),
])
def test_font_indices_map_to_css(font, weight, style):
	doc = document([{"op": "attr", "prim": "text", "font": font},
	                {"op": "text", "text": "x", "x": 0.5, "y": 0.5, "just": "CC", "up": [0, 1]}])
	element = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}text")
	assert element.get("font-weight") == weight
	assert element.get("font-style") == style


@pytest.mark.parametrize("just, anchor, baseline", [
	("TL", "start", "hanging"), ("CC", "middle", "middle"),
	("BR", "end", "alphabetic"), ("ML", "start", "text-after-edge"),
])
def test_justification_maps_to_svg(just, anchor, baseline):
	doc = document([{"op": "text", "text": "x", "x": 0.5, "y": 0.5, "just": just, "up": [0, 1]}])
	element = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}text")
	assert element.get("text-anchor") == anchor
	assert element.get("dominant-baseline") == baseline


def test_horizontal_text_is_not_rotated_and_vertical_text_is():
	horizontal = document([{"op": "text", "text": "x", "x": 0.5, "y": 0.5, "just": "CC", "up": [0, 1]}])
	assert ElementTree.fromstring(render_svg(horizontal)).find(f".//{SVG_NS}text").get("transform") is None
	vertical = document([{"op": "text", "text": "x", "x": 0.5, "y": 0.5, "just": "CC", "up": [-1, 0]}])
	transform = ElementTree.fromstring(render_svg(vertical)).find(f".//{SVG_NS}text").get("transform")
	assert transform is not None and transform.startswith("rotate(-90")


def test_text_carries_textlength_from_the_metrics_table():
	''' §4.5: pins browser text width to the layout AST computed. '''
	doc = document([{"op": "text", "text": "2:00", "x": 0.5, "y": 0.5, "just": "CC", "up": [0, 1]}])
	element = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}text")
	assert float(element.get("textLength")) > 0


def test_text_em_size_follows_the_single_scale_rule():
	'''
	§4.2: every size is a gbox quantity times the ONE scale s — never the raw
	canvas width. At the defaults that is ~9.8 px text.
	'''
	doc = document([{"op": "text", "text": "x", "x": 0.5, "y": 0.5, "just": "CC", "up": [0, 1]}])
	element = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}text")
	renderer = SVGRenderer(doc)
	# the output is rounded to 3 decimal places for byte-stability, so compare
	# within that quantum rather than exactly
	assert float(element.get("font-size")) == pytest.approx(0.02 * renderer._scale, abs=0.001)
	assert float(element.get("font-size")) == pytest.approx(9.85, abs=0.05)


@pytest.mark.parametrize("code, tag", [
	(1, "circle"), (2, "path"), (3, "path"), (4, "circle"), (5, "path"),
	(6, "circle"), (7, "polygon"), (8, "polygon"), (9, "polygon"), (10, "polygon"),
])
def test_every_frozen_marker_code_has_geometry(code, tag):
	doc = document([{"op": "mark", "x": [0.5], "y": [0.5], "marker": code}])
	root = ElementTree.fromstring(render_svg(doc))
	assert root.find(f".//{SVG_NS}{tag}") is not None


def test_dot_marker_is_filled_and_circle_marker_is_stroked():
	def marker(code):
		doc = document([{"op": "mark", "x": [0.5], "y": [0.5], "marker": code}])
		return ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}circle")
	assert marker(6).get("fill") != "none"
	assert marker(4).get("fill") == "none"


def test_clipping_policy():
	'''
	§4.4: grid lines, marks and images clip to the gbox; text is unclipped,
	because AST deliberately places labels outside it.
	'''
	doc = document([
		{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]},
		{"op": "mark", "x": [0.5], "y": [0.5], "marker": 4},
		{"op": "text", "text": "x", "x": 0.5, "y": -0.1, "just": "CC", "up": [0, 1]},
	])
	root = ElementTree.fromstring(render_svg(doc))
	assert root.find(f".//{SVG_NS}polyline").get("clip-path") is not None
	assert root.find(f".//{SVG_NS}text").get("clip-path") is None


# ---------------------------------------------------------------------------
# the frozen fallbacks (§3.3) — behaviours, not rejections
# ---------------------------------------------------------------------------

FALLBACK_DOCUMENT = document([
	{"op": "attr", "prim": "line", "colour": 999, "style": 9},
	{"op": "line", "x": [0.1, 0.9], "y": [0.1, 0.9]},
	{"op": "attr", "prim": "text", "font": 42},
	{"op": "text", "text": "x", "x": 0.5, "y": 0.5, "just": "CC", "up": [0, 1]},
	{"op": "mark", "x": [0.5], "y": [0.5], "marker": 99},
	{"op": "hologram", "intensity": 3},
])


def test_fallback_fixture_still_renders_everything():
	'''
	Out-of-range values degrade VISIBLY — foreground / solid / normal / dot —
	and the unknown op is skipped. Nothing is dropped without saying so.
	'''
	root = ElementTree.fromstring(render_svg(FALLBACK_DOCUMENT))
	assert root.find(f".//{SVG_NS}polyline") is not None
	assert root.find(f".//{SVG_NS}text") is not None
	assert root.find(f".//{SVG_NS}circle") is not None      # marker 99 -> DOT
	polyline = root.find(f".//{SVG_NS}polyline")
	assert polyline.get("stroke") == "#000000"              # colour 999 -> foreground
	assert polyline.get("stroke-dasharray") is None         # style 9 -> solid


@pytest.mark.parametrize("fragment", [
	"colour index 999", "line style 9", "font index 42", "marker code 99", "unknown command op",
])
def test_every_fallback_emits_both_diagnostic_channels(fragment, caplog):
	'''
	§3.7 defines "with a diagnostic" ONCE for every fallback, as two channels:
	an SVG comment and a logged warning. Independent consumers must not invent
	incompatible channels, and a silent fallback is indistinguishable from a bug.
	'''
	with caplog.at_level(logging.WARNING, logger="cornish"):
		svg = render_svg(FALLBACK_DOCUMENT)
	assert f"<!-- {fragment}" in svg or fragment in svg
	assert any(fragment in record.getMessage() for record in caplog.records), \
		f"no logged warning mentioned {fragment!r}"


def test_comment_channel_cannot_break_the_xml():
	''' A diagnostic must never produce an invalid comment. '''
	ElementTree.fromstring(render_svg(FALLBACK_DOCUMENT))


# ---------------------------------------------------------------------------
# the image op and its PNG writer (§2.4, §3.3)
# ---------------------------------------------------------------------------

def decode_png(data_uri: str):
	''' Decode an 8-bit greyscale PNG data URI back to rows, with stdlib zlib. '''
	raw = base64.b64decode(data_uri.split(",", 1)[1])
	assert raw[:8] == b"\x89PNG\r\n\x1a\n"
	width, height, depth, colour_type = struct.unpack(">IIBB", raw[16:26])
	assert (depth, colour_type) == (8, 0)
	start = raw.index(b"IDAT") + 4
	length = struct.unpack(">I", raw[start - 8:start - 4])[0]
	scanlines = zlib.decompress(raw[start:start + length])
	stride = width + 1
	return [list(scanlines[row * stride + 1:(row + 1) * stride]) for row in range(height)]


def test_image_data_uri_round_trips_and_is_flipped_vertically():
	'''
	§2.4 orientation rule, pinned with a vertically-asymmetric array: the input
	follows the FITS/bbox convention (row 0 = BOTTOM), and the PNG writer flips
	so PNG row 0 is the TOP.
	'''
	numpy = pytest.importorskip("numpy")
	data = numpy.array([[0, 0, 0, 0],
	                    [0, 0, 0, 0],
	                    [0, 0, 0, 0],
	                    [255, 255, 255, 255]], dtype=float)   # bright row is the TOP of the field
	rows = decode_png(array_to_png_data_uri(data))
	assert rows[0] == [255, 255, 255, 255]
	assert rows[-1] == [0, 0, 0, 0]


def test_image_scaling_uses_vmin_vmax():
	numpy = pytest.importorskip("numpy")
	data = numpy.array([[0.0, 1.0], [2.0, 3.0]])
	rows = decode_png(array_to_png_data_uri(data, vmin=0.0, vmax=3.0))
	assert rows[-1] == [0, 85]          # bottom row of the field is input row 0


def test_constant_image_maps_to_mid_scale():
	'''
	§2.4: vmin == vmax maps every pixel to 128 — a defined render, not a
	divide-by-zero and not a crash. A uniform image is legitimate data.
	'''
	numpy = pytest.importorskip("numpy")
	rows = decode_png(array_to_png_data_uri(numpy.full((2, 2), 7.0)))
	assert rows == [[128, 128], [128, 128]]


def test_reversed_colour_map():
	numpy = pytest.importorskip("numpy")
	data = numpy.array([[0.0, 1.0]])
	assert decode_png(array_to_png_data_uri(data, cmap="greys"))[0] == [0, 255]
	assert decode_png(array_to_png_data_uri(data, cmap="greys_r"))[0] == [255, 0]


def test_non_finite_pixels_are_reported_not_silently_substituted(caplog):
	'''
	NaNs are routine in astronomical data, so they are tolerated — but the
	substitution is announced, because an 8-bit greyscale PNG has no alpha to
	mark them blank and a silent black pixel is indistinguishable from data.
	'''
	numpy = pytest.importorskip("numpy")
	data = numpy.array([[1.0, float("nan")], [2.0, 3.0]])
	with caplog.at_level(logging.WARNING, logger="cornish"):
		array_to_png_data_uri(data)
	assert any("non-finite" in record.getMessage() for record in caplog.records)


@pytest.mark.parametrize("factory, match", [
	(lambda numpy: numpy.array([[float("nan")]]), "no finite values"),
	(lambda numpy: numpy.zeros((0, 0)), "empty"),
	(lambda numpy: numpy.zeros((2, 2, 2)), "2-D"),
])
def test_image_data_failures_raise(factory, match):
	numpy = pytest.importorskip("numpy")
	with pytest.raises(ValueError, match=match):
		array_to_png_data_uri(factory(numpy))


def test_unknown_colour_map_raises():
	numpy = pytest.importorskip("numpy")
	with pytest.raises(ValueError, match="stage 0"):
		array_to_png_data_uri(numpy.zeros((2, 2)), cmap="viridis")


def test_png_encoder_refuses_ragged_input():
	with pytest.raises(ValueError, match="same length"):
		encode_greyscale_png([b"\x00\x01", b"\x00"])
	with pytest.raises(ValueError):
		encode_greyscale_png([])


def test_image_op_renders_an_image_element():
	doc = document([{"op": "image", "x1": 0.0, "y1": 0.0, "x2": 1.0, "y2": 1.0,
	                 "href": "data:image/png;base64,iVBORw0KGgo="}])
	element = ElementTree.fromstring(render_svg(doc)).find(f".//{SVG_NS}image")
	assert element is not None
	assert element.get("clip-path") is not None
	assert float(element.get("width")) > 0
