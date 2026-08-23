
'''
The Python SVG writer (SPEC-09A §4.2–4.5) and the stdlib PNG encoder (§2.4).

``render_svg`` is a **pure function**: no I/O, no globals, deterministic output
(element order = render order, attribute order fixed), so ``render_svg(doc)``
is byte-stable for a given input. The Flask/daemon ETag pattern depends on
that.

It implements the same mapping tables as ``cornish-view.js`` — one table, two
implementations, cross-checked by the shared conformance fixtures (§7.6). The
renderer does no layout: AST already placed every label server-side, and this
just replays the frozen geometry.

Fallbacks (§3.3) are *behaviours*, not rejections: an out-of-range colour
renders as foreground, an unknown line style as solid, an unknown font as
normal, an unknown marker as a dot, and an unknown op is skipped. Every one of
them emits the §3.7 two-channel diagnostic — an SVG comment *and* a logged
warning — so a degraded render is always visible and never silent.
'''

import logging
import math
import struct
import zlib
from typing import Optional

from .envelope import ATTR_DEFAULTS, GrfOp, PrimType, validate_document
from .metrics import FontMetrics
from .palette import DEFAULT_PALETTE

__all__ = ['SVGRenderer', 'render_svg', 'encode_greyscale_png', 'array_to_png_data_uri']

logger = logging.getLogger("cornish")

#: Line style -> SVG stroke-dasharray pattern, in multiples of the stroke width
#: (§4.2; PGPLOT semantics, matching astplot-js's pixel values at width 1).
DASH_PATTERNS = {
	1: (),                                  # solid
	2: (8, 4),                              # dashed
	3: (8, 4, 2, 4),                        # dash-dot
	4: (2, 4),                              # dotted
	5: (8, 4, 2, 4, 2, 4),                  # dash-dot-dot-dot
}

#: Font index -> (CSS font-weight, CSS font-style) (§3.3, §4.2).
FONT_STYLES = {
	1: ("normal", "normal"),
	2: ("bold", "normal"),
	3: ("normal", "italic"),
	4: ("bold", "italic"),
	5: ("500", "normal"),
	6: ("600", "normal"),
}

#: Justification -> SVG presentation attributes (§4.2).
BASELINES = {"T": "hanging", "C": "middle", "B": "alphabetic", "M": "text-after-edge"}
ANCHORS = {"L": "start", "C": "middle", "R": "end"}

#: Nominal marker radius and stroke width, in gbox units per unit attr (§4.2).
MARKER_RADIUS_UNIT = 0.006
STROKE_WIDTH_UNIT = 0.002
#: Below this a hairline stops being visible at all.
MINIMUM_STROKE_PX = 0.5

CLIP_PATH_ID = "cornish-gbox"

#: Decimal places for coordinates in the output. Enough that a 4096 px render
#: is sub-milli-pixel accurate, few enough to keep the bytes stable.
_PLACES = 3


# ---------------------------------------------------------------------------
# formatting helpers (deterministic by construction)
# ---------------------------------------------------------------------------

def _number(value: float) -> str:
	''' Format a coordinate deterministically, without a trailing ".0". '''
	text = f"{float(value):.{_PLACES}f}"
	if "." in text:
		text = text.rstrip("0").rstrip(".")
	return text or "0"


def _escape_text(text: str) -> str:
	return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _escape_attribute(text: str) -> str:
	return _escape_text(str(text)).replace('"', "&quot;")


def _comment(text: str) -> str:
	'''
	An XML comment carrying a diagnostic. ``--`` is illegal inside a comment, so
	it is spaced out rather than dropped — the message stays readable.
	'''
	return f"<!-- {_escape_text(str(text)).replace('--', '- -')} -->"


# ---------------------------------------------------------------------------
# PNG encoding (stdlib only)
# ---------------------------------------------------------------------------

def _chunk(kind: bytes, payload: bytes) -> bytes:
	return (struct.pack(">I", len(payload)) + kind + payload
	        + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF))


def encode_greyscale_png(rows) -> bytes:
	'''
	Encode 8-bit greyscale *rows* (a sequence of equal-length ``bytes``) as a
	PNG, using only ``zlib`` and ``struct``.

	Row 0 of *rows* becomes the TOP row of the image, per the PNG convention.

	:raises ValueError: if there are no rows, or the rows are ragged
	'''
	rows = [bytes(row) for row in rows]
	if not rows:
		raise ValueError("A PNG needs at least one row of pixels.")
	width = len(rows[0])
	if width == 0:
		raise ValueError("A PNG needs at least one pixel per row.")
	for index, row in enumerate(rows):
		if len(row) != width:
			raise ValueError(f"PNG rows must all be the same length; row {index} has "
			                 f"{len(row)} bytes, row 0 has {width}.")

	header = struct.pack(">IIBBBBB", width, len(rows), 8, 0, 0, 0, 0)
	# filter type 0 (None) on every scanline
	raw = b"".join(b"\x00" + row for row in rows)
	return (b"\x89PNG\r\n\x1a\n"
	        + _chunk(b"IHDR", header)
	        + _chunk(b"IDAT", zlib.compress(raw, 9))
	        + _chunk(b"IEND", b""))


def array_to_png_data_uri(data, *, vmin: Optional[float] = None,
                          vmax: Optional[float] = None, cmap: str = "greys") -> str:
	'''
	Scale a 2-D array to 8-bit greyscale and return it as a ``data:`` URI (§2.4).

	Orientation: the input follows the FITS/bbox convention (**row 0 is the
	BOTTOM row**, lowest y); this function performs the vertical flip so that
	PNG row 0 is the top, as the §3.3 image op requires.

	Degenerate scaling is defined, not an error: ``vmin == vmax`` (a constant
	image, or explicit equal bounds) maps every pixel to mid-scale 128. A
	uniform image is legitimate data — nothing failed.

	Non-finite pixels (NaN/Inf, routine in astronomical data) are excluded from
	the automatic range and rendered at the low end, **with a logged diagnostic
	naming how many** — an explicit, auditable decision rather than a silent
	substitution. An 8-bit greyscale PNG has no alpha channel to express
	"blank" with; a masked rendering is stage-3 work.

	:param cmap: ``'greys'`` (low values dark) or ``'greys_r'`` (low values light)
	:raises ValueError: if the array is not 2-D, is empty, is entirely non-finite, or *cmap* is unknown
	'''
	import numpy as np

	if cmap not in ("greys", "greys_r"):
		raise ValueError(f"cmap '{cmap}' is not supported in stage 0 "
		                 "(only 'greys' and 'greys_r'; a colour LUT is stage-3 work).")

	values = np.asarray(data, dtype=float)
	if values.ndim != 2:
		raise ValueError(f"Image data must be a 2-D array; got {values.ndim} dimension(s).")
	if values.size == 0:
		raise ValueError("Image data is empty.")

	finite = np.isfinite(values)
	if not finite.any():
		raise ValueError("Image data contains no finite values, so it cannot be scaled.")
	non_finite_count = int((~finite).sum())
	if non_finite_count:
		logger.warning("addImage: %d of %d pixels are non-finite; they are excluded from the "
		               "automatic scaling range and rendered at the low end (an 8-bit greyscale "
		               "PNG has no alpha channel to mark them blank).",
		               non_finite_count, values.size)

	low = float(np.min(values[finite])) if vmin is None else float(vmin)
	high = float(np.max(values[finite])) if vmax is None else float(vmax)
	if not (math.isfinite(low) and math.isfinite(high)):
		raise ValueError(f"vmin/vmax must be finite; got vmin={low!r}, vmax={high!r}.")

	if low == high:
		scaled = np.full(values.shape, 128, dtype=np.uint8)
	else:
		if high < low:
			raise ValueError(f"vmax ({high}) must not be less than vmin ({low}).")
		normalized = (values - low) / (high - low)
		normalized = np.where(finite, normalized, 0.0)
		scaled = np.clip(np.rint(normalized * 255.0), 0, 255).astype(np.uint8)

	if cmap == "greys_r":
		scaled = 255 - scaled

	# input row 0 is the BOTTOM row; PNG row 0 is the TOP row
	rows = [scaled[index].tobytes() for index in range(scaled.shape[0] - 1, -1, -1)]

	import base64
	encoded = base64.b64encode(encode_greyscale_png(rows)).decode("ascii")
	return f"data:image/png;base64,{encoded}"


# ---------------------------------------------------------------------------
# the renderer
# ---------------------------------------------------------------------------

class SVGRenderer:
	'''
	Renders a draw-stream document to an SVG string (§4.2–4.5).

	:param document: a full document (the ``layers`` form) — never a bare envelope
	:param width: output pixel width
	:param indent: pretty-print the element tree (does not change the geometry)
	:param metrics: font metrics for ``textLength``; ``None`` loads the shipped table
	:raises EnvelopeError: if *document* is not a valid document
	:raises ValueError: if *width* is not a positive integer
	'''

	def __init__(self, document: dict, *, width: int = 640, indent: bool = False,
	             metrics: Optional[FontMetrics] = None) -> None:
		validate_document(document)
		if isinstance(width, bool) or not isinstance(width, int):
			raise TypeError(f"width must be an integer number of pixels; got {type(width).__name__}.")
		if width <= 0:
			raise ValueError(f"width must be a positive integer; got {width}.")

		self._document = document
		self._width = width
		self._indent = indent
		self._metrics = metrics if metrics is not None else FontMetrics()

		canvas = document["canvas"]
		self._canvas = canvas
		gbox = list(canvas.get("gbox", [0.0, 0.0, 1.0, 1.0]))
		self._x1, self._y1, self._x2, self._y2 = (float(value) for value in gbox)
		gbox_width = self._x2 - self._x1
		gbox_height = self._y2 - self._y1
		margin = float(canvas.get("margin", 0.15))

		# §4.3 — margins are equal in ABSOLUTE units on all four sides
		self._margin = margin * max(gbox_width, gbox_height)
		self._scale = self._width / (gbox_width + 2 * self._margin)
		self._height = self._scale * (gbox_height + 2 * self._margin)

		self._char_height = float(canvas.get("char_height", 0.02))
		self._font_family = str(canvas.get("font_family", "DejaVu Sans"))
		self._palette = list(canvas.get("palette", DEFAULT_PALETTE))
		self._background = canvas.get("background", None)

	# -- coordinate mapping (§4.3) -------------------------------------

	def px(self, x: float) -> float:
		return self._scale * (float(x) - self._x1 + self._margin)

	def py(self, y: float) -> float:
		''' Note the flip: envelope coordinates are Y-up, SVG pixels are Y-down. '''
		return self._height - self._scale * (float(y) - self._y1 + self._margin)

	# -- attribute resolution (§3.3 fallbacks) -------------------------

	def _colour(self, value, diagnostics: list) -> str:
		'''
		Resolve an attr colour to a CSS colour.

		An out-of-range palette index renders as foreground **with a
		diagnostic** — this covers both a malformed stream and the legal but
		hazardous case of a canvas re-theme replacing the palette with a shorter
		one, stranding stored integer references. The primitive still renders,
		visibly; it is never skipped.
		'''
		if isinstance(value, str):
			return value
		index = int(value)
		if 0 <= index < len(self._palette):
			return str(self._palette[index])
		message = (f"colour index {index} is outside the {len(self._palette)}-entry palette; "
		           "rendering as foreground (index 1)")
		diagnostics.append(message)
		logger.warning("render_svg: %s", message)
		return str(self._palette[1]) if len(self._palette) > 1 else "#000000"

	def _dash(self, style, stroke_px: float, diagnostics: list) -> str:
		pattern = DASH_PATTERNS.get(int(style))
		if pattern is None:
			message = f"line style {style} is outside the frozen 1-5 table; rendering solid"
			diagnostics.append(message)
			logger.warning("render_svg: %s", message)
			pattern = DASH_PATTERNS[1]
		if not pattern:
			return ""
		return ",".join(_number(value * stroke_px) for value in pattern)

	def _font(self, font, diagnostics: list):
		style = FONT_STYLES.get(int(font))
		if style is None:
			message = f"font index {font} is outside the frozen 1-6 table; rendering as normal"
			diagnostics.append(message)
			logger.warning("render_svg: %s", message)
			style = FONT_STYLES[1]
		return style

	def _stroke_px(self, width_attr: float) -> float:
		return max(STROKE_WIDTH_UNIT * float(width_attr) * self._scale, MINIMUM_STROKE_PX)

	# -- markers (§3.3 geometry table) ---------------------------------

	def _marker_elements(self, code, centre_x: float, centre_y: float, radius: float,
	                     stroke_attributes: str, colour: str, diagnostics: list) -> list:
		'''
		The frozen 10-code marker geometry. An unknown code renders as a DOT
		with a diagnostic — a visible fallback, never a skip.
		'''
		code = int(code)
		if code not in range(1, 11):
			message = f"marker code {code} is outside the frozen 1-10 table; rendering as a dot"
			diagnostics.append(message)
			logger.warning("render_svg: %s", message)
			code = 6

		if code == 6:                                   # DOT — filled
			return [f'<circle cx="{_number(centre_x)}" cy="{_number(centre_y)}" '
			        f'r="{_number(0.25 * radius)}" fill="{_escape_attribute(colour)}"/>']
		if code in (1, 4):                              # circles — stroked
			factor = 0.6 if code == 1 else 1.0
			return [f'<circle cx="{_number(centre_x)}" cy="{_number(centre_y)}" '
			        f'r="{_number(factor * radius)}" fill="none" {stroke_attributes}/>']

		#: half-length-R segments through the centre, by bearing in degrees
		segments = {
			2: (0.0, 90.0),                             # CROSS '+'
			5: (45.0, 135.0),                           # X
			3: (0.0, 60.0, 120.0),                      # STAR
		}.get(code)

		if segments:
			path = []
			for angle_degrees in segments:
				angle = angle_degrees * math.pi / 180.0
				delta_x = radius * math.cos(angle)
				delta_y = radius * math.sin(angle)
				path.append(f"M{_number(centre_x - delta_x)},{_number(centre_y - delta_y)}"
				            f"L{_number(centre_x + delta_x)},{_number(centre_y + delta_y)}")
			return [f'<path d="{"".join(path)}" fill="none" {stroke_attributes}/>']

		# triangles 7-10: equilateral, circumradius R, apex per the code.
		# Angles are in SVG's Y-DOWN frame, so "apex up" is -y.
		apex = {7: -90.0, 8: 90.0, 9: 180.0, 10: 0.0}[code]
		points = []
		for offset in (0.0, 120.0, 240.0):
			angle = (apex + offset) * math.pi / 180.0
			points.append(f"{_number(centre_x + radius * math.cos(angle))},"
			              f"{_number(centre_y + radius * math.sin(angle))}")
		return [f'<polygon points="{" ".join(points)}" fill="none" {stroke_attributes}/>']

	# -- rendering -----------------------------------------------------

	def render(self) -> str:
		''' The SVG document as a string. '''
		lines = []
		lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
		             f'viewBox="0 0 {_number(self._width)} {_number(self._height)}" '
		             f'width="100%" preserveAspectRatio="xMidYMid meet">')

		# grid lines, marks and images clip to the gbox; text is deliberately
		# unclipped, because AST places labels OUTSIDE it (that is what margin
		# is for) — matching astplot-js's bridge (§4.4)
		lines.append(f'<defs><clipPath id="{CLIP_PATH_ID}">'
		             f'<rect x="{_number(self.px(self._x1))}" y="{_number(self.py(self._y2))}" '
		             f'width="{_number(self._scale * (self._x2 - self._x1))}" '
		             f'height="{_number(self._scale * (self._y2 - self._y1))}"/>'
		             f'</clipPath></defs>')

		if self._background is not None:
			lines.append(f'<rect x="0" y="0" width="{_number(self._width)}" '
			             f'height="{_number(self._height)}" '
			             f'fill="{_escape_attribute(self._background)}"/>')

		for layer in self._ordered_layers():
			lines.extend(self._render_layer(layer))

		lines.append("</svg>")
		separator = "\n" if self._indent else ""
		return separator.join(lines)

	def _ordered_layers(self) -> list:
		'''
		§3.5: layers render sorted by (z, creation sequence). The array position
		is the creation sequence, and Python's sort is stable, so this is
		deterministic for any input ordering.
		'''
		return sorted(self._document["layers"], key=lambda layer: int(layer.get("z", 0)))

	def _render_layer(self, layer: dict) -> list:
		'''
		One ``<g data-layer-id=...>`` group. There is deliberately **no ``id``
		attribute** (§3.5): the legal id charset contains ``.`` and ``-``, which
		break ``#...`` CSS selectors — lookup is by attribute selector instead.
		'''
		attributes = f'data-layer-id="{_escape_attribute(layer["id"])}"'
		if "kind" in layer:
			# a layer with no kind renders with the attribute omitted entirely,
			# never as data-kind=""
			attributes += f' data-kind="{_escape_attribute(layer["kind"])}"'
		lines = [f"<g {attributes}>"]

		state = {PrimType.LINE.value: dict(ATTR_DEFAULTS),
		         PrimType.MARK.value: dict(ATTR_DEFAULTS),
		         PrimType.TEXT.value: dict(ATTR_DEFAULTS)}

		for command in layer["commands"]:
			lines.extend(self._render_command(command, state))

		lines.append("</g>")
		return lines

	def _render_command(self, command: dict, state: dict) -> list:
		diagnostics = []
		operation = command.get("op")
		lines = []

		if operation == GrfOp.ATTR.value:
			# attr state is per-layer and resets at every layer head (§3.4)
			primitive = command["prim"]
			for key in ("colour", "width", "style", "size", "font"):
				if key in command:
					state[primitive][key] = command[key]
			return lines

		if operation == GrfOp.LINE.value:
			lines.extend(self._render_line(command, state[PrimType.LINE.value], diagnostics))
		elif operation == GrfOp.MARK.value:
			lines.extend(self._render_mark(command, state[PrimType.MARK.value], diagnostics))
		elif operation == GrfOp.TEXT.value:
			lines.extend(self._render_text(command, state[PrimType.TEXT.value], diagnostics))
		elif operation == GrfOp.IMAGE.value:
			lines.extend(self._render_image(command))
		else:
			# §3.7: an unknown op is SKIPPED WITH A DIAGNOSTIC, never silently
			# dropped and never a hard failure — that is how a 1.0 consumer
			# survives a 1.1 producer
			message = f"unknown command op '{operation}' skipped (this consumer implements schema 1.0)"
			diagnostics.append(message)
			logger.warning("render_svg: %s", message)

		return [_comment(message) for message in diagnostics] + lines

	def _render_line(self, command: dict, state: dict, diagnostics: list) -> list:
		colour = self._colour(state["colour"], diagnostics)
		stroke_px = self._stroke_px(state["width"])
		dash = self._dash(state["style"], stroke_px, diagnostics)
		points = " ".join(f"{_number(self.px(x))},{_number(self.py(y))}"
		                  for x, y in zip(command["x"], command["y"]))
		dash_attribute = f' stroke-dasharray="{dash}"' if dash else ""
		return [f'<polyline points="{points}" fill="none" stroke="{_escape_attribute(colour)}" '
		        f'stroke-width="{_number(stroke_px)}"{dash_attribute} '
		        f'clip-path="url(#{CLIP_PATH_ID})"/>']

	def _render_mark(self, command: dict, state: dict, diagnostics: list) -> list:
		colour = self._colour(state["colour"], diagnostics)
		stroke_px = self._stroke_px(state["width"])
		radius = MARKER_RADIUS_UNIT * float(state["size"]) * self._scale
		stroke_attributes = (f'stroke="{_escape_attribute(colour)}" '
		                     f'stroke-width="{_number(stroke_px)}"')
		elements = []
		for x, y in zip(command["x"], command["y"]):
			elements.extend(self._marker_elements(command["marker"], self.px(x), self.py(y),
			                                      radius, stroke_attributes, colour, diagnostics))
		if not elements:
			return []
		return [f'<g clip-path="url(#{CLIP_PATH_ID})">'] + elements + ["</g>"]

	def _render_text(self, command: dict, state: dict, diagnostics: list) -> list:
		colour = self._colour(state["colour"], diagnostics)
		weight, style = self._font(state["font"], diagnostics)
		em = self._char_height * float(state["size"])
		font_size = em * self._scale

		anchor_x = self.px(command["x"])
		anchor_y = self.py(command["y"])
		just = command["just"]
		up_x, up_y = command["up"]

		# §4.2: in SVG's Y-down frame the rotation is atan2(upx, upy). The
		# 180/pi factor is written out rather than using math.degrees so this
		# stays clear of the package's coordinate-conversion gate.
		angle = math.atan2(float(up_x), float(up_y)) * 180.0 / math.pi
		transform = ""
		if abs(angle) > 1e-9:
			transform = (f' transform="rotate({_number(angle)},'
			             f'{_number(anchor_x)},{_number(anchor_y)})"')

		text = str(command["text"])
		# textLength pins the browser's text width to the layout AST computed,
		# which is what keeps label spacing consistent across viewer fonts
		length = self._metrics.string_width(text, font_size) if text else 0.0
		length_attribute = f' textLength="{_number(length)}"' if length > 0 else ""

		return [f'<text x="{_number(anchor_x)}" y="{_number(anchor_y)}" '
		        f'fill="{_escape_attribute(colour)}" '
		        f'font-family="{_escape_attribute(self._font_family)}, sans-serif" '
		        f'font-size="{_number(font_size)}" font-weight="{weight}" font-style="{style}" '
		        f'text-anchor="{ANCHORS[just[1]]}" dominant-baseline="{BASELINES[just[0]]}"'
		        f'{length_attribute}{transform}>{_escape_text(text)}</text>']

	def _render_image(self, command: dict) -> list:
		left = self.px(min(command["x1"], command["x2"]))
		right = self.px(max(command["x1"], command["x2"]))
		top = self.py(max(command["y1"], command["y2"]))
		bottom = self.py(min(command["y1"], command["y2"]))
		return [f'<image x="{_number(left)}" y="{_number(top)}" '
		        f'width="{_number(right - left)}" height="{_number(bottom - top)}" '
		        f'preserveAspectRatio="none" clip-path="url(#{CLIP_PATH_ID})" '
		        f'href="{_escape_attribute(command["href"])}"/>']


def render_svg(document: dict, *, width: int = 640, indent: bool = False) -> str:
	'''
	Render a draw-stream **document** to an SVG string (§4.5).

	This takes a document dict, never a bare envelope or a list of them: a
	caller holding envelopes assembles a document first (``SkyPlot.document()``
	exists precisely for that).

	Pure and deterministic — the same input always produces the same bytes,
	which is what lets the daemon use ``envelope_hash`` as an ETag.

	:raises EnvelopeError: if *document* is not a valid document
	'''
	return SVGRenderer(document, width=width, indent=indent).render()
