
'''
``SkyPlot`` v2 — the frozen plotting API (SPEC-09A §2).

There is **no backend option** (D17). ``SkyPlot`` always drives ``Ast.Plot``
through a private :class:`~cornish.plot.serial.GrfSequence`; SVG, JSON, the
browser and matplotlib are all *consumers of the captured stream*. One Grf
implementation means one metrics source and one geometry truth, so the same
plot is identical (modulo font rendering) in a file, a notebook, a browser and
a matplotlib window.

``cornish/plot/matplotlib.py`` (SkyPlot v1) is frozen and superseded; nothing
here imports it.
'''

import logging
import math
import os
from typing import Optional, Sequence, Union

import astropy.units as u
import numpy as np
import starlink.Ast as Ast

from .. import _pyast_bridge as bridge
from ..enums import FontStyle, LineStyle, MARKER_INPUT_NAMES, Marker, SkySystem
from ..exc import IncompleteHeader
from ..mapping.frame.frame_set import ASTFrameSet
from ..region.circle import ASTCircle
from ..region.region import ASTRegion
from .cornish import CornishPlot
from .envelope import (ID_PATTERN, SCHEMA_VERSION, LayerAction, canonical_json,
                       make_envelope, validate_document)
from .metrics import FontMetrics
from .palette import DEFAULT_PALETTE, normalize_colour, normalize_palette_entry
from .serial import GrfSequence
from .svg import array_to_png_data_uri, render_svg

__all__ = ['SkyPlot']

logger = logging.getLogger("cornish")

#: Sky systems with a FITS CTYPE encoding, and the cards each one needs (§2.2).
#: Others (AZEL, GAPPT, …) have no FITS representation; the escape hatch is
#: building your own frameset and using Ast.Plot directly.
SYSTEM_CARDS = {
	SkySystem.ICRS: ("RA---", "DEC--", {"RADESYS": "ICRS"}),
	SkySystem.FK5: ("RA---", "DEC--", {"RADESYS": "FK5", "EQUINOX": 2000.0}),
	SkySystem.FK4: ("RA---", "DEC--", {"RADESYS": "FK4", "EQUINOX": 1950.0}),
	SkySystem.GALACTIC: ("GLON-", "GLAT-", {}),
	SkySystem.ECLIPTIC: ("ELON-", "ELAT-", {}),
}

#: Whole-sky projections verified this design cycle. A zenithal or conic code
#: cannot represent the full sky in the 2:1 equatorial parameterization.
ALL_SKY_PROJECTIONS = ("AIT", "MOL")

#: gbox aspect bounds (§2.2): outside this the 1/100 quantization error exceeds
#: 10%, and naxis2 < 5 is a degenerate shape.
MINIMUM_ASPECT = 0.05
MAXIMUM_ASPECT = 20.0

#: Default z (render order) by layer kind (§2.3).
DEFAULT_Z = {"image": -20, "grid": -10}

#: Layers whose id is fixed — calling the method twice replaces (§2.3).
FIXED_LAYER_IDS = {"grid": "grid", "compass": "compass", "scalebar": "scalebar"}

#: The §2.6 colour defaults, resolved in exactly one place.
REGION_COLOUR = "#4a7f7b"      # cornish's established teal
POINTS_COLOUR = 4              # palette blue, matching astplot-js Colour(Markers)
GRID_COLOUR = 11               # palette grey
FOREGROUND_COLOUR = 1


def _set_plot_attribute(ast_plot: Ast.Plot, name: str, value) -> None:
	'''
	Route a user value to an ``Ast.Plot`` attribute.

	Mirrors :meth:`cornish.ast_object.ASTObject._setAttribute` exactly:
	``astSet`` parses its argument as a COMMA-SEPARATED settings list, so the
	value is double-quoted (an unquoted value containing ``", Attr=Val"`` would
	silently truncate the value *and* set the other attribute), and an
	``Ast.AstError`` is chained into a ``ValueError`` per the type-vs-value
	policy. ``SkyPlot`` is not an ``ASTObject``, but the convention is the same.
	'''
	try:
		ast_plot.set(f'{name}="{value}"')
	except Ast.AstError as error:
		raise ValueError(f"AST could not set the plot attribute '{name}' to {value!r}.") from error


def _validate_common_parameters(*, char_height: float, margin: float, width: int,
                                palette: Optional[Sequence[str]], background: Optional[str]):
	'''
	The construction-time parameter checks shared by ``__init__`` and
	``fromFITSHeader`` (§2.1) — bypassing the synthetic-header builder must not
	bypass parameter validation.

	:returns: ``(char_height, margin, width, palette_or_None, background_or_None)``
	'''
	char_height = float(char_height)
	if not math.isfinite(char_height) or char_height <= 0:
		raise ValueError(f"char_height must be a positive finite number; got {char_height!r}.")

	margin = float(margin)
	if not math.isfinite(margin) or not 0 <= margin <= 0.5:
		raise ValueError(f"margin must lie in [0, 0.5]; got {margin!r}.")

	if isinstance(width, bool) or not isinstance(width, int):
		raise TypeError(f"width must be an integer number of pixels; got {type(width).__name__}.")
	if width <= 0:
		raise ValueError(f"width must be a positive integer; got {width}.")

	if palette is not None:
		# normalize_palette_entry, NOT normalize_colour: slot 0 is legitimately
		# 'rgba(0,0,0,0)'/'transparent', which the user-colour grammar rejects
		palette = [normalize_palette_entry(entry) for entry in palette]
		if len(palette) < 2:
			raise ValueError("A palette needs at least 2 entries "
			                 f"(0 = background, 1 = foreground); got {len(palette)}.")

	if background is not None:
		background = normalize_palette_entry(background)

	return char_height, margin, width, palette, background


class SkyPlot(CornishPlot):
	'''
	A sky plot that records its own drawing as a replayable command stream.

	Exactly one of *extent* / ``all_sky=True`` is required.

	:param extent: an :class:`~cornish.region.region.ASTRegion` on a sky frame covering the area to plot
	:param all_sky: plot the whole sky instead of a finite extent
	:param center: all-sky recentre only; a ``SkyCoord`` or ``(lon_deg, lat_deg)``
	:param projection: FITS projection code; ``None`` resolves to ``TAN`` (finite) or ``AIT`` (all-sky)
	:param system: display system; ``None`` inherits the extent's own
	:param rotation: sky rotation; a bare number means degrees
	:param aspect: gbox height/width; ``None`` means 1.0, or 0.5 all-sky
	:param width: advisory SVG pixel width
	:param char_height: height of size-1.0 text in gbox units
	:param margin: fraction of the larger gbox side reserved for labels
	:param palette: base palette; ``None`` uses :data:`~cornish.plot.palette.DEFAULT_PALETTE`
	:param background: canvas fill; ``None`` is transparent
	:param grid: draw the grid layer on construction
	:param plot_id: id used by the viewer; default ``plot-<pid>-<n>``
	'''

	#: process-wide counter behind the default plot_id
	_plot_counter = 0

	def __init__(self,
	             extent: Optional[ASTRegion] = None, *,
	             all_sky: bool = False,
	             center=None,
	             projection: Optional[str] = None,
	             system: Optional[Union[SkySystem, str]] = None,
	             rotation: Union[u.Quantity, float] = 0 * u.deg,
	             aspect: Optional[float] = None,
	             width: int = 640,
	             char_height: float = 0.02,
	             margin: float = 0.15,
	             palette: Optional[Sequence[str]] = None,
	             background: Optional[str] = None,
	             grid: bool = True,
	             plot_id: Optional[str] = None) -> None:

		if (extent is None) == (not all_sky):
			raise ValueError("Exactly one of 'extent' or 'all_sky=True' is required "
			                 f"(got extent={extent!r}, all_sky={all_sky!r}).")
		if center is not None and not all_sky:
			# with an extent, a centre would contradict the extent's own centre
			raise ValueError("'center' is only legal with all_sky=True; a finite extent "
			                 "already defines its own centre.")

		char_height, margin, width, palette, background = _validate_common_parameters(
			char_height=char_height, margin=margin, width=width,
			palette=palette, background=background)

		rotation_deg = self._rotation_degrees(rotation)

		if all_sky:
			cards, effective_aspect, naxis = self._all_sky_header(
				center=center, projection=projection, system=system,
				aspect=aspect, rotation_deg=rotation_deg)
		else:
			cards, effective_aspect, naxis = self._finite_header(
				extent=extent, projection=projection, system=system,
				aspect=aspect, rotation_deg=rotation_deg)

		self._cards = cards
		self._build(cards=cards, aspect=effective_aspect, naxis=naxis, width=width,
		            char_height=char_height, margin=margin, palette=palette,
		            background=background, plot_id=plot_id, grid=grid)

	# ------------------------------------------------------------------
	# construction helpers
	# ------------------------------------------------------------------

	@staticmethod
	def _rotation_degrees(rotation) -> float:
		''' A ``Quantity`` angle, or a bare number meaning degrees (§2.1). '''
		if isinstance(rotation, u.Quantity):
			try:
				return float(rotation.to_value(u.deg))
			except u.UnitConversionError as error:
				raise ValueError(f"'rotation' must be an angle; got {rotation!r}.") from error
		value = float(rotation)
		if not math.isfinite(value):
			raise ValueError(f"'rotation' must be finite; got {rotation!r}.")
		return value

	@staticmethod
	def _resolve_aspect(aspect: Optional[float]) -> float:
		if aspect is None:
			return 1.0
		aspect = float(aspect)
		if not math.isfinite(aspect):
			raise ValueError(f"'aspect' must be a finite number; got {aspect!r}.")
		if not MINIMUM_ASPECT <= aspect <= MAXIMUM_ASPECT:
			raise ValueError(f"'aspect' must lie in [{MINIMUM_ASPECT}, {MAXIMUM_ASPECT}]; "
			                 f"got {aspect}.")
		return aspect

	@staticmethod
	def _resolve_system(system, extent_system: Optional[str]) -> SkySystem:
		'''
		Resolve the display system. ``None`` inherits the extent's own; an
		extent on a system with no FITS encoding must be given one explicitly
		rather than being silently relabelled.
		'''
		if system is None:
			if extent_system is None:
				return SkySystem.ICRS
			try:
				candidate = SkySystem(extent_system)
			except ValueError:
				candidate = None
			if candidate is None or candidate not in SYSTEM_CARDS:
				raise ValueError(
					f"The extent is on the '{extent_system}' system, which has no FITS CTYPE "
					f"encoding, so a display system cannot be inherited. Pass system= "
					f"explicitly (supported: {[str(value) for value in SYSTEM_CARDS]}).")
			return candidate

		if isinstance(system, SkySystem):
			candidate = system
		else:
			try:
				candidate = SkySystem(str(system).upper())
			except ValueError as error:
				raise ValueError(f"'{system}' is not a known sky system "
				                 f"(supported here: {[str(value) for value in SYSTEM_CARDS]}).") from error
		if candidate not in SYSTEM_CARDS:
			raise ValueError(f"System '{candidate}' has no FITS CTYPE encoding and cannot be "
			                 f"plotted; supported: {[str(value) for value in SYSTEM_CARDS]}. "
			                 "Build your own frameset and use Ast.Plot directly if you need it.")
		return candidate

	def _finite_header(self, *, extent, projection, system, aspect, rotation_deg):
		''' The §2.2 finite-extent synthetic header. '''
		if not isinstance(extent, ASTRegion):
			raise TypeError("'extent' must be an ASTRegion covering the area to plot; "
			                f"got {type(extent).__name__}.")

		circle = extent if isinstance(extent, ASTCircle) else extent.boundingCircle()
		if not bridge.is_sky(circle.astObject):
			raise ValueError("SkyPlot requires an extent defined in a sky frame (got a region "
			                 f"on a '{bridge.current_frame_of(circle.astObject).Class}' frame).")

		centre = (float(circle.center[0]), float(circle.center[1]))
		radius_deg = float(circle.radius.to_value(u.deg))
		if not all(math.isfinite(value) for value in (*centre, radius_deg)):
			# v1 logged a warning here and proceeded to a guaranteed downstream
			# failure; that is silent-skip-adjacent, so v2 crashes at the source
			raise ValueError(f"The extent's bounding circle is not finite (centre={centre}, "
			                 f"radius={radius_deg}), so no plot can be built from it.")
		if radius_deg <= 0:
			raise ValueError(f"The extent's bounding circle has a non-positive radius ({radius_deg} deg).")

		extent_frame = bridge.current_frame_of(circle.astObject)
		extent_system = getattr(extent_frame, "System", None)
		display_system = self._resolve_system(system, extent_system)
		centre = self._convert_centre(centre, extent_frame, display_system)

		projection = (projection or "TAN").upper()
		if len(projection) != 3:
			raise ValueError(f"'projection' must be a 3-character FITS code; got '{projection}'.")

		aspect = self._resolve_aspect(aspect)
		naxis1 = 100
		naxis2 = int(round(naxis1 * aspect))
		# the effective aspect is quantized to 1/100 and REPORTED in the canvas:
		# Scales() = (1, 1) and Qch's isotropy claim are honest only if the gbox
		# aspect equals naxis2*CDELT2 / (naxis1*CDELT1) exactly
		effective_aspect = naxis2 / naxis1

		delta = 2.1 * radius_deg / min(naxis1, naxis2)
		longitude, latitude, extra = SYSTEM_CARDS[display_system]
		cards = {
			"NAXES": 2, "NAXIS1": naxis1, "NAXIS2": naxis2,
			"CRVAL1": centre[0], "CRVAL2": centre[1],
			"CTYPE1": f"{longitude}{projection}", "CTYPE2": f"{latitude}{projection}",
			"CRPIX1": naxis1 / 2 + 0.5, "CRPIX2": naxis2 / 2 + 0.5,
			"CDELT1": -delta,           # RA/longitude increases to the LEFT
			"CDELT2": delta,
		}
		cards.update(extra)
		if rotation_deg != 0:
			cards["CROTA2"] = rotation_deg
		return cards, effective_aspect, (naxis1, naxis2)

	def _all_sky_header(self, *, center, projection, system, aspect, rotation_deg):
		''' The §2.2 all-sky header: a fixed 200x100 grid, aspect exactly 0.5. '''
		projection = (projection or "AIT").upper()
		if projection not in ALL_SKY_PROJECTIONS:
			raise ValueError(f"'{projection}' cannot represent the whole sky in this "
			                 f"parameterization; all_sky supports {list(ALL_SKY_PROJECTIONS)}.")

		if aspect is not None and float(aspect) != 0.5:
			raise ValueError("With all_sky=True the field is inherently 2:1, so 'aspect' must be "
			                 f"None or exactly 0.5; got {aspect}.")

		display_system = self._resolve_system(system, None)
		centre_longitude = 0.0
		if center is not None:
			centre_longitude, centre_latitude = self._all_sky_centre(center, display_system)
			if abs(centre_latitude) > 1e-9:
				raise ValueError("An AIT/MOL fiducial point must lie on the equator, so the "
				                 f"all-sky 'center' latitude must be 0; got {centre_latitude} deg "
				                 f"in the {display_system} system.")

		naxis1, naxis2 = 200, 100
		longitude, latitude, extra = SYSTEM_CARDS[display_system]
		cards = {
			"NAXES": 2, "NAXIS1": naxis1, "NAXIS2": naxis2,
			"CRVAL1": centre_longitude, "CRVAL2": 0.0,
			"CTYPE1": f"{longitude}{projection}", "CTYPE2": f"{latitude}{projection}",
			"CRPIX1": naxis1 / 2 + 0.5, "CRPIX2": naxis2 / 2 + 0.5,
			"CDELT1": -360.0 / naxis1,
			"CDELT2": 180.0 / naxis2,
		}
		cards.update(extra)
		if rotation_deg != 0:
			cards["CROTA2"] = rotation_deg
		return cards, 0.5, (naxis1, naxis2)

	def _all_sky_centre(self, center, display_system: SkySystem):
		'''
		Resolve an all-sky ``center`` to ``(lon_deg, lat_deg)`` in the display
		system. A bare sequence is BY DEFINITION already in the display system;
		only a ``SkyCoord`` carries its own frame and gets converted.
		'''
		from astropy.coordinates import SkyCoord

		if isinstance(center, SkyCoord):
			target = Ast.SkyFrame(f"System={display_system.value}")
			converted = bridge.to_frame_units(center, target, squeeze=True)
			degrees = bridge.from_frame_units(converted, target)
			return float(degrees[0]), float(degrees[1])

		try:
			values = list(center)
		except TypeError as error:
			raise TypeError("'center' must be a SkyCoord or a 2-sequence "
			                f"(lon_deg, lat_deg); got {type(center).__name__}.") from error
		if len(values) != 2:
			raise ValueError("'center' must be a SkyCoord or a 2-sequence (lon_deg, lat_deg); "
			                 f"got {len(values)} value(s).")
		longitude, latitude = (float(values[0]), float(values[1]))
		if not (math.isfinite(longitude) and math.isfinite(latitude)):
			raise ValueError(f"'center' must be finite; got ({longitude}, {latitude}).")
		return longitude, latitude

	@staticmethod
	def _convert_centre(centre, extent_frame, display_system: SkySystem):
		'''
		Convert the bounding-circle centre into the display system.

		The centre comes back in the extent frame's own system. Writing e.g.
		``RA---TAN`` cards around a galactic-longitude CRVAL would be silent
		mislabelling — the bug class this constructor refuses everywhere.
		'''
		extent_system = getattr(extent_frame, "System", None)
		if extent_system == display_system.value:
			return centre

		target = Ast.SkyFrame(f"System={display_system.value}")
		conversion = extent_frame.convert(target)
		if conversion is None:
			raise ValueError(f"AST found no conversion from the extent's '{extent_system}' frame "
			                 f"to the requested display system '{display_system}'.")
		radians = bridge.to_frame_units(np.asarray(centre, dtype=float), extent_frame, squeeze=True)
		converted = conversion.tran(radians.reshape(2, 1), True)
		degrees = bridge.from_frame_units(converted.reshape(2), target)
		return float(degrees[0]), float(degrees[1])

	def _build(self, *, cards, aspect, naxis, width, char_height, margin,
	           palette, background, plot_id, grid) -> None:
		''' Shared final construction, used by ``__init__`` and ``fromFITSHeader``. '''
		naxis1, naxis2 = naxis

		self._frameset = ASTFrameSet.fromFITSHeader(fits_header=dict(cards))
		self._aspect = aspect
		self._width = width
		self._char_height = char_height
		self._margin = margin
		self._background = background
		self._base_palette = list(palette) if palette is not None else list(DEFAULT_PALETTE)

		self._grf = GrfSequence(char_height=char_height, palette=self._base_palette,
		                        metrics=FontMetrics())
		self._gbox = (0.0, 0.0, 1.0, aspect)
		bbox = (0.5, 0.5, naxis1 + 0.5, naxis2 + 0.5)

		# the frozen grey is an index into the DEFAULT palette; a shorter custom
		# palette resolves it the same way §3.3 resolves a stranded index
		grey = GRID_COLOUR if GRID_COLOUR < len(self._base_palette) else 1
		options = ",".join([
			f"Grid={1 if grid else 0}",
			"Border=1",                        # explicit, though it is AST's default:
			                                   # with Border=1 grid() draws the frame itself
			f"Colour(Grid)={grey}", "Width(Grid)=0.5",
			f"Colour(Border)={grey}", "Width(Border)=2",
			"TickAll=0", "MinTickLen(1)=0", "MinTickLen(2)=0",
		])
		# every argument positional: pyast constructors silently swallow kwargs,
		# and the options string is the 5th positional (D15 comment convention)
		self.astPlot = Ast.Plot(self._frameset.astObject, self._gbox, bbox, self._grf, options)

		self._layers = {}
		self._layer_sequence = 0

		if plot_id is None:
			SkyPlot._plot_counter += 1
			plot_id = f"plot-{os.getpid()}-{SkyPlot._plot_counter}"
		elif not ID_PATTERN.match(str(plot_id)):
			raise ValueError(f"'{plot_id}' is not a legal plot id (must match {ID_PATTERN.pattern}).")
		self._plot_id = str(plot_id)

		if grid:
			self.grid()

	# ------------------------------------------------------------------
	# layer plumbing
	# ------------------------------------------------------------------

	def _next_layer_id(self, kind: str, name: Optional[str]) -> str:
		'''
		Resolve a layer id (§2.3). Auto ids are ``<kind>-<seq>`` with ONE
		plot-wide counter shared across kinds — astplot-js's actual scheme.
		'''
		if kind in FIXED_LAYER_IDS:
			if name is not None:
				raise ValueError(f"The '{kind}' layer has a fixed id and takes no 'name'.")
			return FIXED_LAYER_IDS[kind]

		if name is not None:
			name = str(name)
			if not ID_PATTERN.match(name):
				raise ValueError(f"'{name}' is not a legal layer name (must match {ID_PATTERN.pattern}).")
			existing = self._layers.get(name)
			if existing is not None and existing["kind"] != kind:
				raise ValueError(f"Layer '{name}' already exists with kind '{existing['kind']}'; "
				                 f"refusing to replace it with a '{kind}' layer.")
			return name

		self._layer_sequence += 1
		return f"{kind}-{self._layer_sequence}"

	def _capture(self, layer_id: str, kind: str, draw) -> str:
		'''
		Run *draw* and record everything it caused AST to emit as one layer.

		Every draw method touches exactly one layer, which is what makes
		live-mode's one-envelope-per-call POST atomic.
		'''
		self._grf.take_commands()           # discard anything left over; a layer starts clean
		draw()
		commands = self._grf.take_commands()
		self._warn_if_text_escapes_the_viewport(layer_id, commands)

		existing = self._layers.get(layer_id)
		sequence = existing["seq"] if existing is not None else self._next_creation_sequence()
		self._layers[layer_id] = {
			"kind": kind,
			"z": DEFAULT_Z.get(kind, 0),
			"commands": commands,
			"seq": sequence,
		}
		return layer_id

	def _next_creation_sequence(self) -> int:
		return 1 + max((layer["seq"] for layer in self._layers.values()), default=0)

	def _warn_if_text_escapes_the_viewport(self, layer_id: str, commands) -> None:
		'''
		§4.4 fail-loud guard. The viewport crop stays (a *dynamic* margin would
		make the §4.3 viewport a function of content, breaking byte-stable
		rendering and the fixed aspect), but a label cropped out of existence is
		never silent: the layer and the text are named, and the user's remedies
		are ``margin=`` or shorter/re-anchored text.
		'''
		x1, y1, x2, y2 = self._gbox
		bleed = self._margin * max(x2 - x1, y2 - y1)
		for command in commands:
			if command.get("op") != "text":
				continue
			x, y = command["x"], command["y"]
			if not (x1 - bleed <= x <= x2 + bleed and y1 - bleed <= y <= y2 + bleed):
				logger.warning("Layer '%s': the text %r is anchored at (%.4f, %.4f), outside the "
				               "gbox+margin viewport; it will be cropped. Increase margin= or "
				               "shorten/re-anchor the text.", layer_id, command["text"], x, y)

	def _resolve_colour(self, colour, color, default) -> int:
		'''
		Resolve the ``colour``/``color`` synonym pair to the palette index AST
		needs (§1.1, §2.6).

		Both given and different is an error, never a silent preference — v1's
		``if color: colour = color`` silently favoured the synonym. Every colour
		parameter defaults to the ``None`` sentinel precisely so that "was this
		explicitly given" stays decidable.

		A CSS colour is turned into an index through the recorder's
		``ColToInt``, because ``astSet`` accepts **only** an integer for a
		``Colour(...)`` attribute — verified: ``Colour(Border)=red`` and
		``Colour(Border)=#4a7f7b`` are both rejected (AST falls through to the
		SkyFrame and reports the setting invalid), while ``Colour(Border)=11``
		is accepted. Growing the palette here is also exactly what makes the
		colour serialize as a literal CSS string in the layer (D18).
		'''
		if colour is not None and color is not None and colour != color:
			raise ValueError(f"'colour' and 'color' are synonyms but were given different values "
			                 f"({colour!r} and {color!r}); pass only one.")
		value = colour if colour is not None else color
		if value is None:
			return self._default_colour_index(default)
		normalized = normalize_colour(value, self._base_palette)
		if isinstance(normalized, int):
			return normalized
		return self._resolve_css_colour(normalized)

	def _default_colour_index(self, default) -> int:
		'''
		Resolve one of the §2.6 *defaults* to a palette index.

		The frozen defaults are stated as indices into the 16-entry
		``DEFAULT_PALETTE`` (grid grey 11, points blue 4, foreground 1). A caller
		who supplies a SHORTER custom palette has not asked for any particular
		colour, so refusing to draw would be hostile — but silently picking an
		arbitrary one would be worse.

		The resolution mirrors §3.3's frozen rule for exactly this condition (an
		out-of-range colour index renders as foreground, with a diagnostic),
		applied at production time rather than render time. Index 1 always
		exists: a palette needs at least 2 entries.

		A USER-supplied out-of-range index still raises — that is a stated
		colour the plot cannot honour, not an unstated default.
		'''
		if isinstance(default, str):
			return self._resolve_css_colour(default)
		index = int(default)
		if 0 <= index < len(self._base_palette):
			return index
		logger.warning("The default colour index %d is outside this plot's %d-entry palette; "
		               "using the foreground (index 1) instead. Pass colour= explicitly to "
		               "choose one.", index, len(self._base_palette))
		return 1

	def _resolve_css_colour(self, css: str) -> int:
		''' A CSS colour string as a palette index, growing the palette (D18). '''
		index = self._grf.ColToInt(css)
		if index < 0:
			raise ValueError(f"The recorder could not resolve the colour {css!r} "
			                 "to a palette index.")
		return index

	def _set_around(self, settings, draw):
		'''
		Apply AST plot attributes, draw, and restore them — the v1
		set-around-restore pattern, with the values routed through the quoting
		seam.
		'''
		saved = {}
		for name, value in settings.items():
			saved[name] = self.astPlot.get(name)
			_set_plot_attribute(self.astPlot, name, value)
		try:
			draw()
		finally:
			for name, value in saved.items():
				_set_plot_attribute(self.astPlot, name, value)

	@staticmethod
	def _line_style(style) -> int:
		if isinstance(style, LineStyle):
			return int(style)
		value = int(style)
		if value not in tuple(int(member) for member in LineStyle):
			raise ValueError(f"'{style}' is not a known line style; expected one of "
			                 f"{[int(member) for member in LineStyle]} (see cornish.enums.LineStyle).")
		return value

	# ------------------------------------------------------------------
	# drawing methods (§2.4)
	# ------------------------------------------------------------------

	def grid(self, visible: bool = True, *, colour=None, color=None,
	         style: LineStyle = LineStyle.SOLID, width: float = 0.5) -> str:
		'''
		(Re)draw the coordinate grid **and the frame border** into the one layer
		``'grid'``.

		With ``Border=1`` set in the construction options (also AST's default),
		a single ``astPlot.grid()`` call already draws the frame border itself —
		an explicit ``.border()`` call after it emits a byte-identical duplicate
		border polyline, so there is deliberately none. A single layer also
		keeps the invariant that every draw method touches exactly one layer.

		*colour*/*style*/*width* apply to the **grid lines only**; the border
		keeps its §2.6 construction defaults. Border restyling is a power-user
		job via ``.astPlot``.

		``visible=False`` is DECLARATIVE: it ensures the layer is absent and is a
		no-op when it already is. The imperative :meth:`removeLayer` is the call
		that raises on absence.

		:returns: the layer id, ``'grid'``
		'''
		if not visible:
			self._layers.pop("grid", None)
			return "grid"

		resolved = self._resolve_colour(colour, color, GRID_COLOUR)
		settings = {
			"Grid": 1,
			"Colour(Grid)": resolved,
			"Style(Grid)": self._line_style(style),
			"Width(Grid)": float(width),
		}
		return self._capture("grid", "grid",
		                     lambda: self._set_around(settings, self.astPlot.grid))

	def addRegionOutline(self, region: Union[ASTRegion, Ast.Region], *,
	                     colour: Union[str, int, None] = None, color=None,
	                     style: LineStyle = LineStyle.SOLID, width: float = 1.5,
	                     name: Optional[str] = None) -> str:
		'''
		Overlay the outline of *region*.

		:param region: an :class:`ASTRegion` (including :class:`ASTMoc`) or a raw ``Ast.Region``
		:returns: the layer id
		'''
		if isinstance(region, ASTRegion):
			ast_region = region.astObject
		elif isinstance(region, Ast.Region):
			ast_region = region
		else:
			raise TypeError("'region' must be an ASTRegion or a starlink.Ast.Region; "
			                f"got {type(region).__name__}.")

		resolved = self._resolve_colour(colour, color, REGION_COLOUR)
		settings = {
			"Colour(Border)": resolved,
			"Style(Border)": self._line_style(style),
			"Width(Border)": float(width),
		}
		layer_id = self._next_layer_id("region", name)
		return self._capture(layer_id, "region",
		                     lambda: self._set_around(settings,
		                                              lambda: self.astPlot.regionoutline(ast_region)))

	def addMoc(self, moc, *, colour: Union[str, int, None] = None, color=None,
	           style: LineStyle = LineStyle.SOLID, width: float = 1.0,
	           name: Optional[str] = None) -> str:
		'''
		Overlay a MOC's coverage outline. Sugar over :meth:`addRegionOutline`
		with its own layer kind, so consumers can style and toggle coverage maps
		apart from exact regions.
		'''
		ast_region = moc.astObject if isinstance(moc, ASTRegion) else moc
		if not isinstance(ast_region, Ast.Region):
			raise TypeError(f"'moc' must be an ASTMoc or an Ast.Region; got {type(moc).__name__}.")

		resolved = self._resolve_colour(colour, color, REGION_COLOUR)
		settings = {
			"Colour(Border)": resolved,
			"Style(Border)": self._line_style(style),
			"Width(Border)": float(width),
		}
		layer_id = self._next_layer_id("moc", name)
		return self._capture(layer_id, "moc",
		                     lambda: self._set_around(settings,
		                                              lambda: self.astPlot.regionoutline(ast_region)))

	def addPoints(self, points=None, *, ra=None, dec=None,
	              marker: Union[Marker, str] = Marker.CIRCLE, size: float = 1.0,
	              colour: Union[str, int, None] = None, color=None,
	              name: Optional[str] = None) -> str:
		'''
		Mark a set of sky positions.

		Every accepted input form is the SPEC-04 bridge contract: an
		``Ast.Plot`` *is* a FrameSet with the sky frame current, so SkyCoord /
		Quantity / bare-degrees / parallel arrays / single points are all
		handled uniformly.

		:param points: any bridge point form
		:param ra: longitudes, paired with *dec* (scalars make a single point; arrays are parallel)
		:param dec: latitudes, paired with *ra*
		:param marker: a :class:`~cornish.enums.Marker` or one of its string names
		'''
		if all(value is None for value in (points, ra, dec)):
			raise ValueError("'ra' and 'dec', OR 'points', must be provided.")
		if points is not None and (ra is not None or dec is not None):
			raise ValueError("Provide 'ra' and 'dec' OR 'points', not both.")
		if (ra is None) != (dec is None):
			raise ValueError("'ra' and 'dec' must be given together "
			                 f"(got ra={'None' if ra is None else 'a value'}, "
			                 f"dec={'None' if dec is None else 'a value'}).")

		marker_code = self._marker_code(marker)
		resolved = self._resolve_colour(colour, color, POINTS_COLOUR)

		if ra is not None:
			if isinstance(ra, u.Quantity) or isinstance(dec, u.Quantity):
				stacked = u.Quantity([u.Quantity(ra, u.deg), u.Quantity(dec, u.deg)])
			else:
				stacked = np.stack([np.asarray(ra), np.asarray(dec)])
			coordinates = bridge.to_frame_units(stacked, self.astPlot,
			                                    parallel_axes=(True if np.ndim(ra) > 0 else None))
		else:
			try:
				empty = (len(points) == 0)
			except TypeError:
				empty = False               # scalar SkyCoord/Quantity forms have no len()
			if empty:
				raise ValueError("No points were provided to plot.")
			coordinates = bridge.to_frame_units(points, self.astPlot)

		if coordinates.ndim == 1:
			coordinates = coordinates.reshape(2, 1)

		settings = {"Colour(Markers)": resolved, "Size(Markers)": float(size)}
		layer_id = self._next_layer_id("points", name)
		# ONE mark() call with the (2, n) array — not one call per point
		return self._capture(layer_id, "points",
		                     lambda: self._set_around(settings,
		                                              lambda: self.astPlot.mark(coordinates, marker_code)))

	@staticmethod
	def _marker_code(marker) -> int:
		if isinstance(marker, Marker):
			return int(marker)
		if isinstance(marker, str):
			try:
				return int(MARKER_INPUT_NAMES[marker.strip().lower()])
			except KeyError as error:
				raise ValueError(f"'{marker}' is not a known marker name; expected one of "
				                 f"{sorted(MARKER_INPUT_NAMES)}.") from error
		code = int(marker)
		if code not in tuple(int(member) for member in Marker):
			raise ValueError(f"{marker!r} is not a known marker code; expected one of "
			                 f"{[int(member) for member in Marker]}.")
		return code

	def addText(self, text: str, position, *, just: str = "CC",
	            up: tuple = (0.0, 1.0), size: float = 1.0,
	            font=None, colour: Union[str, int, None] = None, color=None,
	            name: Optional[str] = None) -> str:
		'''
		Draw a text string at a sky position.

		:param position: any single-point bridge form
		:param just: 2-character justification, vertical ``T``/``C``/``B``/``M`` then horizontal ``L``/``C``/``R``
		'''
		if not isinstance(text, str):
			raise TypeError(f"'text' must be a str; got {type(text).__name__}.")
		just = str(just).upper()
		if len(just) != 2 or just[0] not in "TCBM" or just[1] not in "LCR":
			raise ValueError(f"'{just}' is not a legal justification (vertical T/C/B/M then "
			                 "horizontal L/C/R).")
		up_vector = [float(value) for value in up]
		if len(up_vector) != 2 or (up_vector[0] == 0 and up_vector[1] == 0):
			raise ValueError(f"'up' must be a non-zero 2-vector; got {up!r}.")

		font_code = int(FontStyle.NORMAL if font is None else font)
		if font_code not in tuple(int(member) for member in FontStyle):
			raise ValueError(f"{font!r} is not a known font index; expected one of "
			                 f"{[int(member) for member in FontStyle]}.")

		coordinates = bridge.to_frame_units(position, self.astPlot, squeeze=True)
		resolved = self._resolve_colour(colour, color, FOREGROUND_COLOUR)
		settings = {"Colour(Strings)": resolved, "Size(Strings)": float(size),
		            "Font(Strings)": font_code}
		layer_id = self._next_layer_id("text", name)
		return self._capture(
			layer_id, "text",
			lambda: self._set_around(settings,
			                         lambda: self.astPlot.text(text, coordinates, up_vector, just)))

	def addImage(self, data, *, vmin: Optional[float] = None, vmax: Optional[float] = None,
	             cmap: str = "greys", name: Optional[str] = None) -> str:
		'''
		Underlay a greyscale image.

		**Stage-0 constraint (frozen)**: the image pixel grid IS the plot's own
		bbox grid — ``data.shape`` must equal ``(naxis2, naxis1)`` of the
		synthetic header. There is deliberately no ``wcs`` parameter yet: a
		frameset argument would need a value-equality check pyast cannot
		honestly provide, and "must equal the plot's own frameset" conveys zero
		information. Stage 3 adds ``wcs=`` together with the
		``Mapping.resample`` regridding that makes it meaningful.

		The input follows the FITS/bbox convention (row 0 = bottom).
		'''
		values = np.asarray(data)
		expected = (self._naxis2, self._naxis1)
		if values.shape != expected:
			raise ValueError(
				f"addImage in stage 0 requires the image grid to BE the plot's bbox grid: "
				f"expected shape {expected} (naxis2, naxis1), got {values.shape}. "
				"Regridding an image in its own WCS onto the plot is SPEC-09 stage-3 work "
				"(Mapping.resample).")

		href = array_to_png_data_uri(values, vmin=vmin, vmax=vmax, cmap=cmap)
		layer_id = self._next_layer_id("image", name)
		x1, y1, x2, y2 = self._gbox
		command = {"op": "image", "x1": x1, "y1": y1, "x2": x2, "y2": y2, "href": href}

		sequence = (self._layers[layer_id]["seq"] if layer_id in self._layers
		            else self._next_creation_sequence())
		self._layers[layer_id] = {"kind": "image", "z": DEFAULT_Z["image"],
		                          "commands": [command], "seq": sequence}
		return layer_id

	def addCompass(self, *, position=None, length: Optional[u.Quantity] = None,
	               colour: Union[str, int, None] = None, color=None) -> str:
		'''
		Draw a north/east compass (R25).

		The two arms are **geodesics whose initial bearings are exactly north
		and east**, computed with ``Frame.offset2``. Naive coordinate offsets
		(``alpha + delta_alpha``) would not point east anywhere but the equator.
		'''
		frame = bridge.current_frame_of(self.astPlot)
		anchor, length_radians = self._compass_anchor(position, length, frame)

		# offset2's position angle is measured from +axis-2 toward +axis-1:
		# 0 is due north, pi/2 due east [verified]
		_, north = frame.offset2(anchor, 0.0, length_radians)
		_, east = frame.offset2(anchor, math.pi / 2, length_radians)

		resolved = self._resolve_colour(colour, color, FOREGROUND_COLOUR)
		settings = {"Colour(Curves)": resolved, "Colour(Strings)": resolved}

		def draw():
			self.astPlot.curve(anchor, north)
			self.astPlot.curve(anchor, east)
			self.astPlot.text("N", north, [0.0, 1.0], "CC")
			self.astPlot.text("E", east, [0.0, 1.0], "CC")

		return self._capture("compass", "compass", lambda: self._set_around(settings, draw))

	def _compass_anchor(self, position, length, frame):
		'''
		Resolve the compass anchor and arm length.

		The default anchor is 12% in from the top-right gbox corner, inverse
		mapped to the sky; a field where the arms leave the gbox falls back to
		the top-left corner and then gives up — a compass on a field containing
		the pole is ill-defined, so it raises rather than drawing something
		wrong.
		'''
		x1, y1, x2, y2 = self._gbox
		if length is None:
			shorter = min(x2 - x1, y2 - y1)
			length_radians = self._gbox_span_to_radians(0.12 * shorter)
		else:
			if not isinstance(length, u.Quantity):
				raise TypeError(f"'length' must be an angular Quantity; got {type(length).__name__}.")
			length_radians = bridge.to_frame_distance(length, frame)
			if length_radians <= 0:
				raise ValueError(f"'length' must be positive; got {length}.")

		if position is not None:
			return bridge.to_frame_units(position, self.astPlot, squeeze=True), length_radians

		for fraction_x, corner in ((0.88, "top-right"), (0.12, "top-left")):
			anchor = self._sky_at_gbox(x1 + fraction_x * (x2 - x1), y1 + 0.88 * (y2 - y1))
			if anchor is None:
				continue
			if self._arms_stay_inside(anchor, length_radians, frame):
				return anchor, length_radians

		raise ValueError("A compass could not be placed: the north/east arms leave the plot at "
		                 "both candidate corners. A compass on a field containing the pole is "
		                 "ill-defined — pass position= explicitly if you know where you want it.")

	#: How far past the gbox edge a compass tip may land and still count as
	#: inside, as a fraction of the gbox extent.
	#:
	#: The default placement is a knife-edge BY CONSTRUCTION: the anchor sits 12%
	#: in from the corner and the arms are 12% of the field long, so a tip lands
	#: exactly ON the boundary. Projection curvature then pushes it a hair past —
	#: measured at y = 1.00012 for the golden TAN field, i.e. 0.012%. A strictly
	#: exclusive test would therefore reject the compass on ordinary fields,
	#: which inverts the spec's intent (ordinary fields get a compass; fields
	#: containing the pole raise). This tolerance admits the boundary while
	#: staying far too small to rescue a genuinely escaping arm, which lands
	#: wildly outside or fails to map at all.
	_COMPASS_EDGE_TOLERANCE = 0.005

	def _arms_stay_inside(self, anchor, length_radians, frame) -> bool:
		x1, y1, x2, y2 = self._gbox
		slack_x = self._COMPASS_EDGE_TOLERANCE * (x2 - x1)
		slack_y = self._COMPASS_EDGE_TOLERANCE * (y2 - y1)
		for angle in (0.0, math.pi / 2):
			_, tip = frame.offset2(anchor, angle, length_radians)
			graphics = self._gbox_at_sky(tip)
			if graphics is None:
				return False
			if not (x1 - slack_x <= graphics[0] <= x2 + slack_x
			        and y1 - slack_y <= graphics[1] <= y2 + slack_y):
				return False
		return True

	def _gbox_span_to_radians(self, span: float) -> float:
		''' The angular length a gbox-unit span subtends at the field centre. '''
		x1, y1, x2, y2 = self._gbox
		centre_x, centre_y = (x1 + x2) / 2, (y1 + y2) / 2
		first = self._sky_at_gbox(centre_x, centre_y)
		second = self._sky_at_gbox(centre_x + span, centre_y)
		if first is None or second is None:
			raise ValueError("The plot centre does not map to the sky, so a default compass "
			                 "or scale-bar length cannot be derived; pass length= explicitly.")
		frame = bridge.current_frame_of(self.astPlot)
		return float(frame.distance(first, second))

	def _sky_at_gbox(self, x: float, y: float):
		''' Inverse-map a gbox point to sky coordinates (radians), or ``None``. '''
		result = self.astPlot.tran(np.array([[float(x)], [float(y)]]), True).reshape(2)
		if not all(math.isfinite(value) and value != Ast.BAD for value in result):
			return None
		return result

	def _gbox_at_sky(self, position):
		''' Forward-map sky coordinates (radians) to a gbox point, or ``None``. '''
		result = self.astPlot.tran(np.asarray(position, dtype=float).reshape(2, 1), False).reshape(2)
		if not all(math.isfinite(value) and value != Ast.BAD for value in result):
			return None
		return result

	def addScaleBar(self, length: u.Quantity, *, label: Optional[str] = None,
	                colour: Union[str, int, None] = None, color=None) -> str:
		'''
		Draw a scale bar of exactly *length* on the sky (R25).

		``astPlot.curve`` draws great circles, so the endpoints are chosen to
		make the **geodesic** length exact rather than the parallel arc: at
		declination δ the equal-declination great-circle inversion gives
		``Δα = 2·asin( sin(length/2) / cos δ )``.
		'''
		if not isinstance(length, u.Quantity):
			raise TypeError(f"'length' must be an angular Quantity; got {type(length).__name__}.")
		frame = bridge.current_frame_of(self.astPlot)
		length_radians = bridge.to_frame_distance(length, frame)
		if not math.isfinite(length_radians) or length_radians <= 0:
			raise ValueError(f"'length' must be a positive angle; got {length}.")

		x1, y1, x2, y2 = self._gbox
		start, finish, midpoint = self.scaleBarEndpoints(length_radians)

		field_width = self._gbox_span_to_radians(x2 - x1)
		if length_radians > field_width:
			raise ValueError(f"A {length} bar is longer than the plotted field "
			                 f"({field_width * 180 / math.pi:.4g} deg across).")

		if label is None:
			label = format(length)

		resolved = self._resolve_colour(colour, color, FOREGROUND_COLOUR)
		settings = {"Colour(Curves)": resolved, "Colour(Strings)": resolved}

		def draw():
			self.astPlot.curve(start, finish)
			self.astPlot.text(str(label), midpoint, [0.0, 1.0], "BC")

		return self._capture("scalebar", "scalebar", lambda: self._set_around(settings, draw))

	def scaleBarEndpoints(self, length_radians: float):
		'''
		The two sky endpoints of a scale bar of *length_radians*, and its
		midpoint, all in the current frame's radians.

		``astPlot.curve`` draws GREAT CIRCLES, so a bar drawn between two points
		of equal declination is a geodesic, not the parallel arc — and the two
		differ. Choosing the endpoints by the equal-declination great-circle
		inversion is what makes the drawn bar's length exactly the requested
		angle::

		    Δα = 2·asin( sin(length/2) / cos δ )

		Exposed (rather than buried in :meth:`addScaleBar`) because "the bar is
		exactly 30 arcmin long" is the property worth testing directly.

		:raises ValueError: if the anchor does not map to the sky, or the bar cannot fit at that latitude
		'''
		x1, y1, x2, y2 = self._gbox
		anchor = self._sky_at_gbox(x1 + 0.08 * (x2 - x1), y1 + 0.08 * (y2 - y1))
		if anchor is None:
			raise ValueError("The scale-bar anchor does not map to the sky "
			                 "(the bottom-left corner of this projection is blank).")

		longitude, latitude = float(anchor[0]), float(anchor[1])
		cos_latitude = math.cos(latitude)
		half = math.sin(float(length_radians) / 2.0)
		if cos_latitude <= 0 or half > cos_latitude:
			raise ValueError("That bar cannot fit at this latitude: it would need to span more "
			                 "than the full small circle of constant declination there.")
		half_span = math.asin(half / cos_latitude)

		return (np.array([longitude - half_span, latitude]),
		        np.array([longitude + half_span, latitude]),
		        np.array([longitude, latitude]))

	# ------------------------------------------------------------------
	# construction from an existing WCS
	# ------------------------------------------------------------------

	@classmethod
	def fromFITSHeader(cls, header, *, width: int = 640, char_height: float = 0.02,
	                   margin: float = 0.15, palette: Optional[Sequence[str]] = None,
	                   background: Optional[str] = None, grid: bool = True,
	                   plot_id: Optional[str] = None) -> "SkyPlot":
		'''
		Construct a SkyPlot in an **existing complete WCS**, bypassing the
		synthetic-header builder (§2.4).

		``extent``/``projection``/``system``/``rotation``/``aspect`` are
		meaningless here and deliberately absent. Parameter validation is *not*
		bypassed: the same shared helper runs.

		Headers with non-square pixels are REFUSED rather than rendered subtly
		distorted — gbox isotropy, and therefore the honesty of ``Scales()`` and
		``Qch()``, is impossible for them under this rule.

		:raises cornish.exc.IncompleteHeader: if NAXIS1/NAXIS2 are missing
		:raises ValueError: if the pixels are not square, or the NAXIS ratio is outside the aspect domain
		'''
		from ..channel import ASTFITSChannel

		char_height, margin, width, palette, background = _validate_common_parameters(
			char_height=char_height, margin=margin, width=width,
			palette=palette, background=background)

		channel = ASTFITSChannel(header=header)
		dimensions = channel.dimensions        # raises IncompleteHeader when NAXIS/NAXISn missing
		if len(dimensions) < 2:
			raise IncompleteHeader(f"The header declares {len(dimensions)} axis/axes; a sky plot "
			                       "needs NAXIS1 and NAXIS2.")
		naxis1, naxis2 = int(dimensions[0]), int(dimensions[1])
		if naxis1 <= 0 or naxis2 <= 0:
			raise ValueError(f"NAXIS1/NAXIS2 must be positive; got {naxis1} and {naxis2}.")

		aspect = naxis2 / naxis1
		if not MINIMUM_ASPECT <= aspect <= MAXIMUM_ASPECT:
			raise ValueError(f"The header's NAXIS2/NAXIS1 ratio is {aspect}, outside the "
			                 f"supported aspect domain [{MINIMUM_ASPECT}, {MAXIMUM_ASPECT}].")

		plot = cls.__new__(cls)
		plot._cards = dict(header) if isinstance(header, dict) else header
		plot._build(cards=header, aspect=aspect, naxis=(naxis1, naxis2), width=width,
		            char_height=char_height, margin=margin, palette=palette,
		            background=background, plot_id=plot_id, grid=False)
		plot._refuse_non_square_pixels()
		if grid:
			plot.grid()
		return plot

	def _refuse_non_square_pixels(self) -> None:
		'''
		Detect non-square pixels in a way that is **card-form independent**:
		CD/PC-matrix headers may carry no CDELT cards at all, so reading cards
		would silently miss exactly the case that matters. Instead transform the
		bbox-centre pixel and its two unit neighbours through the frameset and
		measure both steps with the current frame's own ``distance()`` — pure
		AST, uniform across CDELT/CD/PC forms.
		'''
		frameset = self._frameset.astObject
		naxis1, naxis2 = self._naxis1, self._naxis2
		centre = np.array([[naxis1 / 2.0], [naxis2 / 2.0]])
		steps = np.hstack([centre,
		                   centre + np.array([[1.0], [0.0]]),
		                   centre + np.array([[0.0], [1.0]])])
		transformed = frameset.tran(steps, True)
		frame = bridge.current_frame_of(frameset)
		first = float(frame.distance(transformed[:, 0], transformed[:, 1]))
		second = float(frame.distance(transformed[:, 0], transformed[:, 2]))
		if not (math.isfinite(first) and math.isfinite(second)) or max(first, second) == 0:
			raise ValueError("The header's pixel scales could not be measured through its own "
			                 "WCS, so square pixels cannot be confirmed.")
		if abs(first - second) / max(first, second) > 0.01:
			degrees = 180.0 / math.pi
			raise ValueError(
				f"This header has non-square pixels ({first * degrees:.6g} deg/pixel along axis 1 "
				f"vs {second * degrees:.6g} deg/pixel along axis 2). Stage 0 refuses them rather "
				"than rendering subtly distorted text layout: gbox isotropy — and therefore the "
				"honesty of Scales() and Qch() — is impossible for them.")

	# ------------------------------------------------------------------
	# output
	# ------------------------------------------------------------------

	@property
	def frameset(self) -> ASTFrameSet:
		''' The pixel->sky frameset backing the plot. '''
		return self._frameset

	@property
	def plot_id(self) -> str:
		''' The id used by the viewer. '''
		return self._plot_id

	@property
	def layers(self) -> list:
		''' Live layer ids, in render order. '''
		return [layer["id"] for layer in self._ordered_layers()]

	@property
	def _naxis1(self) -> int:
		return int(self._cards["NAXIS1"])

	@property
	def _naxis2(self) -> int:
		return int(self._cards["NAXIS2"])

	def removeLayer(self, name: str) -> None:
		'''
		Remove a layer. Unknown ids raise — this is the imperative call, unlike
		``grid(visible=False)``, which is declarative.

		:raises KeyError: if no such layer exists
		'''
		if name not in self._layers:
			raise KeyError(f"No layer '{name}' on this plot; live layers are {sorted(self._layers)}.")
		del self._layers[name]

	def _ordered_layers(self) -> list:
		''' Layers sorted by (z, creation sequence) — §3.5. '''
		entries = [dict(layer, id=layer_id) for layer_id, layer in self._layers.items()]
		return sorted(entries, key=lambda layer: (layer["z"], layer["seq"]))

	def canvas(self) -> dict:
		'''
		The §3.2 canvas block.

		The palette written out is the BASE palette: colours grown by
		``ColToInt`` are inlined into their layers as literal CSS strings (D18),
		so the canvas stays a themeable table rather than a growing list.
		'''
		block = {
			"gbox": list(self._gbox),
			"aspect": self._aspect,
			"margin": self._margin,
			"char_height": self._char_height,
			"font_family": self._grf.metrics.font_family,
			"palette": list(self._base_palette),
			"width": self._width,
		}
		block["background"] = self._background
		return block

	def document(self) -> dict:
		''' The full draw-stream document — canvas plus every live layer in render order. '''
		document = {
			"cornish_view": SCHEMA_VERSION,
			"type": "draw_stream",
			"plot_id": self._plot_id,
			"canvas": self.canvas(),
			"layers": [{"id": layer["id"], "kind": layer["kind"], "z": layer["z"],
			            "commands": layer["commands"]}
			           for layer in self._ordered_layers()],
		}
		validate_document(document)
		return document

	def envelopes(self) -> list:
		''' The same content as a per-layer envelope list, canvas on the first. '''
		ordered = self._ordered_layers()
		if not ordered:
			return [make_envelope(canvas=self.canvas(), plot_id=self._plot_id)]
		results = []
		for index, layer in enumerate(ordered):
			results.append(make_envelope(
				layer_id=layer["id"], action=LayerAction.REPLACE, kind=layer["kind"],
				z=layer["z"], commands=layer["commands"], plot_id=self._plot_id,
				canvas=self.canvas() if index == 0 else None))
		return results

	def envelope_for(self, layer_id: str) -> dict:
		''' The single-layer envelope for *layer_id* (what live mode POSTs). '''
		if layer_id not in self._layers:
			raise KeyError(f"No layer '{layer_id}' on this plot.")
		layer = self._layers[layer_id]
		return make_envelope(layer_id=layer_id, action=LayerAction.REPLACE, kind=layer["kind"],
		                     z=layer["z"], commands=layer["commands"], plot_id=self._plot_id)

	def svg(self, *, width: Optional[int] = None) -> str:
		''' Render the plot to an SVG string. '''
		return render_svg(self.document(), width=self._width if width is None else width)

	def save(self, path) -> None:
		'''
		Write the plot to *path*.

		``.svg`` writes the rendered SVG; ``.json`` writes the canonical
		draw-stream document. Any other suffix raises — PNG and friends go
		through :meth:`to_matplotlib`.
		'''
		suffix = os.path.splitext(os.fspath(path))[1].lower()
		if suffix == ".svg":
			content = self.svg()
		elif suffix == ".json":
			content = canonical_json(self.document())
		else:
			raise ValueError(f"save() supports '.svg' and '.json'; got '{suffix or path}'. "
			                 "Raster formats go through to_matplotlib().")
		with open(path, "w", encoding="utf-8") as stream:
			stream.write(content)

	def _repr_mimebundle_(self, include=None, exclude=None) -> dict:
		''' Rich display in Jupyter: the plot renders itself as SVG. '''
		return {"image/svg+xml": self.svg()}

	def to_matplotlib(self, axes=None):
		'''
		Replay the captured document into a matplotlib Axes (imports matplotlib
		here, and only here).

		This is a **replay of frozen layout** (D17): geometry replays verbatim in
		gbox coordinates on an equal-aspect axes, and matplotlib does no label
		placement of its own.

		:returns: the :class:`matplotlib.figure.Figure`
		'''
		from matplotlib import pyplot

		from .envelope import ATTR_DEFAULTS

		document = self.document()
		x1, y1, x2, y2 = (float(value) for value in document["canvas"]["gbox"])
		margin = float(document["canvas"].get("margin", 0.15)) * max(x2 - x1, y2 - y1)
		palette = list(document["canvas"].get("palette", DEFAULT_PALETTE))
		char_height = float(document["canvas"].get("char_height", 0.02))

		if axes is None:
			figure = pyplot.figure()
			axes = figure.add_axes([0, 0, 1, 1])
		else:
			figure = axes.get_figure()
		axes.set_xlim(x1 - margin, x2 + margin)
		axes.set_ylim(y1 - margin, y2 + margin)
		axes.set_aspect("equal")
		axes.set_axis_off()

		def colour_of(value):
			if isinstance(value, str):
				return value
			index = int(value)
			return palette[index] if 0 <= index < len(palette) else palette[1]

		height = (y2 - y1) + 2 * margin
		for layer in document["layers"]:
			state = {"line": dict(ATTR_DEFAULTS), "mark": dict(ATTR_DEFAULTS),
			         "text": dict(ATTR_DEFAULTS)}
			for command in layer["commands"]:
				operation = command["op"]
				if operation == "attr":
					state[command["prim"]].update(
						{key: command[key] for key in
						 ("colour", "width", "style", "size", "font") if key in command})
				elif operation == "line":
					axes.plot(command["x"], command["y"],
					          color=colour_of(state["line"]["colour"]),
					          linewidth=float(state["line"]["width"]))
				elif operation == "mark":
					axes.plot(command["x"], command["y"], linestyle="none", marker="o",
					          color=colour_of(state["mark"]["colour"]),
					          markersize=3 * float(state["mark"]["size"]))
				elif operation == "text":
					# em size in points, from the axes height in points, so text
					# scales with the figure. Taken from the axes geometry rather
					# than a live renderer: this is a replay of frozen layout, and
					# consulting a renderer is exactly what makes grf_matplotlib's
					# metrics backend-dependent (D17).
					axes_height_points = axes.get_position().height * figure.get_figheight() * 72
					points = char_height * float(state["text"]["size"]) / height * axes_height_points
					axes.text(command["x"], command["y"], command["text"],
					          color=colour_of(state["text"]["colour"]),
					          fontsize=max(points, 1.0),
					          ha={"L": "left", "C": "center", "R": "right"}[command["just"][1]],
					          va={"T": "top", "C": "center", "B": "baseline",
					              "M": "bottom"}[command["just"][0]])
		return figure
