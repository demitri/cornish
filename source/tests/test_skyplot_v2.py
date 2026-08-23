
'''
``SkyPlot`` v2 behaviour (SPEC-09A §7.7).

The constructor is where most of this spec's fail-loud decisions live, so most
of what follows is a validation matrix: every ambiguous or under-specified
input must raise rather than be quietly resolved into a plot that lies about
what it shows.
'''

import logging
import math

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from cornish import ASTCircle, ASTICRSFrame, _pyast_bridge as bridge
from cornish.exc import IncompleteHeader
from cornish.enums import LineStyle, Marker
from cornish.mapping.frame.sky_frame import ASTSkyFrame
from cornish.plot.envelope import validate_document, validate_envelope
from cornish.plot.palette import DEFAULT_PALETTE
from cornish.plot.sky import SkyPlot

#: A complete, square-pixel WCS header for the fromFITSHeader tests.
SQUARE_HEADER = {
	"NAXIS": 2, "NAXIS1": 64, "NAXIS2": 32,
	"CRVAL1": 30.0, "CRVAL2": 45.0,
	"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
	"CRPIX1": 32.5, "CRPIX2": 16.5,
	"CDELT1": -0.01, "CDELT2": 0.01, "RADESYS": "ICRS",
}


def circle(centre=(30.0, 45.0), radius=2.0, frame=None):
	return ASTCircle(frame=frame or ASTICRSFrame(), center=list(centre), radius=radius * u.deg)


@pytest.fixture
def plot():
	return SkyPlot(extent=circle())


# ---------------------------------------------------------------------------
# constructor validation matrix (§2.1)
# ---------------------------------------------------------------------------

def test_extent_or_all_sky_is_required():
	with pytest.raises(ValueError, match="Exactly one"):
		SkyPlot()
	with pytest.raises(ValueError, match="Exactly one"):
		SkyPlot(extent=circle(), all_sky=True)


def test_center_is_all_sky_only():
	''' With an extent, a centre would contradict the extent's own centre. '''
	with pytest.raises(ValueError, match="only legal with all_sky"):
		SkyPlot(extent=circle(), center=[10.0, 0.0])


def test_extent_must_be_a_region():
	with pytest.raises(TypeError):
		SkyPlot(extent="a circle, honest")


def test_non_sky_extent_is_refused():
	''' A basic-frame extent has no sky coordinates to build WCS cards from. '''
	import starlink.Ast as Ast
	from cornish.mapping.frame.frame import ASTFrame
	basic = ASTCircle(frame=ASTFrame(naxes=2), center=[1.0, 1.0], radius=0.5)
	with pytest.raises(ValueError, match="sky frame"):
		SkyPlot(extent=basic)


def test_non_finite_bounding_circle_raises():
	'''
	v1 logged a warning here and continued to a guaranteed downstream failure —
	silent-skip-adjacent. v2 crashes at the source.
	'''
	class NaNCircle(ASTCircle):
		@property
		def center(self):
			return (float("nan"), float("nan"))

	bad = NaNCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg)
	with pytest.raises(ValueError, match="not finite"):
		SkyPlot(extent=bad)


@pytest.mark.parametrize("aspect", [0.04, 20.1, -1.0, float("nan")])
def test_aspect_outside_the_domain_raises(aspect):
	with pytest.raises(ValueError):
		SkyPlot(extent=circle(), aspect=aspect)


def test_aspect_is_quantized_and_reported():
	'''
	§2.2: the requested aspect is quantized to 1/100 so that gbox aspect equals
	naxis2·CDELT2 / (naxis1·CDELT1) EXACTLY — which is what makes Scales()=(1,1)
	and Qch's isotropy claim honest. The quantized value is reported, not hidden.
	'''
	plot = SkyPlot(extent=circle(), aspect=0.333)
	assert plot.document()["canvas"]["aspect"] == 0.33
	assert plot._cards["NAXIS2"] == 33


def test_all_sky_defaults():
	plot = SkyPlot(all_sky=True)
	cards = plot._cards
	assert (cards["NAXIS1"], cards["NAXIS2"]) == (200, 100)
	assert cards["CTYPE1"] == "RA---AIT"
	assert cards["CRVAL1"] == 0.0 and cards["CRVAL2"] == 0.0
	assert cards["CDELT1"] == -1.8 and cards["CDELT2"] == 1.8
	assert plot.document()["canvas"]["aspect"] == 0.5


def test_all_sky_projection_is_restricted():
	'''
	A zenithal code cannot represent the whole sky in this parameterization.
	Note the sentinel pattern: an OMITTED projection resolves to AIT, while an
	explicit TAN raises — a literal "TAN" default could not support both.
	'''
	assert SkyPlot(all_sky=True)._cards["CTYPE1"].endswith("AIT")
	with pytest.raises(ValueError, match="whole sky"):
		SkyPlot(all_sky=True, projection="TAN")
	assert SkyPlot(all_sky=True, projection="MOL")._cards["CTYPE1"].endswith("MOL")


@pytest.mark.parametrize("aspect", [1.0, 0.25, 2.0])
def test_all_sky_aspect_must_be_half(aspect):
	''' The full-sky field is inherently 2:1; another aspect could only lie. '''
	with pytest.raises(ValueError, match="2:1"):
		SkyPlot(all_sky=True, aspect=aspect)


def test_all_sky_accepts_exactly_half():
	assert SkyPlot(all_sky=True, aspect=0.5).document()["canvas"]["aspect"] == 0.5


def test_all_sky_centre_must_be_on_the_equator():
	''' The AIT/MOL fiducial point must lie on the equator — never silently zeroed. '''
	assert SkyPlot(all_sky=True, center=[120.0, 0.0])._cards["CRVAL1"] == 120.0
	with pytest.raises(ValueError, match="equator"):
		SkyPlot(all_sky=True, center=[120.0, 30.0])


def test_all_sky_centre_latitude_is_checked_after_conversion():
	'''
	§2.1, order of operations: the constraint is on CRVAL2, which lives in the
	DISPLAY system. A galactic b=0 SkyCoord displayed in ICRS therefore
	correctly FAILS, because its ICRS declination is not zero.
	'''
	galactic_equator = SkyCoord(l=120 * u.deg, b=0 * u.deg, frame="galactic")
	with pytest.raises(ValueError, match="equator"):
		SkyPlot(all_sky=True, center=galactic_equator, system="ICRS")
	# ...and succeeds when the display system is the one it is on the equator of
	SkyPlot(all_sky=True, center=galactic_equator, system="GALACTIC")


@pytest.mark.parametrize("center", [[1.0], [1.0, 2.0, 3.0], "middle"])
def test_all_sky_centre_shape_is_pinned(center):
	with pytest.raises((ValueError, TypeError)):
		SkyPlot(all_sky=True, center=center)


def test_unsupported_system_raises():
	with pytest.raises(ValueError, match="AZEL|no known sky system|FITS CTYPE"):
		SkyPlot(extent=circle(), system="AZEL")
	with pytest.raises(ValueError):
		SkyPlot(extent=circle(), system="nonsense")


def test_bad_projection_is_refused_by_ast():
	'''
	§2.1: no hand-maintained projection whitelist for the finite case — the
	synthetic header is built and AST is the source of truth.
	'''
	with pytest.raises(ValueError):
		SkyPlot(extent=circle(), projection="ZZZ")


@pytest.mark.parametrize("kwargs", [
	{"char_height": 0}, {"char_height": -1}, {"margin": 0.9}, {"margin": -0.1},
	{"width": 0}, {"palette": ["#000000"]}, {"palette": ["#000000", "nope"]},
	{"background": "not-a-colour"},
])
def test_common_parameter_domains(kwargs):
	with pytest.raises((ValueError, TypeError)):
		SkyPlot(extent=circle(), **kwargs)


def test_palette_uses_the_wider_entry_grammar():
	'''
	§2.1: palettes are validated with normalize_palette_entry, NOT
	normalize_colour — slot 0 is legitimately transparent, which the
	user-colour grammar rejects.
	'''
	SkyPlot(extent=circle(), palette=["transparent", "#000000", "#ff0000"])
	SkyPlot(extent=circle(), palette=list(DEFAULT_PALETTE))


def test_bad_plot_id_raises():
	with pytest.raises(ValueError, match="legal plot id"):
		SkyPlot(extent=circle(), plot_id="not a plot id")


# ---------------------------------------------------------------------------
# synthetic-header math (§2.2)
# ---------------------------------------------------------------------------

def test_cdelt1_is_negative():
	'''
	RA/longitude increases to the LEFT — the astronomical convention. v1 had
	CDELT1 positive, i.e. every v1 plot was mirrored east-west; this is a v1 bug
	being fixed, not a compatible choice.
	'''
	cards = SkyPlot(extent=circle())._cards
	assert cards["CDELT1"] < 0 < cards["CDELT2"]
	assert cards["CDELT1"] == -cards["CDELT2"]


def test_shorter_axis_spans_2_1_radii():
	cards = SkyPlot(extent=circle(radius=2.0))._cards
	assert cards["CDELT2"] * min(cards["NAXIS1"], cards["NAXIS2"]) == pytest.approx(2.1 * 2.0)


def test_crota2_present_exactly_when_rotation_is_nonzero():
	assert "CROTA2" not in SkyPlot(extent=circle())._cards
	assert "CROTA2" not in SkyPlot(extent=circle(), rotation=0)._cards
	assert SkyPlot(extent=circle(), rotation=30 * u.deg)._cards["CROTA2"] == 30.0
	assert SkyPlot(extent=circle(), rotation=30)._cards["CROTA2"] == 30.0   # bare number = degrees


def test_rotation_must_be_an_angle():
	with pytest.raises(ValueError):
		SkyPlot(extent=circle(), rotation=30 * u.m)


@pytest.mark.parametrize("system, ctype1, ctype2, extra", [
	("ICRS", "RA---TAN", "DEC--TAN", {"RADESYS": "ICRS"}),
	("FK5", "RA---TAN", "DEC--TAN", {"RADESYS": "FK5", "EQUINOX": 2000.0}),
	("FK4", "RA---TAN", "DEC--TAN", {"RADESYS": "FK4", "EQUINOX": 1950.0}),
	("GALACTIC", "GLON-TAN", "GLAT-TAN", {}),
	("ECLIPTIC", "ELON-TAN", "ELAT-TAN", {}),
])
def test_per_system_cards(system, ctype1, ctype2, extra):
	cards = SkyPlot(extent=circle(), system=system)._cards
	assert cards["CTYPE1"] == ctype1 and cards["CTYPE2"] == ctype2
	for key, value in extra.items():
		assert cards[key] == value


def test_display_system_is_inherited_from_the_extent():
	galactic = circle(centre=(120.0, 30.0), radius=1.0, frame=ASTSkyFrame(system="GALACTIC"))
	cards = SkyPlot(extent=galactic)._cards
	assert cards["CTYPE1"] == "GLON-TAN"
	assert cards["CRVAL1"] == pytest.approx(120.0)
	assert cards["CRVAL2"] == pytest.approx(30.0)


def test_extent_centre_is_converted_not_relabelled():
	'''
	§2.1: the centre comes back in the extent frame's system. Writing RA---TAN
	cards around a galactic-longitude CRVAL would be silent mislabelling — the
	bug class this constructor refuses everywhere. Checked against astropy.
	'''
	galactic = circle(centre=(120.0, 30.0), radius=1.0, frame=ASTSkyFrame(system="GALACTIC"))
	cards = SkyPlot(extent=galactic, system="ICRS")._cards
	expected = SkyCoord(l=120 * u.deg, b=30 * u.deg, frame="galactic").icrs
	assert cards["CTYPE1"] == "RA---TAN"
	# the centre was CONVERTED, not relabelled: 120 would mean relabelling
	assert cards["CRVAL1"] != pytest.approx(120.0, abs=1e-3)
	# AST and astropy differ on the galactic pole by ~1e-5 deg (0.04 arcsec),
	# so agreement is asserted at a tolerance that still excludes relabelling
	assert cards["CRVAL1"] == pytest.approx(expected.ra.deg, abs=1e-4)
	assert cards["CRVAL2"] == pytest.approx(expected.dec.deg, abs=1e-4)


def test_gbox_is_isotropic_and_scales_is_honest():
	plot = SkyPlot(extent=circle(), aspect=0.5)
	assert plot._gbox == (0.0, 0.0, 1.0, 0.5)
	assert plot._grf.Scales() == (1.0, 1.0)


# ---------------------------------------------------------------------------
# layers and naming (§2.3)
# ---------------------------------------------------------------------------

def test_layer_ids_use_one_plot_wide_counter(plot):
	''' §2.3: one counter shared across kinds — astplot-js's actual scheme. '''
	assert plot.addRegionOutline(circle()) == "region-1"
	assert plot.addPoints(ra=[30.0], dec=[45.0]) == "points-2"
	assert plot.addText("M51", [30.0, 45.0]) == "text-3"


def test_fixed_layer_ids(plot):
	assert plot.grid() == "grid"
	assert plot.addCompass() == "compass"
	assert plot.addScaleBar(30 * u.arcmin) == "scalebar"


def test_fixed_layers_replace_rather_than_accumulate(plot):
	plot.addCompass()
	plot.addCompass()
	assert plot.layers.count("compass") == 1


def test_named_layer_is_used_verbatim(plot):
	assert plot.addRegionOutline(circle(), name="fov") == "fov"
	assert "fov" in plot.layers


def test_named_layer_of_the_same_kind_replaces(plot):
	plot.addRegionOutline(circle(), name="fov")
	plot.addRegionOutline(circle(radius=1.0), name="fov")
	assert plot.layers.count("fov") == 1


def test_named_layer_collision_across_kinds_raises(plot):
	plot.addRegionOutline(circle(), name="thing")
	with pytest.raises(ValueError, match="already exists with kind"):
		plot.addPoints(ra=[30.0], dec=[45.0], name="thing")


def test_illegal_layer_name_raises(plot):
	with pytest.raises(ValueError, match="legal layer name"):
		plot.addRegionOutline(circle(), name="not a name")


def test_fixed_id_methods_expose_no_name_parameter():
	'''
	§2.3: ``grid``/``compass``/``scalebar`` have fixed ids, so those methods
	take no ``name`` at all — the invariant is enforced by the signature rather
	than by a runtime check.
	'''
	import inspect
	for method in (SkyPlot.grid, SkyPlot.addCompass, SkyPlot.addScaleBar):
		assert "name" not in inspect.signature(method).parameters, method.__name__
	for method in (SkyPlot.addRegionOutline, SkyPlot.addMoc, SkyPlot.addPoints,
	               SkyPlot.addText, SkyPlot.addImage):
		assert "name" in inspect.signature(method).parameters, method.__name__


def test_default_z_by_kind(plot):
	plot.addRegionOutline(circle())
	plot.addImage(np.zeros((100, 100)))
	by_id = {layer["id"]: layer["z"] for layer in plot.document()["layers"]}
	assert by_id["image-2"] == -20
	assert by_id["grid"] == -10
	assert by_id["region-1"] == 0


def test_render_order_puts_image_behind_grid_behind_the_rest(plot):
	plot.addRegionOutline(circle())
	plot.addImage(np.zeros((100, 100)))
	assert plot.layers == ["image-2", "grid", "region-1"]


def test_remove_layer_is_imperative(plot):
	plot.addRegionOutline(circle(), name="fov")
	plot.removeLayer("fov")
	assert "fov" not in plot.layers
	with pytest.raises(KeyError):
		plot.removeLayer("fov")


def test_grid_visible_false_is_declarative(plot):
	'''
	§2.4: ``grid(visible=False)`` ensures the layer is ABSENT and is a no-op
	when it already is; ``removeLayer`` is the call that raises on absence.
	'''
	plot.grid(visible=False)
	assert "grid" not in plot.layers
	plot.grid(visible=False)                    # no-op, must not raise
	assert "grid" not in plot.layers


def test_grid_false_at_construction_starts_empty():
	assert SkyPlot(extent=circle(), grid=False).layers == []


def test_grid_draws_the_border_itself():
	'''
	§2.4 [verified round 11]: with Border=1, grid() already draws the frame
	border, and an explicit border() call after it emits a byte-identical
	duplicate. There is therefore no separate border layer and no border() call.
	'''
	plot = SkyPlot(extent=circle())
	commands = plot.document()["layers"][0]["commands"]
	# the border polyline is the closed 4-corner-ish outline drawn at width 2
	assert any(command["op"] == "attr" and command.get("width") == 2.0
	           for command in commands), "no width-2 border stroke in the grid layer"
	assert plot.layers == ["grid"], "there must be exactly one grid layer, border included"


# ---------------------------------------------------------------------------
# colour handling (§1.1, §2.6)
# ---------------------------------------------------------------------------

def colour_taking_calls(plot):
	'''
	Every draw method that takes a colour — ``grid()`` included, not only the
	``add*`` family.
	'''
	return [
		("grid", lambda **kwargs: plot.grid(**kwargs)),
		("addRegionOutline", lambda **kwargs: plot.addRegionOutline(circle(), **kwargs)),
		("addMoc", lambda **kwargs: plot.addMoc(circle().astObject, **kwargs)),
		("addPoints", lambda **kwargs: plot.addPoints(ra=[30.0], dec=[45.0], **kwargs)),
		("addText", lambda **kwargs: plot.addText("M51", [30.0, 45.0], **kwargs)),
		("addCompass", lambda **kwargs: plot.addCompass(**kwargs)),
		("addScaleBar", lambda **kwargs: plot.addScaleBar(30 * u.arcmin, **kwargs)),
	]


@pytest.mark.parametrize("index", range(7))
def test_colour_and_color_contradiction_raises(plot, index):
	'''
	§1.1: both given and different is an error, never a silent preference —
	v1's ``if color: colour = color`` silently favoured the synonym.
	'''
	name, call = colour_taking_calls(plot)[index]
	with pytest.raises(ValueError, match="synonyms"):
		call(colour="red", color="blue")


@pytest.mark.parametrize("index", range(7))
def test_color_synonym_alone_works(plot, index):
	'''
	The None-sentinel class rule in action: ``color=`` alone must not spuriously
	conflict with the default, which a literal default would have caused.
	'''
	name, call = colour_taking_calls(plot)[index]
	call(color="red")


@pytest.mark.parametrize("index", range(7))
def test_matching_colour_and_color_is_accepted(plot, index):
	name, call = colour_taking_calls(plot)[index]
	call(colour="red", color="red")


def test_default_colours_match_the_frozen_table(plot):
	''' §2.6 resolution table, checked through the emitted stream. '''
	def first_colour(layer_id):
		layer = [item for item in plot.document()["layers"] if item["id"] == layer_id][0]
		for command in layer["commands"]:
			if command["op"] == "attr" and "colour" in command:
				return command["colour"]
		return None

	plot.grid()
	assert first_colour("grid") == 11                       # grid grey

	plot.addRegionOutline(circle(), name="r")
	assert first_colour("r") == "#4a7f7b"                   # cornish teal, inlined as CSS

	plot.addPoints(ra=[30.0], dec=[45.0], name="p")
	assert first_colour("p") == 4                           # palette blue


def test_unknown_colour_raises(plot):
	with pytest.raises(ValueError):
		plot.addRegionOutline(circle(), colour="octarine")


def test_out_of_range_palette_index_raises(plot):
	''' A USER-supplied index is validated up front, unlike a stored one in a stream. '''
	with pytest.raises(ValueError):
		plot.addPoints(ra=[30.0], dec=[45.0], colour=99)


# ---------------------------------------------------------------------------
# individual draw methods
# ---------------------------------------------------------------------------

def test_add_region_outline_accepts_raw_ast_regions(plot):
	plot.addRegionOutline(circle().astObject)


def test_add_region_outline_rejects_other_types(plot):
	with pytest.raises(TypeError):
		plot.addRegionOutline("a circle")


def test_moc_layer_has_its_own_kind(plot):
	from cornish import ASTMoc
	moc = circle().toMoc(max_order=6) if hasattr(circle(), "toMoc") else None
	if moc is None:
		pytest.skip("toMoc unavailable")
	layer_id = plot.addMoc(moc)
	kinds = {layer["id"]: layer["kind"] for layer in plot.document()["layers"]}
	assert kinds[layer_id] == "moc"


def test_add_points_input_forms(plot):
	''' Every accepted form is the SPEC-04 bridge contract. '''
	plot.addPoints(ra=[30.0, 30.5], dec=[45.0, 45.2], name="arrays")
	plot.addPoints(ra=30.0, dec=45.0, name="scalars")
	plot.addPoints(points=SkyCoord(ra=[30.0] * u.deg, dec=[45.0] * u.deg), name="skycoord")
	plot.addPoints(points=[[30.0, 45.0]], name="bare")


def test_add_points_emits_a_single_mark_command(plot):
	''' §2.4: ONE mark((2,n)) call, not one per point — v1 looped. '''
	layer_id = plot.addPoints(ra=[30.0, 30.2, 30.4], dec=[45.0, 45.1, 45.2])
	layer = [item for item in plot.document()["layers"] if item["id"] == layer_id][0]
	marks = [command for command in layer["commands"] if command["op"] == "mark"]
	assert len(marks) == 1
	assert len(marks[0]["x"]) == 3


@pytest.mark.parametrize("kwargs", [
	{},                                                    # nothing at all
	{"ra": [30.0]},                                        # dec missing
	{"dec": [45.0]},                                       # ra missing
	{"points": [[30.0, 45.0]], "ra": [30.0], "dec": [45.0]},   # both forms
	{"points": []},                                        # empty
])
def test_add_points_argument_errors(plot, kwargs):
	with pytest.raises(ValueError):
		plot.addPoints(**kwargs)


def test_marker_accepts_enum_name_and_code(plot):
	plot.addPoints(ra=[30.0], dec=[45.0], marker=Marker.STAR, name="a")
	plot.addPoints(ra=[30.0], dec=[45.0], marker="triangle down", name="b")
	plot.addPoints(ra=[30.0], dec=[45.0], marker=7, name="c")
	with pytest.raises(ValueError):
		plot.addPoints(ra=[30.0], dec=[45.0], marker="squiggle", name="d")


def test_add_text_justification_grammar(plot):
	plot.addText("M51", [30.0, 45.0], just="TL")
	with pytest.raises(ValueError, match="justification"):
		plot.addText("M51", [30.0, 45.0], just="ZZ")
	with pytest.raises(ValueError, match="justification"):
		plot.addText("M51", [30.0, 45.0], just="C")


def test_add_text_up_vector_must_be_non_zero(plot):
	with pytest.raises(ValueError, match="non-zero"):
		plot.addText("M51", [30.0, 45.0], up=(0.0, 0.0))


def test_add_text_requires_a_string(plot):
	with pytest.raises(TypeError):
		plot.addText(51, [30.0, 45.0])


def test_line_style_validation(plot):
	plot.addRegionOutline(circle(), style=LineStyle.DASHED)
	plot.addRegionOutline(circle(), style=3, name="s3")
	with pytest.raises(ValueError):
		plot.addRegionOutline(circle(), style=9, name="s9")


# ---------------------------------------------------------------------------
# compass and scale bar (§2.4, R25)
# ---------------------------------------------------------------------------

def test_compass_has_two_arms_and_two_labels(plot):
	plot.addCompass()
	layer = [item for item in plot.document()["layers"] if item["id"] == "compass"][0]
	lines = [command for command in layer["commands"] if command["op"] == "line"]
	texts = [command for command in layer["commands"] if command["op"] == "text"]
	assert len(lines) >= 2
	assert {command["text"] for command in texts} == {"N", "E"}


def test_compass_arms_use_true_bearings(plot):
	'''
	§2.4: the arms are geodesics whose INITIAL BEARINGS are north and east,
	from ``Frame.offset2``. A naive coordinate offset would not point east
	anywhere but the equator — at dec 45 the error is large.
	'''
	frame = bridge.current_frame_of(plot.astPlot)
	anchor = np.array([30.0 * math.pi / 180, 45.0 * math.pi / 180])
	length = 0.01
	_, north = frame.offset2(anchor, 0.0, length)
	_, east = frame.offset2(anchor, math.pi / 2, length)
	# due north: longitude unchanged, latitude increased by exactly the length
	assert north[0] == pytest.approx(anchor[0], abs=1e-12)
	assert north[1] - anchor[1] == pytest.approx(length, rel=1e-9)
	# due east: a naive alpha + length would be wrong by 1/cos(dec) ~ 1.41x here
	naive = anchor[0] + length
	assert east[0] != pytest.approx(naive, rel=1e-6)
	assert float(frame.distance(anchor, east)) == pytest.approx(length, rel=1e-9)


def test_compass_is_refused_where_it_cannot_be_placed():
	'''
	§2.4: when the arms leave the plot at both candidate corners the compass is
	refused rather than drawn somewhere wrong. An all-sky projection's corners
	lie outside the sky ellipse entirely.
	'''
	with pytest.raises(ValueError, match="compass could not be placed"):
		SkyPlot(all_sky=True).addCompass()


def test_compass_accepts_an_explicit_position(plot):
	assert plot.addCompass(position=[30.0, 45.0]) == "compass"


def test_scale_bar_length_is_exact_on_the_great_circle(plot):
	'''
	§2.4: astPlot.curve draws GREAT CIRCLES, so the endpoints come from the
	equal-declination great-circle inversion rather than the parallel arc. The
	spec asks for 1 part in 1e6; the closed form delivers far better.
	'''
	frame = bridge.current_frame_of(plot.astPlot)
	requested = bridge.to_frame_distance(30 * u.arcmin, frame)
	start, finish, _ = plot.scaleBarEndpoints(requested)
	assert float(frame.distance(start, finish)) == pytest.approx(requested, rel=1e-6)


def test_scale_bar_layer_contents(plot):
	plot.addScaleBar(30 * u.arcmin)
	layer = [item for item in plot.document()["layers"] if item["id"] == "scalebar"][0]
	assert any(command["op"] == "line" for command in layer["commands"])
	assert any(command["op"] == "text" for command in layer["commands"])


def test_scale_bar_longer_than_the_field_raises(plot):
	with pytest.raises(ValueError, match="longer than the plotted field"):
		plot.addScaleBar(10 * u.deg)


def test_scale_bar_at_extreme_latitude_raises():
	''' A 5° bar cannot fit at δ = 89.9°: it would exceed the small circle there. '''
	plot = SkyPlot(extent=circle(centre=(0.0, 89.9), radius=0.1))
	frame = bridge.current_frame_of(plot.astPlot)
	with pytest.raises(ValueError):
		plot.scaleBarEndpoints(bridge.to_frame_distance(5 * u.deg, frame))


def test_scale_bar_requires_a_positive_angular_quantity(plot):
	with pytest.raises(TypeError):
		plot.addScaleBar(30)
	with pytest.raises(ValueError):
		plot.addScaleBar(-1 * u.deg)


def test_scale_bar_label_defaults_to_the_formatted_length(plot):
	plot.addScaleBar(30 * u.arcmin)
	layer = [item for item in plot.document()["layers"] if item["id"] == "scalebar"][0]
	texts = [command["text"] for command in layer["commands"] if command["op"] == "text"]
	assert any("30" in text for text in texts)
	plot.addScaleBar(30 * u.arcmin, label="half a degree")
	layer = [item for item in plot.document()["layers"] if item["id"] == "scalebar"][0]
	texts = [command["text"] for command in layer["commands"] if command["op"] == "text"]
	assert "half a degree" in texts


# ---------------------------------------------------------------------------
# addImage (§2.4)
# ---------------------------------------------------------------------------

def test_add_image_requires_the_plots_own_grid(plot):
	'''
	The frozen stage-0 constraint. The error names stage 3 as the path for
	regridding, so a user meeting it learns what would make it work.
	'''
	plot.addImage(np.zeros((100, 100)))
	with pytest.raises(ValueError, match="stage-3|stage 3"):
		plot.addImage(np.zeros((50, 50)), name="wrong")


def test_add_image_emits_a_data_uri(plot):
	layer_id = plot.addImage(np.random.RandomState(0).rand(100, 100))
	layer = [item for item in plot.document()["layers"] if item["id"] == layer_id][0]
	command = layer["commands"][0]
	assert command["op"] == "image"
	assert command["href"].startswith("data:image/png;base64,")


# ---------------------------------------------------------------------------
# fromFITSHeader (§2.4)
# ---------------------------------------------------------------------------

def test_from_fits_header_uses_the_naxis_ratio():
	plot = SkyPlot.fromFITSHeader(dict(SQUARE_HEADER))
	assert plot.document()["canvas"]["aspect"] == 32 / 64
	assert plot.layers == ["grid"]


def test_from_fits_header_missing_naxis_raises_incomplete_header():
	''' The existing cornish exception type, not a bare ValueError. '''
	header = dict(SQUARE_HEADER)
	del header["NAXIS1"]
	with pytest.raises(IncompleteHeader):
		SkyPlot.fromFITSHeader(header)


def test_from_fits_header_refuses_non_square_pixels_in_cdelt_form():
	header = dict(SQUARE_HEADER)
	header["CDELT2"] = 0.02
	with pytest.raises(ValueError, match="non-square pixels"):
		SkyPlot.fromFITSHeader(header)


def test_from_fits_header_refuses_non_square_pixels_in_cd_matrix_form():
	'''
	The detection is deliberately card-form INDEPENDENT: a CD/PC header may
	carry no CDELT cards at all, so reading cards would silently miss exactly
	the case that matters. The check transforms pixels through the frameset and
	measures with the frame's own distance().
	'''
	header = {key: value for key, value in SQUARE_HEADER.items()
	          if key not in ("CDELT1", "CDELT2")}
	header.update({"CD1_1": -0.01, "CD1_2": 0.0, "CD2_1": 0.0, "CD2_2": 0.02})
	with pytest.raises(ValueError, match="non-square pixels"):
		SkyPlot.fromFITSHeader(header)


def test_from_fits_header_error_names_both_scales():
	header = dict(SQUARE_HEADER)
	header["CDELT2"] = 0.02
	with pytest.raises(ValueError) as raised:
		SkyPlot.fromFITSHeader(header)
	assert "0.01" in str(raised.value) and "0.02" in str(raised.value)


def test_from_fits_header_still_validates_parameters():
	''' Bypassing the synthetic-header builder must not bypass validation. '''
	with pytest.raises(ValueError):
		SkyPlot.fromFITSHeader(dict(SQUARE_HEADER), char_height=0)
	with pytest.raises(ValueError):
		SkyPlot.fromFITSHeader(dict(SQUARE_HEADER), margin=0.9)


# ---------------------------------------------------------------------------
# output surfaces
# ---------------------------------------------------------------------------

def test_document_round_trips_through_the_validator(plot):
	plot.addRegionOutline(circle())
	plot.addPoints(ra=[30.0], dec=[45.0])
	validate_document(plot.document())


def test_envelopes_carry_the_canvas_on_the_first_only(plot):
	plot.addRegionOutline(circle())
	envelopes = plot.envelopes()
	assert "canvas" in envelopes[0]
	assert all("canvas" not in envelope for envelope in envelopes[1:])
	for envelope in envelopes:
		validate_envelope(envelope)


def test_envelopes_of_an_empty_plot_is_the_canvas_only_form():
	envelopes = SkyPlot(extent=circle(), grid=False).envelopes()
	assert len(envelopes) == 1
	assert "layer" not in envelopes[0] and "canvas" in envelopes[0]
	validate_envelope(envelopes[0])


def test_canvas_palette_is_the_base_palette_not_the_grown_one(plot):
	'''
	D18: colours grown by ColToInt are inlined into their layers as literal CSS,
	so the canvas stays a themeable table of fixed length.
	'''
	plot.addRegionOutline(circle())                 # grows the palette with the teal
	assert plot.document()["canvas"]["palette"] == list(DEFAULT_PALETTE)


def test_save_suffix_dispatch(plot, tmp_path):
	svg_path = tmp_path / "plot.svg"
	json_path = tmp_path / "plot.json"
	plot.save(svg_path)
	plot.save(json_path)
	assert svg_path.read_text().startswith("<svg")
	assert json_path.read_text().startswith("{")
	with pytest.raises(ValueError, match="'.svg' and '.json'"):
		plot.save(tmp_path / "plot.png")


def test_saved_json_is_canonical(plot, tmp_path):
	from cornish.plot.envelope import canonical_json
	path = tmp_path / "plot.json"
	plot.save(path)
	assert path.read_text() == canonical_json(plot.document())


def test_repr_mimebundle_returns_svg(plot):
	bundle = plot._repr_mimebundle_()
	assert set(bundle) == {"image/svg+xml"}
	assert bundle["image/svg+xml"].startswith("<svg")


def test_svg_renders_without_matplotlib_installed(plot):
	'''
	The stage-0 exit criterion. matplotlib is never imported by the render path
	— the purity gate in test_render_svg.py enforces that structurally; here we
	confirm the SVG actually comes out.
	'''
	assert plot.svg().startswith("<svg")


def test_to_matplotlib_replays_the_document(plot):
	matplotlib = pytest.importorskip("matplotlib")
	matplotlib.use("Agg", force=True)
	plot.addPoints(ra=[30.0], dec=[45.0])
	figure = plot.to_matplotlib()
	assert figure.axes
	matplotlib.pyplot.close(figure)


# ---------------------------------------------------------------------------
# the out-of-viewport guard (§4.4)
# ---------------------------------------------------------------------------

def test_text_outside_the_viewport_warns(plot, caplog):
	'''
	§4.4: the viewport crop stays — a dynamic margin would make the viewport a
	function of content, breaking byte-stable rendering — but a cropped label is
	never silent. The warning names the layer and the text.
	'''
	far = plot._sky_at_gbox(3.0, 3.0)
	if far is None:
		pytest.skip("this projection has no sky at that graphics position")
	degrees = bridge.from_frame_units(far, plot.astPlot)
	with caplog.at_level(logging.WARNING, logger="cornish"):
		plot.addText("way out there", [float(degrees[0]), float(degrees[1])], name="far")
	messages = [record.getMessage() for record in caplog.records]
	assert any("far" in message and "way out there" in message for message in messages)


def test_text_inside_the_viewport_stays_silent(plot, caplog):
	with caplog.at_level(logging.WARNING, logger="cornish"):
		plot.addText("M51", [30.0, 45.0], name="near")
	assert not [record for record in caplog.records if "outside the gbox" in record.getMessage()]
