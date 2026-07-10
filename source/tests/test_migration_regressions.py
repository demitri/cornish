
'''
SPEC-04 migration regression tests (SPEC-04A §9 T6/T7) plus the
conversions-outside-bridge gate.

Each T7 test pins one flagged behavior change (Δ) from the SPEC-04A §5
migration map; the test names follow the spec's table.
'''

import os
import glob
import pickle
import re
import warnings

import numpy as np
import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord
import starlink.Ast as Ast

from cornish import (ASTBox, ASTCircle, ASTPolygon, ASTRegion, ASTCompoundRegion,
                     ASTFrameSet, ASTICRSFrame, ASTMoc)
from cornish import ASTFrame as CornishFrame
from cornish.mapping.frame.sky_frame import ASTSkyFrame
from cornish.channel import ASTFITSChannel
from cornish import _pyast_bridge as bridge
from cornish.exc import IncompleteHeader, CornishError


GOLDEN_TAN_HEADER = {
	"NAXIS": 2, "NAXIS1": 512, "NAXIS2": 512,
	"CRPIX1": 256.5, "CRPIX2": 256.5, "CRVAL1": 30.0, "CRVAL2": 45.0,
	"CDELT1": -0.001, "CDELT2": 0.001,
	"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"}


@pytest.fixture
def tan_frameset():
	return ASTFrameSet.fromFITSHeader(fits_header=GOLDEN_TAN_HEADER)

@pytest.fixture
def sky_to_galactic_frameset():
	''' an all-sky frameset: base ICRS, current galactic '''
	icrs = Ast.SkyFrame("System=ICRS")
	galactic = Ast.SkyFrame("System=Galactic")
	fs = Ast.Frame.convert(icrs, galactic)
	assert fs is not None
	return ASTFrameSet(ast_object=fs)


# ---------------------------------------------------------------------------
# T7 — one test per flagged behavior change
# ---------------------------------------------------------------------------

def test_pointInRegion_galactic_skycoord():
	''' M11: a galactic SkyCoord at a known ICRS-region point (was: AttributeError) '''
	circle = ASTCircle(center=[30, 45], radius=2.0)
	galactic_point = SkyCoord(ra=30 * u.deg, dec=45 * u.deg, frame="icrs").galactic
	assert circle.pointInRegion(galactic_point) is True


def test_pointInRegion_basic_frame_no_deg2rad():
	''' M11: GRID-frame circle no longer deg2rad's the test point (was: False) '''
	frame = CornishFrame(naxes=2)
	circle = ASTCircle(frame=frame, center=[50, 50], radius=10.0)
	assert circle.pointInRegion([50, 50]) is True
	assert circle.pointInRegion([65, 50]) is False


def test_convertRaDec_quantity_returns_degrees(sky_to_galactic_frameset):
	''' M25/N2: Quantity callers were receiving radians labeled degrees '''
	ra_out, dec_out = sky_to_galactic_frameset.convertRaDec(266.405 * u.deg, -28.936 * u.deg)
	# (266.405, -28.936) ICRS is the galactic centre: l~0, b~0
	assert isinstance(ra_out, u.Quantity) and ra_out.unit == u.deg
	# radians mislabeled as degrees would put the galactic centre at |l| < 0.02 "deg" only
	# by accident; check the true value in degrees
	l = ra_out.to_value(u.deg) % 360
	assert min(l, 360 - l) < 0.1
	assert abs(dec_out.to_value(u.deg)) < 0.1


def test_convertRaDec_inverse_nonsky_units(tan_frameset):
	''' M25 round-2 P1: Quantity in, non-sky destination -> plain values, never a Quantity '''
	ra_out, dec_out = tan_frameset.convertRaDec(30.0 * u.deg, 45.0 * u.deg, forward=False)
	assert not isinstance(ra_out, u.Quantity)
	assert not isinstance(dec_out, u.Quantity)
	assert ra_out == pytest.approx(256.5, abs=0.5)
	assert dec_out == pytest.approx(256.5, abs=0.5)


@pytest.mark.parametrize("method", ["pix2world", "world2pix", "convertPoints"])
def test_frameset_square_input_follows_D3(tan_frameset, method):
	''' M24/M26/M27: the (2,2) square case reads as PAIRS; parallel_axes=True restores parallel '''
	fs = tan_frameset
	if method == "pix2world":
		call = fs.pix2world
		p1, p2 = [100.0, 200.0], [300.0, 400.0]
	elif method == "world2pix":
		call = fs.world2pix
		p1, p2 = [30.1, 45.1], [29.9, 44.9]
	else:
		call = fs.convertPoints
		p1, p2 = [100.0, 200.0], [300.0, 400.0]

	single1 = call(p1)
	single2 = call(p2)
	square = np.array([p1, p2])

	pairs_reading = call(square)
	np.testing.assert_allclose(pairs_reading, np.array([single1, single2]), atol=1e-9)

	parallel_reading = call(square, parallel_axes=True)
	expected1 = call([square[0][0], square[1][0]])
	expected2 = call([square[0][1], square[1][1]])
	np.testing.assert_allclose(parallel_reading, np.array([expected1, expected2]), atol=1e-9)


def test_pix2world_documented_parallel_form(tan_frameset):
	''' M26/N5: the documented (2, n) input form works (was: crash for n != 2) '''
	xs = np.array([10.0, 100.0, 200.0, 300.0, 400.0])
	ys = np.array([20.0, 110.0, 210.0, 310.0, 410.0])
	out = tan_frameset.pix2world(np.stack([xs, ys]))
	assert out.shape == (5, 2)
	# spot-check one point against the single-point form
	single = tan_frameset.pix2world([xs[2], ys[2]])
	np.testing.assert_allclose(out[2], single, atol=1e-9)
	# output is degrees near the reference point
	assert 25 < out[2][0] < 35 and 40 < out[2][1] < 50


def test_world2pix_pairs_and_skycoord(tan_frameset):
	''' M27: pairs form + galactic SkyCoord single point '''
	pairs = np.array([[30.1, 45.1], [30.0, 45.0], [29.9, 44.9], [30.05, 44.95], [29.95, 45.05]])
	out = tan_frameset.world2pix(pairs)
	assert out.shape == (5, 2)

	galactic_centre_point = SkyCoord(ra=30 * u.deg, dec=45 * u.deg, frame="icrs").galactic
	single = tan_frameset.world2pix(galactic_centre_point)
	assert single.shape == (2,)
	np.testing.assert_allclose(single, [256.5, 256.5], atol=0.5)


def test_angle_basic_frame_native():
	'''
	M22/N4: basic-frame angle no longer deg2rad-mangles its inputs.
	(AST's astAngle is SIGNED on 2D basic frames — positive anticlockwise —
	so the right angle here pins to magnitude.)
	'''
	frame = CornishFrame(naxes=2)
	angle = frame.angle(vertex=[0, 0], points=([1, 0], [0, 1]))
	assert abs(angle.to_value(u.deg)) == pytest.approx(90.0)


def test_distance_basic_frame_quantity_raises():
	''' M21/N8: basic-frame distance with Quantity input raises instead of silent unit-drop '''
	frame = CornishFrame(naxes=2)
	with pytest.raises(ValueError):
		frame.distance([0, 0] * u.deg, [3, 4])
	# bare values still work natively
	assert float(frame.distance([0, 0], [3, 4])) == pytest.approx(5.0)


def test_offset_basic_frame_float_ok_quantity_raises():
	''' M23/N7: basic-frame offset takes floats (new), rejects Quantities (was silent .value drop) '''
	frame = CornishFrame(naxes=2)
	out = frame.offsetAlongGeodesicCurve([0, 0], [10, 0], 5.0)
	np.testing.assert_allclose(np.asarray(out, dtype=float), [5.0, 0.0], atol=1e-9)
	with pytest.raises(ValueError):
		frame.offsetAlongGeodesicCurve([0, 0], [10, 0], 5.0 * u.deg)


def test_uncertainty_reaches_ast():
	''' M8/N1: uncertainty is delivered to AST for the first time '''
	c1 = ASTCircle(center=[30.0, 45.0], radius=1.0, uncertainty=0.05 * u.deg)
	c2 = ASTCircle(center=[30.005, 45.0], radius=1.0, uncertainty=0.05 * u.deg)
	assert c1.isIdenticalTo(c2) is True

	c3 = ASTCircle(center=[30.0, 45.0], radius=1.0, uncertainty=1e-7 * u.deg)
	c4 = ASTCircle(center=[30.005, 45.0], radius=1.0, uncertainty=1e-7 * u.deg)
	assert c3.isIdenticalTo(c4) is False


def test_uncertainty_default_is_arcsec_quantity():
	''' §8: the omitted-parameter default materializes and is delivered (Unc block present) '''
	circle = ASTCircle(center=[30, 45], radius=1.0)
	block = bridge.dump_value(circle.astObject, "Unc")
	assert block.strip().startswith("Begin")
	assert isinstance(circle.uncertainty, u.Quantity)
	assert circle.uncertainty == 1 * u.arcsec


def test_uncertainty_three_way_default():
	''' Δ4: omitted -> ~1" Unc delivered; explicit None -> AST default; explicit value -> delivered '''
	from cornish.channel.channel_io import ListSource

	# omitted: materialized 1 arcsec
	c_omitted = ASTCircle(center=[10.0, 20.0], radius=1.0)
	block = bridge.dump_value(c_omitted.astObject, "Unc")
	unc_region = Ast.Channel(ListSource(block), None).read()
	radius_deg = np.rad2deg(unc_region.circlepars()[1])
	# AST re-derives the stored uncertainty region (~6% at this dec); never compare exactly
	assert radius_deg == pytest.approx(1 / 3600, rel=0.15)

	# explicit None: nothing delivered; AST's internal default governs
	c_none = ASTCircle(center=[10.0, 20.0], radius=1.0, uncertainty=None)
	with pytest.raises(KeyError):
		bridge.dump_value(c_none.astObject, "Unc")
	assert c_none.uncertainty is None

	# explicit value: delivered
	c_explicit = ASTCircle(center=[10.0, 20.0], radius=1.0, uncertainty=0.01 * u.deg)
	block = bridge.dump_value(c_explicit.astObject, "Unc")
	unc_region = Ast.Channel(ListSource(block), None).read()
	assert np.rad2deg(unc_region.circlepars()[1]) == pytest.approx(0.01, rel=0.15)


def test_constructors_no_positional_none_unc():
	''' V18: uncertainty=None means the positional argument is OMITTED, not passed as None '''
	box = ASTBox.fromCorners(frame=ASTICRSFrame(), corners=([10, 20], [11, 21]), uncertainty=None)
	with pytest.raises(KeyError):
		bridge.dump_value(box.astObject, "Unc")
	circle = ASTCircle(center=[30, 45], radius=1.0, uncertainty=None)
	with pytest.raises(KeyError):
		bridge.dump_value(circle.astObject, "Unc")


def test_uncertainty_setter_raises():
	''' §8: the setter remains a loud refusal until a real astSetUnc exposure exists '''
	circle = ASTCircle(center=[30, 45], radius=1.0)
	with pytest.raises(NotImplementedError):
		circle.uncertainty = 0.1 * u.deg


def test_fromFITSHeader_unit_honest_fallback():
	''' M19: the rebuilt (unit-honest) mesh fallback reproduces the fast-path footprint '''
	channel = ASTFITSChannel(header=GOLDEN_TAN_HEADER)
	fast = ASTPolygon._fromParsedWCS(channel.frameSet, channel.dimensions)
	fallback = ASTPolygon._fromParsedWCS(channel.frameSet, channel.dimensions, _force_fallback=True)

	assert fast.overlapType(fallback).name == "IDENTICAL"
	# >= 99.9% mutual edge membership
	mesh = fallback.boundaryPointMesh(npoints=1000)
	inside = sum(fast.pointInRegion(p) for p in mesh)
	assert inside / len(mesh) >= 0.999
	mesh = fast.boundaryPointMesh(npoints=1000)
	inside = sum(fallback.pointInRegion(p) for p in mesh)
	assert inside / len(mesh) >= 0.999


def test_fromFITSHeader_wcs_retention():
	''' SPEC-04 §2: the returned region carries its pixel<->world frame set '''
	polygon = ASTPolygon.fromFITSHeader(GOLDEN_TAN_HEADER)
	assert isinstance(polygon.wcs, ASTFrameSet)
	# and propagates through downsize / toPolygon / regionWithMapping
	assert polygon.downsize(maxerr=1 * u.arcsec, maxvert=50).wcs is polygon.wcs
	assert polygon.toPolygon().wcs is polygon.wcs


def test_region_fromFITSHeader_dispatcher():
	''' SPEC-04 §3.4 / B8: the README front door works again '''
	region = ASTRegion.fromFITSHeader(GOLDEN_TAN_HEADER)
	assert isinstance(region, ASTPolygon)
	assert region.pointInRegion([30.0, 45.0])


def test_fromFITSHeader_missing_naxis_incomplete_header():
	''' SPEC-04 §3.5: headers without NAXIS1/2 raise IncompleteHeader naming the card '''
	header = dict(GOLDEN_TAN_HEADER)
	del header["NAXIS1"]
	with pytest.raises(IncompleteHeader, match="NAXIS1"):
		ASTPolygon.fromFITSHeader(header)


def test_polygon_ctor_frameset_degrees(tan_frameset):
	''' M17: ASTPolygon(frame=<FITS frameset>, points=<degrees>) — the honest replacement of the accident chain '''
	points = np.array([[30.2, 45.2], [30.2, 44.8], [29.8, 44.8], [29.8, 45.2]])
	polygon = ASTPolygon(frame=tan_frameset, points=points)
	assert polygon.pointInRegion([30.0, 45.0])
	assert not polygon.pointInRegion([31.0, 46.0])


def test_box_fromCorners_frameset_frame(tan_frameset):
	''' M3: a FrameSet as the frame argument works (was: latent NameError path) '''
	box = ASTBox.fromCorners(frame=tan_frameset, corners=([30.1, 45.1], [29.9, 44.9]))
	assert box.pointInRegion([30.0, 45.0])


def test_box_fromCentreAndCorner_frameset_frame(tan_frameset):
	''' M2: same for fromCentreAndCorner — a separate copy of the latent bug '''
	box = ASTBox.fromCentreAndCorner(frame=tan_frameset, centre=[30.0, 45.0], corner=[30.1, 45.1])
	assert box.pointInRegion([30.0, 45.0])
	assert not box.pointInRegion([30.5, 45.5])


def test_convertPoints_mixed_frameset(tan_frameset):
	''' M24: pixel<->sky frameset conversion, both directions, pinned against pix2world/world2pix '''
	pixel_points = np.array([[100.0, 200.0], [300.0, 400.0], [256.5, 256.5]])
	forward = tan_frameset.convertPoints(pixel_points)
	np.testing.assert_allclose(forward, tan_frameset.pix2world(pixel_points), atol=1e-12)
	# and values are sane degrees (not deg2rad'd pixels)
	assert np.all((forward[:, 0] > 25) & (forward[:, 0] < 35))

	back = tan_frameset.convertPoints(forward, forward=False)
	np.testing.assert_allclose(back, pixel_points, atol=1e-6)
	np.testing.assert_allclose(back, tan_frameset.world2pix(forward), atol=1e-12)


def test_convertPoints_skycoord_list(sky_to_galactic_frameset):
	''' M24 convenience: mixed-frame SkyCoords in -> destination-system SkyCoords out (never ICRS-by-fiat) '''
	sc_icrs = SkyCoord(ra=266.405 * u.deg, dec=-28.936 * u.deg, frame="icrs")
	sc_galactic = SkyCoord(l=90 * u.deg, b=0 * u.deg, frame="galactic")
	result = sky_to_galactic_frameset.convertPoints([sc_icrs, sc_galactic])
	assert isinstance(result, list) and len(result) == 2
	for sc in result:
		assert sc.frame.name == "galactic"
	# the galactic centre in galactic coordinates
	assert min(result[0].l.deg % 360, 360 - result[0].l.deg % 360) < 0.1
	assert abs(result[0].b.deg) < 0.1
	# a galactic input into a galactic destination round-trips
	assert result[1].l.deg == pytest.approx(90.0, abs=1e-6)
	assert result[1].b.deg == pytest.approx(0.0, abs=1e-6)

	# scalar SkyCoord in -> scalar SkyCoord out; array in -> array out
	scalar = sky_to_galactic_frameset.convertPoints(sc_icrs)
	assert scalar.isscalar
	array_in = SkyCoord(ra=[10, 20] * u.deg, dec=[30, 40] * u.deg)
	array_out = sky_to_galactic_frameset.convertPoints(array_in)
	assert not array_out.isscalar and len(array_out) == 2


def test_convertPoints_skycoord_unmappable_destination():
	''' M24: SkyCoord-form input with an AZEL-current destination raises the §2.3 refusal '''
	icrs = Ast.SkyFrame("System=ICRS")
	azel = Ast.SkyFrame("System=AZEL")
	fs = Ast.FrameSet(icrs)
	fs.addframe(1, Ast.UnitMap(2), azel)
	frameset = ASTFrameSet(ast_object=fs)
	sc = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
	with pytest.raises(ValueError, match="no exact astropy equivalent"):
		frameset.convertPoints(sc)


def test_convert_deprecated_alias(tan_frameset):
	''' M28: convert() is a deprecated alias for convertPoints() '''
	with pytest.warns(DeprecationWarning):
		out = tan_frameset.convert([100.0, 200.0])
	np.testing.assert_allclose(out, tan_frameset.pix2world([100.0, 200.0]), atol=1e-12)


def test_box_area_notimplemented_raises():
	''' M7: the previously-unraised NotImplementedError now raises '''
	frame = CornishFrame(naxes=2) # default Frame: Cartesian system, no GRID domain
	box = ASTBox.fromCorners(frame=frame, corners=([0, 0], [10, 10]))
	with pytest.raises(NotImplementedError):
		box.area


def test_box_area_grid_frame():
	''' M7 sibling: the Cartesian GRID branch returns the full box area '''
	frame = CornishFrame(naxes=2)
	frame.domain = "GRID"
	box = ASTBox.fromCorners(frame=frame, corners=([0, 0], [10, 20]))
	assert box.area.to_value(u.pixel * u.pixel) == pytest.approx(200.0)


def test_string_points_bad_input_valueerror():
	''' M17: string parse failure is a chained ValueError, not a bare Exception '''
	with pytest.raises(ValueError) as excinfo:
		bridge.parse_points_string("not points")
	assert excinfo.value.__cause__ is not None
	with pytest.raises(ValueError):
		ASTPolygon(frame=ASTICRSFrame(), points="not points")


def test_fits_channel_boundingPolygon_delegates():
	''' M31: was B9-dead; now a sky polygon identical to ASTPolygon.fromFITSHeader '''
	channel = ASTFITSChannel(header=GOLDEN_TAN_HEADER)
	from_channel = channel.boundingPolygon()
	from_header = ASTPolygon.fromFITSHeader(GOLDEN_TAN_HEADER)
	assert from_channel.overlapType(from_header).name == "IDENTICAL"
	assert bridge.is_sky(from_channel.astObject)


def test_fits_channel_ast_object_constructor():
	''' review finding: ASTFITSChannel(ast_object=...) previously always raised NotImplementedError '''
	raw = Ast.FitsChan()
	# an EMPTY FitsChan is falsy (card-count length), so this also pins the
	# is-not-None fix: truthiness silently discarded the provided object
	channel = ASTFITSChannel(ast_object=raw)
	assert channel.astObject is raw


def test_fits_channel_boundingCircle():
	''' M32: was B9-dead; the circle contains all four field corners, radius AST-exact '''
	channel = ASTFITSChannel(header=GOLDEN_TAN_HEADER)
	circle = channel.boundingCircle()

	frameset = ASTFrameSet.fromFITSHeader(fits_header=GOLDEN_TAN_HEADER)
	dims = [512, 512]
	pixel_corners = np.array([[0.5, 0.5], [0.5, dims[1] + 0.5],
	                          [dims[0] + 0.5, dims[1] + 0.5], [dims[0] + 0.5, 0.5]])
	sky_corners = frameset.pix2world(pixel_corners)
	for corner in sky_corners:
		assert circle.pointInRegion(corner)

	# radius within 1% of getregiondisc on the field polygon
	polygon = ASTPolygon.fromFITSHeader(GOLDEN_TAN_HEADER)
	centre_rad, radius_rad = polygon.astObject.getregiondisc()
	assert circle.radius.to_value(u.deg) == pytest.approx(np.rad2deg(radius_rad), rel=0.01)


def test_downsize_quantity_and_float_deprecation():
	''' M33: Quantity maxerr works as documented; the sky-frame float form warns, values identical '''
	polygon = ASTCircle(center=[30, 45], radius=2.0).toPolygon(npoints=100)
	downsized_quantity = polygon.downsize(maxerr=1 * u.arcsec, maxvert=0)
	with pytest.warns(DeprecationWarning):
		downsized_float = polygon.downsize(maxerr=4.848e-6, maxvert=0)
	np.testing.assert_allclose(downsized_quantity.points, downsized_float.points, atol=1e-9)


def test_downsize_basic_frame_float_native():
	''' M33 basic branch: floats are native units, no warning; Quantity -> ValueError '''
	flat = CornishFrame(naxes=2)
	flat.setUnitForAxis(axis=1, unit="deg")
	flat.setUnitForAxis(axis=2, unit="deg")
	theta = np.linspace(0, 2 * np.pi, 50, endpoint=False)
	points = np.stack([10 + np.cos(theta), 20 + np.sin(theta)]).T
	polygon = ASTPolygon(frame=flat, points=points)
	with warnings.catch_warnings():
		warnings.simplefilter("error")
		downsized = polygon.downsize(maxerr=0.001, maxvert=0)
	assert len(downsized.points) <= len(points)
	with pytest.raises(ValueError):
		polygon.downsize(maxerr=1 * u.arcsec, maxvert=0)


def test_downsize_invalid_maxerr():
	''' M33 validation: pyast accepts all of these silently — the cornish gate is the only gate '''
	polygon = ASTCircle(center=[30, 45], radius=2.0).toPolygon(npoints=50)
	for bad in [float("nan"), Ast.BAD, -1.0]:
		with pytest.raises(ValueError), warnings.catch_warnings():
			# the -1.0 case legitimately passes the sky-frame bare-float
			# deprecation warning on its way to the negative-value ValueError
			warnings.simplefilter("ignore", DeprecationWarning)
			polygon.downsize(maxerr=bad, maxvert=0)
	with pytest.raises(TypeError):
		polygon.downsize(maxerr=True, maxvert=0)
	# zero stays legal: maxvert governs
	with pytest.warns(DeprecationWarning): # bare float on a sky frame
		downsized = polygon.downsize(maxerr=0.0, maxvert=10)
	assert len(downsized.points) <= 10


def test_circle_basic_frame_radius_native():
	''' M9: a basic-frame circle's radius is a bare native float — not rad2deg'd, not labeled deg '''
	circle = ASTCircle(frame=CornishFrame(naxes=2), center=[50, 50], radius=10.0)
	radius = circle.radius
	assert not isinstance(radius, u.Quantity)
	assert radius == pytest.approx(10.0)
	np.testing.assert_allclose(np.asarray(circle.centre), [50.0, 50.0], atol=1e-9)


def test_circle_basic_frame_operators():
	''' M10: basic-frame circle operators work natively; Quantity operand refused '''
	circle = ASTCircle(frame=CornishFrame(naxes=2), center=[50, 50], radius=10.0)
	assert (circle * 2).radius == pytest.approx(20.0)
	assert (circle / 2).radius == pytest.approx(5.0)
	assert (circle + 5.0).radius == pytest.approx(15.0)
	with pytest.raises(ValueError):
		circle + 2 * u.deg


def test_circle_sky_frame_operators():
	''' M10: sky-frame circle operators keep Quantity semantics '''
	circle = ASTCircle(center=[30, 45], radius=2.0)
	assert (circle + 1.0).radius.to_value(u.deg) == pytest.approx(3.0)
	assert (circle + 30 * u.arcmin).radius.to_value(u.deg) == pytest.approx(2.5)
	assert (circle * 2).radius.to_value(u.deg) == pytest.approx(4.0)
	assert (circle / 2).radius.to_value(u.deg) == pytest.approx(1.0)


@pytest.mark.skipif(not pytest.importorskip("matplotlib"), reason="matplotlib required")
class TestAddPoints:
	''' M29: addPoints via a real SkyPlot with mark calls captured through a stub Grf '''

	@pytest.fixture
	def plot_with_captured_marks(self):
		import matplotlib
		matplotlib.use("Agg", force=True)
		from cornish.plot.matplotlib import SkyPlot
		from cornish.ast_object import ASTObject

		# the pyast Plot's own attributes are read-only, so capture mark()
		# calls through a bridge-transparent wrapper (an ASTObject whose
		# .astObject is the real Ast.Plot) swapped onto the SkyPlot instance
		class PlotCapture(ASTObject):
			def __init__(self, real_plot):
				super().__init__(ast_object=real_plot)
				self.marks = []
			def mark(self, point, style):
				self.marks.append((np.asarray(point), style))
			def __getattr__(self, name):
				return getattr(self.astObject, name)

		circle = ASTCircle(center=[30, 45], radius=2.0)
		skyplot = SkyPlot(extent=circle, figsize=(4, 4))
		capture = PlotCapture(skyplot.astPlot)
		skyplot.astPlot = capture
		return skyplot, capture.marks

	def test_addPoints_quantity_array(self, plot_with_captured_marks):
		''' N6: the (n,2)-Quantity branch was dead on arrival (NameError) '''
		skyplot, marks = plot_with_captured_marks
		points = np.array([[30.0, 45.0], [30.1, 45.1], [29.9, 44.9]]) * u.deg
		skyplot.addPoints(points=points)
		assert len(marks) == 3
		np.testing.assert_allclose(marks[0][0], np.deg2rad([30.0, 45.0]), atol=1e-12)

	def test_addPoints_scalar_ra_dec(self, plot_with_captured_marks):
		''' the codex-caught scalar path: scalar ra/dec plots exactly one mark '''
		skyplot, marks = plot_with_captured_marks
		skyplot.addPoints(ra=30 * u.deg, dec=45 * u.deg)
		assert len(marks) == 1
		np.testing.assert_allclose(marks[0][0], np.deg2rad([30.0, 45.0]), atol=1e-12)
		marks.clear()
		skyplot.addPoints(ra=30.0, dec=45.0)
		assert len(marks) == 1
		np.testing.assert_allclose(marks[0][0], np.deg2rad([30.0, 45.0]), atol=1e-12)

	def test_addPoints_skycoord_list(self, plot_with_captured_marks):
		skyplot, marks = plot_with_captured_marks
		galactic = SkyCoord(ra=30 * u.deg, dec=45 * u.deg, frame="icrs").galactic
		skyplot.addPoints(points=[galactic, SkyCoord(ra=30.1 * u.deg, dec=45.1 * u.deg)])
		assert len(marks) == 2
		np.testing.assert_allclose(marks[0][0], np.deg2rad([30.0, 45.0]), atol=1e-10)

	def test_addPoints_scalar_skycoord(self, plot_with_captured_marks):
		''' review finding: a scalar SkyCoord (single-point form P1) must not die in a len() probe '''
		skyplot, marks = plot_with_captured_marks
		skyplot.addPoints(points=SkyCoord(ra=30 * u.deg, dec=45 * u.deg))
		assert len(marks) == 1
		np.testing.assert_allclose(marks[0][0], np.deg2rad([30.0, 45.0]), atol=1e-12)

	def test_addPoints_square_follows_D3(self, plot_with_captured_marks):
		''' Δ7: the (2,2) square case reads as PAIRS at the M29 site too '''
		skyplot, marks = plot_with_captured_marks
		skyplot.addPoints(points=np.array([[30.0, 45.0], [30.1, 45.1]]))
		assert len(marks) == 2
		np.testing.assert_allclose(marks[0][0], np.deg2rad([30.0, 45.0]), atol=1e-12)
		np.testing.assert_allclose(marks[1][0], np.deg2rad([30.1, 45.1]), atol=1e-12)

	def test_addPoints_both_forms_rejected(self, plot_with_captured_marks):
		''' review finding: the old both-ra/dec-and-points guard was vacuous '''
		skyplot, marks = plot_with_captured_marks
		with pytest.raises(ValueError):
			skyplot.addPoints(points=[[30.0, 45.0]], ra=30 * u.deg, dec=45 * u.deg)

	def test_skyplot_nonsky_extent_refused(self):
		''' review finding: a basic-frame extent must be a loud ValueError, not an AttributeError '''
		import matplotlib
		matplotlib.use("Agg", force=True)
		from cornish.plot.matplotlib import SkyPlot
		box = ASTBox.fromCorners(frame=CornishFrame(naxes=2), corners=([0, 0], [10, 10]), uncertainty=None)
		with pytest.raises(ValueError, match="sky frame"):
			SkyPlot(extent=box, figsize=(4, 4))


# ---------------------------------------------------------------------------
# T6 — pickle round-trips (D6/D16)
# ---------------------------------------------------------------------------

class TestT6Pickle:

	def roundtrip(self, obj):
		return pickle.loads(pickle.dumps(obj))

	def test_circle_sky(self):
		circle = ASTCircle(center=[30, 45], radius=2.0, uncertainty=0.01 * u.deg)
		restored = self.roundtrip(circle)
		assert restored.isIdenticalTo(circle)
		assert restored.uncertainty == 0.01 * u.deg          # Python-side state preserved
		bridge.dump_value(restored.astObject, "Unc")         # delivered unc preserved (no KeyError)

	def test_circle_basic(self):
		circle = ASTCircle(frame=CornishFrame(naxes=2), center=[50, 50], radius=10.0)
		restored = self.roundtrip(circle)
		assert restored.radius == pytest.approx(10.0)
		np.testing.assert_allclose(np.asarray(restored.centre), [50.0, 50.0], atol=1e-9)
		assert restored.astObject.overlap(circle.astObject) == 5
		assert not restored.astObject.getregionframe().isaskyframe()

	def test_polygon(self):
		polygon = ASTPolygon(frame=ASTICRSFrame(), points=[[10, 20], [10, 21], [11, 21], [11, 20]])
		restored = self.roundtrip(polygon)
		assert restored.isIdenticalTo(polygon)

	def test_negated_polygon(self):
		'''
		Negation SEMANTICS round-trip. (Verified: AST folds the dump's
		``Negate`` component into the region definition on read, so the
		read-back ``Negated`` attribute reads 0 while membership and overlap
		are exactly preserved — pin the behavior, not the attribute.)
		'''
		polygon = ASTPolygon(frame=ASTICRSFrame(), points=[[10, 20], [10, 21], [11, 21], [11, 20]])
		polygon.negate()
		restored = self.roundtrip(polygon)
		assert restored.isIdenticalTo(polygon)
		assert restored.pointInRegion([10.5, 20.5]) is False # inside the vertices = outside the negated region
		assert restored.pointInRegion([50.0, 50.0]) is True

	def test_box(self):
		box = ASTBox.fromCorners(frame=ASTICRSFrame(), corners=([10, 20], [11, 21]))
		restored = self.roundtrip(box)
		assert restored.isIdenticalTo(box)

	def test_frameset_fits_read(self, tan_frameset):
		restored = self.roundtrip(tan_frameset)
		np.testing.assert_allclose(restored.pix2world([100.0, 200.0]),
		                           tan_frameset.pix2world([100.0, 200.0]), atol=1e-12)

	def test_skyframe_fk5_equinox(self):
		frame = ASTSkyFrame(system="FK5", equinox="1975.0")
		restored = self.roundtrip(frame)
		assert restored.system == "FK5"
		assert restored.equinox == "1975.0"
		assert restored.naxes == 2
		assert restored.astObject.get("Epoch") == frame.astObject.get("Epoch")

	def test_moc(self):
		circle = ASTCircle(center=[30, 45], radius=2.0)
		moc = circle.toMoc(max_order=8)
		restored = self.roundtrip(moc)
		assert restored.astObject.overlap(moc.astObject) == 5
		assert restored.maxOrder == moc.maxOrder

	def test_compound_region(self):
		c1 = ASTCircle(center=[30, 45], radius=2.0)
		c2 = ASTCircle(center=[32, 44], radius=1.0)
		compound = ASTCompoundRegion(regions=[c1, c2], operation=Ast.OR)
		restored = self.roundtrip(compound)
		assert restored.astObject.overlap(compound.astObject) == 5

	def test_wcs_attribute_rides_state(self):
		polygon = ASTPolygon.fromFITSHeader(GOLDEN_TAN_HEADER)
		restored = self.roundtrip(polygon)
		assert isinstance(restored.wcs, ASTFrameSet)
		np.testing.assert_allclose(restored.wcs.pix2world([100.0, 200.0]),
		                           polygon.wcs.pix2world([100.0, 200.0]), atol=1e-12)

	def test_reconstruct_garbage_valueerror(self):
		with pytest.raises(ValueError):
			bridge.reconstruct(ASTCircle, "garbage")

	def test_copy_and_deepcopy(self):
		''' SPEC-04 §4: copy() is independent of the original '''
		import copy as copy_module
		circle = ASTCircle(center=[30, 45], radius=2.0)
		for cloned in (circle.copy(), copy_module.deepcopy(circle)):
			assert cloned.isIdenticalTo(circle)
			cloned.negate()
			assert cloned.isNegated and not circle.isNegated


# ---------------------------------------------------------------------------
# Gate — no coordinate conversions outside the bridge
# ---------------------------------------------------------------------------

# angle-of-result conversions: reading out an *angle* value, not converting
# coordinates (SPEC-04A §5 tail / M22)
CONVERSION_ALLOWLIST = {
	("mapping/frame/frame.py", "return np.rad2deg(angle_rad) * u.deg"),
	("region/polygon.py", "angles.append(angle.to(u.rad).value)"),
}

CONVERSION_TOKENS = re.compile(r"deg2rad|rad2deg|\.to\(u\.rad\)|math\.radians|math\.degrees")


def test_no_conversions_outside_bridge():
	'''
	The SPEC-04 exit criterion (HANDOFF build table row 1), machine-enforced:
	no coordinate unit conversions exist outside _pyast_bridge.py, save the
	explicitly allowlisted angle-of-result lines.
	'''
	package_dir = os.path.join(os.path.dirname(__file__), "..", "cornish")
	violations = []
	for path in glob.glob(os.path.join(package_dir, "**", "*.py"), recursive=True):
		rel_path = os.path.relpath(path, package_dir)
		if rel_path == "_pyast_bridge.py":
			continue
		with open(path) as source_file:
			for line_number, line in enumerate(source_file, start=1):
				if CONVERSION_TOKENS.search(line):
					if (rel_path, line.strip()) in CONVERSION_ALLOWLIST:
						continue
					violations.append(f"{rel_path}:{line_number}: {line.strip()}")
	assert not violations, ("coordinate conversions found outside the bridge "
	                        "(add to _pyast_bridge or, for angle-of-result lines only, "
	                        "the allowlist):\n" + "\n".join(violations))
