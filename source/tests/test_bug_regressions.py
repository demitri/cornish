'''
Regression tests for bugs found and fixed in the 2026-07-09 review.
Each test names the bug it guards against.
'''

import io
import contextlib

import numpy as np
import pytest
import astropy.units as u
import starlink.Ast as Ast

from cornish import ASTCircle, ASTPolygon, ASTICRSFrame, ASTCompoundRegion, ASTFrameSet
from cornish.mapping import ASTMapping


def test_circle_from_ast_object():
	'''B1: ASTCircle(ast_object=...) raised unconditionally due to inverted validation.'''
	c1 = ASTCircle(center=[30, 45], radius=2.0)
	c2 = ASTCircle(ast_object=c1.astObject)
	assert c2.radius.to(u.deg).value == pytest.approx(2.0)


def test_compound_region_component_regions():
	'''B2: componentRegions() called the list it was returning.'''
	c1 = ASTCircle(center=[30, 45], radius=2.0)
	c2 = ASTCircle(center=[32, 44], radius=1.0)
	compound = ASTCompoundRegion(regions=[c1, c2], operation=Ast.AND)
	components = compound.componentRegions()
	assert len(components) == 2
	assert all(isinstance(c, ASTCircle) for c in components)


def test_compound_region_does_not_consume_callers_list():
	'''B3: the constructor pop()ed the regions out of the caller's list.'''
	c1 = ASTCircle(center=[30, 45], radius=2.0)
	c2 = ASTCircle(center=[32, 44], radius=1.0)
	regions = [c1, c2]
	ASTCompoundRegion(regions=regions, operation=Ast.OR)
	assert len(regions) == 2


def test_polygon_area_does_not_print():
	'''B5: ASTPolygon.area printed debug output for every vertex.'''
	polygon = ASTPolygon(frame=ASTICRSFrame(), points="((0,0), (90,0), (90,90))")
	buffer = io.StringIO()
	with contextlib.redirect_stdout(buffer):
		area = polygon.area
	assert buffer.getvalue() == ""
	assert area.to(u.sr).value == pytest.approx(np.pi / 2)


def test_frameset_from_raw_ast_frames():
	'''B7: fromFrames() referenced the nonexistent starlink.Ast.AstFrame.'''
	frame_set = ASTFrameSet.fromFrames(Ast.SkyFrame("System=ICRS"), Ast.SkyFrame("System=Galactic"))
	assert frame_set is not None


def test_frameset_invalid_ast_object_raises():
	'''frame_set.py built exceptions for invalid input without raising them.'''
	with pytest.raises(Exception):
		ASTFrameSet(ast_object="not a frame set")


def test_inverse_mapping():
	'''B10: inverseMapping() wrapped the None returned by pyast's in-place invert().'''
	zoom = ASTMapping(ast_object=Ast.ZoomMap(2, 2.0))
	inverse = zoom.inverseMapping()
	assert inverse.astObject is not None
	out = inverse.astObject.tran([[2.0], [4.0]], True)
	assert out[0][0] == pytest.approx(1.0)
	assert out[1][0] == pytest.approx(2.0)


def test_exception_name_spelling():
	'''B14: NoWCSFoumd typo; new name with compatibility alias.'''
	from cornish.exc import NoWCSFound, NoWCSFoumd
	assert NoWCSFoumd is NoWCSFound


def test_skyplot_constructs_on_current_ast():
	'''B16: SkyPlot passed an options string AST 9.x rejects; the entire plotting module failed.'''
	matplotlib = pytest.importorskip("matplotlib")
	matplotlib.use("Agg")
	from cornish.plot.matplotlib import SkyPlot
	circle = ASTCircle(center=[30, 45], radius=2.0 * u.deg)
	plot = SkyPlot(extent=circle, figsize=(4, 4))
	plot.addRegionOutline(circle)


def test_skyframe_from_fits_header():
	'''B12: ASTSkyFrame.fromFITSHeader referenced undefined names (ASTFrameSet, exc).'''
	from cornish import ASTSkyFrame
	header = {"NAXIS": 2, "NAXIS1": 100, "NAXIS2": 100,
	          "CRPIX1": 50.5, "CRPIX2": 50.5, "CRVAL1": 30.0, "CRVAL2": 45.0,
	          "CDELT1": -0.001, "CDELT2": 0.001,
	          "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"}
	frame = ASTSkyFrame.fromFITSHeader(header)
	assert frame.isSkyFrame


def test_skyframe_from_headerless_dict_raises_nowcsfound():
	'''fromFITSHeader promises NoWCSFound; a WCS-less header must not leak FrameNotFoundException.'''
	from cornish import ASTSkyFrame
	from cornish.exc import NoWCSFound
	with pytest.raises(NoWCSFound):
		ASTSkyFrame.fromFITSHeader({"NAXIS": 0})
