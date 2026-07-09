'''
Tests for ASTMoc and STC-S serialization (reference implementation of SPEC-01/SPEC-02 subsets).
'''

import numpy as np
import pytest
import astropy.units as u
import starlink.Ast as Ast

from cornish import ASTCircle, ASTMoc, ASTPolygon, ASTRegion, ASTCompoundRegion


@pytest.fixture
def circle():
	return ASTCircle(center=[30.0, 45.0], radius=2.0 * u.deg)


def test_moc_from_region_area(circle):
	moc = ASTMoc.fromRegion(circle, max_order=10)
	analytic = 2 * np.pi * (1 - np.cos(np.deg2rad(2.0))) * (180 / np.pi) ** 2
	assert moc.area.to(u.deg * u.deg).value == pytest.approx(analytic, rel=2e-3)


def test_moc_containment(circle):
	moc = circle.toMoc(max_order=12)
	assert moc.pointInRegion([30.0, 45.0])
	assert not moc.pointInRegion([35.0, 45.0])


def test_moc_string_round_trip(circle):
	moc = circle.toMoc(max_order=8)
	serialized = moc.toString()
	restored = ASTMoc.fromString(serialized)
	assert restored.toString() == serialized
	assert restored.area.value == pytest.approx(moc.area.value)


def test_moc_json_form(circle):
	moc = circle.toMoc(max_order=8)
	json_string = moc.toString(json=True)
	assert json_string.startswith("{")


def test_moc_union_of_regions(circle):
	other = ASTCircle(center=[130.0, -20.0], radius=1.0 * u.deg)
	moc = ASTMoc.fromRegions([circle, other], max_order=10)
	assert moc.pointInRegion([30.0, 45.0])
	assert moc.pointInRegion([130.0, -20.0])
	area_circle = 2 * np.pi * (1 - np.cos(np.deg2rad(2.0))) * (180 / np.pi) ** 2
	area_other = 2 * np.pi * (1 - np.cos(np.deg2rad(1.0))) * (180 / np.pi) ** 2
	assert moc.area.to(u.deg * u.deg).value == pytest.approx(area_circle + area_other, rel=5e-3)


def test_moc_bounding_circle_of_compound_region(circle):
	'''The wish-list scenario: a finite bounding circle for a compound region, via MOC.'''
	other = ASTCircle(center=[32.0, 44.0], radius=1.0 * u.deg)
	compound = ASTCompoundRegion(regions=[circle, other], operation=Ast.OR)
	bounding = compound.toMoc(max_order=9).boundingCircle()
	assert np.isfinite(bounding.radius.value)
	assert bounding.radius < 5 * u.deg


def test_stcs_write_circle(circle):
	stcs = circle.toSTCS()
	assert stcs.startswith("Circle ICRS")
	assert "30" in stcs and "45" in stcs and "2" in stcs


def test_stcs_read_circle():
	region = ASTRegion.fromSTCS("Circle ICRS 30.0 45.0 2.0")
	assert isinstance(region, ASTCircle)
	assert region.radius.to(u.deg).value == pytest.approx(2.0)
	assert region.pointInRegion([30.0, 45.0])


def test_stcs_round_trip_polygon():
	polygon = ASTPolygon.fromFITSHeader({
		"NAXIS": 2, "NAXIS1": 512, "NAXIS2": 512,
		"CRPIX1": 256.5, "CRPIX2": 256.5, "CRVAL1": 30.0, "CRVAL2": 45.0,
		"CDELT1": -0.001, "CDELT2": 0.001,
		"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"})
	stcs = polygon.toSTCS()
	assert stcs.startswith("Polygon ICRS")
	restored = ASTRegion.fromSTCS(stcs)
	assert restored.isIdenticalTo(polygon)


def test_stcs_compound_region(circle):
	other = ASTCircle(center=[32.0, 44.0], radius=1.0 * u.deg)
	compound = ASTCompoundRegion(regions=[circle, other], operation=Ast.OR)
	stcs = compound.toSTCS()
	assert stcs.startswith("Union ICRS")
	restored = ASTRegion.fromSTCS(stcs)
	assert restored.pointInRegion([30.0, 45.0])
	assert restored.pointInRegion([32.0, 44.0])


def test_stcs_invalid_string_raises():
	with pytest.raises(ValueError):
		ASTRegion.fromSTCS("This is not a region")


def test_moc_to_polygon_not_implemented(circle):
	'''MOC boundary meshes are unordered; toPolygon must refuse rather than return bad geometry.'''
	with pytest.raises(NotImplementedError):
		circle.toMoc(max_order=8).toPolygon()


def test_moc_from_invalid_string_raises():
	with pytest.raises(ValueError):
		ASTMoc.fromString("this is not a MOC")


def test_moc_add_non_sky_region_raises():
	'''Adding a region in a non-sky frame to a MOC raises NotASkyRegion, not a raw AST error.'''
	from cornish.exc import NotASkyRegion
	from cornish.mapping import ASTFrame
	import starlink.Ast as Ast
	pixel_frame = Ast.Frame(2)
	pixel_circle = Ast.Circle(pixel_frame, 1, [50, 50], [10])
	with pytest.raises(NotASkyRegion):
		ASTMoc(max_order=8).add(pixel_circle)


def test_moc_from_string_none_raises_typeerror():
	with pytest.raises(TypeError):
		ASTMoc.fromString(None)


def test_stcs_from_none_raises_typeerror():
	with pytest.raises(TypeError):
		ASTRegion.fromSTCS(None)


def test_moc_max_order_must_be_integer():
	with pytest.raises(TypeError):
		ASTMoc(max_order=2.5)


def test_empty_moc_bounding_circle_raises():
	'''An empty MOC has no bounding circle; must raise rather than return NaN values.'''
	c1 = ASTCircle(center=[30.0, 45.0], radius=1.0 * u.deg)
	c2 = ASTCircle(center=[130.0, -40.0], radius=1.0 * u.deg)
	empty = ASTMoc(max_order=8)
	empty.add(c1)
	empty.add(c2, operation=Ast.AND)  # AND of disjoint regions -> empty coverage
	with pytest.raises(ValueError):
		empty.boundingCircle()


def test_box_stcs_documented_limitation():
	'''AST 9.3's STC-S reader cannot parse the PositionInterval form its writer emits for
	a Box; cornish surfaces this as ValueError (documented in toSTCS). If this test ever
	fails because the round-trip WORKS, AST fixed it upstream: update toSTCS's docstring.'''
	from cornish import ASTBox, ASTICRSFrame
	box = ASTBox.fromCorners(frame=ASTICRSFrame(), corners=([10, 10], [12, 12]))
	stcs = box.toSTCS()
	assert stcs.startswith("PositionInterval")
	with pytest.raises(ValueError):
		ASTRegion.fromSTCS(stcs)


def test_moc_max_order_rejects_bool():
	with pytest.raises(TypeError):
		ASTMoc(max_order=True)


def test_moc_uncertainty_property_accessible():
	'''Opus F1: ASTMoc skipped super().__init__, leaving inherited .uncertainty broken.'''
	assert ASTMoc(max_order=8).uncertainty is None
	assert ASTMoc().uncertainty is None


def test_moc_add_wrong_type_raises_typeerror():
	'''Opus F2: wrong type -> TypeError (policy: TypeError for types, ValueError for values).'''
	with pytest.raises(TypeError):
		ASTMoc(max_order=8).add(123)
	with pytest.raises(TypeError):
		ASTMoc(ast_object=123)


def test_empty_moc_max_order_is_none():
	assert ASTMoc().maxOrder is None
	assert ASTMoc(max_order=9).maxOrder == 9
