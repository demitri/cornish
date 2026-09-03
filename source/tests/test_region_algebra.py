'''
SPEC-03A tests: the type registry (§1), the operator algebra (§2), n-ary
operations (§3), and the ASTCompoundRegion repair (§4).

Test-plan row numbers (T1…T13) from SPEC-03A §6 appear in the section headers.
'''

import importlib
import logging
import warnings

import numpy as np
import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord
import starlink.Ast as Ast

import cornish
from cornish import (ASTBox, ASTCircle, ASTPolygon, ASTRegion, ASTCompoundRegion,
                     ASTMoc, ASTFrameSet, ASTICRSFrame)
from cornish import ASTFrame as CornishFrame
from cornish.mapping.frame.sky_frame import ASTSkyFrame
from cornish.mapping.frame.compound_frame import ASTCompoundFrame
from cornish._registry import (region_from_ast, frame_from_ast, register_region,
                               register_frame, _AST_REGION_CLASSES, _AST_FRAME_CLASSES)
from cornish.enums import RegionOperation, OverlapType
from cornish.exc import (UnsupportedASTClass, CoordinateSystemsCouldNotBeMapped,
                         RegionNotConnected, SerializationNotPossible)


# ---------------------------------------------------------------------------
# fixtures — the SPEC-03A §6 probe geometry
# ---------------------------------------------------------------------------

@pytest.fixture
def c1():
	''' unit circle at (10, 20) ICRS '''
	return ASTCircle(center=[10, 20], radius=1.0)

@pytest.fixture
def c2():
	''' unit circle at (11, 20) ICRS — overlaps c1 '''
	return ASTCircle(center=[11, 20], radius=1.0)

@pytest.fixture
def far():
	''' unit circle at (40, 20) ICRS — disjoint from both '''
	return ASTCircle(center=[40, 20], radius=1.0)

#: probe points: in both, c1 only, c2 only, neither
IN_BOTH = [10.5, 20.0]
C1_ONLY = [9.5, 20.0]
C2_ONLY = [11.5, 20.0]
NEITHER = [40.0, 40.0]

@pytest.fixture
def basic_circle():
	''' a circle on a plain 2-axis (non-sky) frame '''
	return ASTCircle(frame=CornishFrame(naxes=2), center=[50, 50], radius=10.0)


class _MinimalGrf:
	''' the smallest object Ast.Plot accepts as a graphics back end (test-only) '''
	def Attr(self, *args): return 0.0
	def BBuf(self, *args): pass
	def EBuf(self, *args): pass
	def Cap(self, *args): return 0
	def Flush(self, *args): pass
	def Line(self, *args): pass
	def Mark(self, *args): pass
	def Qch(self, *args): return (1.0, 1.0)
	def Scales(self, *args): return (1.0, 1.0)
	def Text(self, *args): pass
	def TxExt(self, *args): return ([0, 0, 0, 0], [0, 0, 0, 0])


def _ellipse() -> Ast.Ellipse:
	sky = Ast.SkyFrame("System=ICRS")
	return Ast.Ellipse(sky, 1, [np.radians(10), np.radians(20)],
	                   [np.radians(0.5), np.radians(0.2)], [0.0, 0.0])


# ===========================================================================
# T6 — the registry (§1)
# ===========================================================================

@pytest.mark.parametrize("build,expected", [
	(lambda: ASTCircle(center=[10, 20], radius=1.0), ASTCircle),
	(lambda: ASTBox.fromCorners(frame=ASTICRSFrame(), corners=([10, 20], [11, 21])), ASTBox),
	(lambda: ASTCircle(center=[10, 20], radius=1.0).toPolygon(npoints=20), ASTPolygon),
	(lambda: ASTCircle(center=[10, 20], radius=1.0).toMoc(max_order=6), ASTMoc),
	(lambda: ASTCompoundRegion(regions=[ASTCircle(center=[10, 20], radius=1.0),
	                                    ASTCircle(center=[11, 20], radius=1.0)],
	                           operation=RegionOperation.OR), ASTCompoundRegion),
], ids=["Circle", "Box", "Polygon", "Moc", "CmpRegion"])
def test_region_from_ast_round_trips_every_registered_class(build, expected):
	''' T6: every registered region class round-trips through the registry '''
	region = build()
	wrapped = region_from_ast(region.astObject)
	assert type(wrapped) is expected
	assert wrapped.astObject is region.astObject # by reference, per the contract


@pytest.mark.parametrize("build", [
	_ellipse,
	lambda: Ast.Interval(Ast.Frame(2), [0.0, 0.0], [1.0, 1.0]),
	lambda: Ast.NullRegion(Ast.Frame(2)),
	lambda: Ast.PointList(Ast.Frame(2), [[0.0], [1.0]]),
], ids=["Ellipse", "Interval", "NullRegion", "PointList"])
def test_region_from_ast_unregistered_class_raises(build):
	''' T6: an unwrapped region class raises loudly and names SPEC-05 '''
	with pytest.raises(UnsupportedASTClass) as excinfo:
		region_from_ast(build())
	assert "SPEC-05" in str(excinfo.value)


def test_region_from_ast_rejects_cornish_wrapper(c1):
	''' T6: a cornish wrapper is a TypeError with the already-wrapped message '''
	with pytest.raises(TypeError) as excinfo:
		region_from_ast(c1)
	assert "already a cornish wrapper" in str(excinfo.value)
	assert ".astObject" in str(excinfo.value)


def test_region_from_ast_rejects_non_region():
	''' T6: a non-region is a TypeError naming the received type '''
	with pytest.raises(TypeError) as excinfo:
		region_from_ast(Ast.UnitMap(2))
	assert "UnitMap" in str(excinfo.value)
	with pytest.raises(TypeError):
		region_from_ast("not a region")


def test_frame_from_ast_rejects_region(c1):
	''' T6: a Region satisfies isinstance(x, Ast.Frame) by MRO — refuse it explicitly '''
	with pytest.raises(TypeError) as excinfo:
		frame_from_ast(c1.astObject)
	assert "getregionframe" in str(excinfo.value)


def test_frame_from_ast_rejects_non_frame_and_wrapper():
	''' T6: non-frames and cornish wrappers are TypeErrors '''
	with pytest.raises(TypeError):
		frame_from_ast(Ast.UnitMap(2))
	with pytest.raises(TypeError) as excinfo:
		frame_from_ast(ASTICRSFrame())
	assert "already a cornish wrapper" in str(excinfo.value)


def test_registry_decorator_contract():
	'''
	T6: the §1.1 promises — a DIFFERENT class under a taken name is a ValueError;
	re-registering the SAME class is a no-op (module re-imports are harmless).
	'''
	class _Dummy: pass

	register_region("Interval")(_Dummy)
	try:
		assert _AST_REGION_CLASSES["Interval"] is _Dummy
		# same class again: no-op, no raise
		register_region("Interval")(_Dummy)
		assert _AST_REGION_CLASSES["Interval"] is _Dummy

		class _Other: pass
		with pytest.raises(ValueError) as excinfo:
			register_region("Interval")(_Other)
		assert "already registered" in str(excinfo.value)
	finally:
		del _AST_REGION_CLASSES["Interval"]

	# the frame table enforces the same contract
	register_frame("TimeFrame")(_Dummy)
	try:
		class _Other2: pass
		with pytest.raises(ValueError):
			register_frame("TimeFrame")(_Other2)
	finally:
		del _AST_FRAME_CLASSES["TimeFrame"]

	# only classes may be registered
	with pytest.raises(TypeError):
		register_region("Interval")(lambda: None)
	with pytest.raises(TypeError):
		register_region("")


def test_registry_tables_exact_contents():
	'''
	T6: the §1.2 tables, asserted exactly — a module that silently stops
	registering (or registers something new without updating SPEC-03A §1.2)
	fails here.
	'''
	assert _AST_REGION_CLASSES == {
		"Box": ASTBox,
		"Circle": ASTCircle,
		"CmpRegion": ASTCompoundRegion,
		"Moc": ASTMoc,
		"Polygon": ASTPolygon,
	}
	assert _AST_FRAME_CLASSES == {
		"SkyFrame": ASTSkyFrame,
		"CmpFrame": ASTCompoundFrame,
		"FrameSet": ASTFrameSet,
		"Plot": ASTFrameSet,   # a Plot IS a FrameSet (Δ2)
		"Frame": CornishFrame, # also the unregistered-class fallback
	}


@pytest.mark.parametrize("build,expected", [
	(lambda: Ast.SkyFrame("System=ICRS"), ASTSkyFrame),
	(lambda: Ast.CmpFrame(Ast.Frame(1), Ast.Frame(1)), ASTCompoundFrame),
	(lambda: Ast.Frame(2), CornishFrame),
	(lambda: Ast.FrameSet(Ast.Frame(2)), ASTFrameSet),                 # Δ2
	(lambda: Ast.Plot(Ast.Frame(2), [0, 0, 1, 1], [0, 0, 1, 1], _MinimalGrf()), ASTFrameSet),  # Δ2
	(lambda: Ast.FluxFrame(), CornishFrame),                           # Δ1
	(lambda: Ast.SpecFluxFrame(Ast.SpecFrame(), Ast.FluxFrame()), CornishFrame), # Δ1
	(lambda: Ast.SpecFrame(), CornishFrame),                           # Δ7
	(lambda: Ast.TimeFrame(), CornishFrame),                           # Δ7
	(lambda: Ast.DSBSpecFrame(), CornishFrame),                        # Δ7 / isinstance-order pitfall
], ids=["SkyFrame", "CmpFrame", "Frame", "FrameSet", "Plot", "FluxFrame",
        "SpecFluxFrame", "SpecFrame", "TimeFrame", "DSBSpecFrame"])
def test_frame_from_ast_dispatch(build, expected):
	''' T6: frame dispatch incl. the Δ1/Δ2/Δ7 fallbacks '''
	frame = frame_from_ast(build())
	assert type(frame) is expected


def test_frame_fallback_is_logged_not_silent(caplog):
	''' T6: the generic-ASTFrame fallback records the class name at DEBUG '''
	with caplog.at_level(logging.DEBUG, logger="cornish"):
		frame_from_ast(Ast.FluxFrame())
	assert any("FluxFrame" in record.getMessage() for record in caplog.records)


def test_delta7_master_conduct_pins():
	'''
	T6: Δ7's premises, pinned so SPEC-06's repair of these two wrappers is
	NOTICED (and this table updated) rather than silently diverging.
	'''
	with pytest.raises(ModuleNotFoundError):
		importlib.import_module("cornish.mapping.frame.spec_frame")
	from cornish.mapping.frame.time_frame import ASTTimeFrame
	with pytest.raises(NotImplementedError):
		ASTTimeFrame(ast_object=Ast.TimeFrame())


def test_fromAstObject_subclass_validation(c1):
	''' T6/§1.3: fromAstObject on a subclass validates the registry's answer '''
	assert type(ASTRegion.fromAstObject(c1.astObject)) is ASTCircle
	assert type(ASTCircle.fromAstObject(c1.astObject)) is ASTCircle
	with pytest.raises(TypeError) as excinfo:
		ASTPolygon.fromAstObject(c1.astObject)
	assert "ASTPolygon" in str(excinfo.value) and "ASTCircle" in str(excinfo.value)


# --- L1..L6 migration regressions ------------------------------------------

@pytest.mark.parametrize("stcs,expected", [
	("Circle ICRS 30.0 45.0 2.0", ASTCircle),
	("Polygon ICRS 10 20 10 21 11 21 11 20", ASTPolygon),
	("Union ICRS ( Circle 30 45 2 Circle 32 44 1 )", ASTCompoundRegion),
], ids=["circle", "polygon", "union"])
def test_L1_fromSTCS_returns_registry_classes(stcs, expected):
	''' L1: fromSTCS dispatches through the registry '''
	assert type(ASTRegion.fromSTCS(stcs)) is expected


def test_L2_regionWithMapping_class_and_wcs(c1):
	''' L2: regionWithMapping wraps via the registry and keeps propagating .wcs '''
	marker = object()
	c1.wcs = marker
	same_frame = c1.frame()
	mapped = c1.regionWithMapping(map=Ast.UnitMap(2), frame=same_frame)
	assert type(mapped) is ASTCircle
	assert mapped.wcs is marker

	# Δ3's sibling: a class with no registered wrapper raises loudly rather than
	# being silently genericized (mapping an ICRS circle into galactic yields an
	# Ast.Ellipse, which SPEC-05 has not wrapped yet)
	galactic = ASTSkyFrame(system="GALACTIC")
	with pytest.raises(UnsupportedASTClass):
		c1.regionWithMapping(map=c1.frame().astObject.convert(galactic.astObject),
		                     frame=galactic)


def test_L3_componentRegions_deprecated_delegate(c1, c2):
	''' L3: componentRegions() still works, now with a DeprecationWarning '''
	compound = ASTCompoundRegion(regions=[c1, c2], operation=RegionOperation.OR)
	with pytest.warns(DeprecationWarning):
		components = compound.componentRegions()
	assert isinstance(components, list) and len(components) == 2
	assert all(isinstance(r, ASTRegion) for r in components)


@pytest.mark.parametrize("build,expected", [
	(lambda: Ast.SkyFrame("System=ICRS"), ASTSkyFrame),
	(lambda: Ast.CmpFrame(Ast.Frame(1), Ast.Frame(1)), ASTCompoundFrame),
	(lambda: Ast.SpecFrame(), CornishFrame),
	(lambda: Ast.TimeFrame(), CornishFrame),
], ids=["SkyFrame", "CmpFrame", "SpecFrame", "TimeFrame"])
def test_L4_frameFromAstObject_delegates_to_registry(build, expected):
	''' L4: the public staticmethod is kept and now returns registry answers '''
	assert type(CornishFrame.frameFromAstObject(build())) is expected


def test_L4_frameFromAstObject_still_validates():
	''' L4: the None / wrong-type contracts are unchanged '''
	with pytest.raises(ValueError):
		CornishFrame.frameFromAstObject(None)
	with pytest.raises(TypeError):
		CornishFrame.frameFromAstObject(Ast.UnitMap(2))


def test_L6_frameAtIndex_registry_consistent():
	''' L6: frameAtIndex returns the same wrapper the registry would '''
	frame_set = ASTFrameSet.fromFITSHeader(fits_header={
		"NAXIS": 2, "NAXIS1": 512, "NAXIS2": 512,
		"CRPIX1": 256.5, "CRPIX2": 256.5, "CRVAL1": 30.0, "CRVAL2": 45.0,
		"CDELT1": -0.001, "CDELT2": 0.001,
		"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"})
	current = frame_set.frameAtIndex(Ast.CURRENT)
	assert type(current) is ASTSkyFrame
	base = frame_set.frameAtIndex(Ast.BASE)
	assert type(base) is type(frame_from_ast(frame_set.astObject.getframe(Ast.BASE)))


# ===========================================================================
# T1 — operator truth tables and operand handling (§2.1)
# ===========================================================================

def _membership(region):
	return tuple(bool(region.pointInRegion(p)) for p in (IN_BOTH, C1_ONLY, C2_ONLY, NEITHER))


@pytest.mark.parametrize("operator,expected", [
	("|", (True, True, True, False)),
	("&", (True, False, False, False)),
	("^", (False, True, True, False)),
	("-", (False, True, False, False)),
], ids=["union", "intersection", "symmetric_difference", "difference"])
def test_T1_operator_truth_table(c1, c2, operator, expected):
	''' T1: the four operators at the four probe points '''
	result = {"|": c1 | c2, "&": c1 & c2, "^": c1 ^ c2, "-": c1 - c2}[operator]
	assert isinstance(result, ASTCompoundRegion)
	assert _membership(result) == expected


@pytest.mark.parametrize("operator,method", [
	("|", "union"), ("&", "intersection"),
	("^", "symmetricDifference"), ("-", "difference"),
])
def test_T1_named_methods_match_operators(c1, c2, operator, method):
	''' T1: each named method is identical to its operator form '''
	operator_form = {"|": c1 | c2, "&": c1 & c2, "^": c1 ^ c2, "-": c1 - c2}[operator]
	named_form = getattr(c1, method)(c2)
	assert operator_form.overlapType(named_form) == OverlapType.IDENTICAL


@pytest.mark.parametrize("operator", ["|", "&", "^", "-"])
def test_T1_reflected_forms_preserve_written_order(c1, c2, operator):
	'''
	T1: raw pyast regions define no bitwise dunders, so a raw LEFT operand always
	reaches the reflected method — which must build what the written order means.
	'''
	raw = c1.astObject
	reflected = {"|": raw | c2, "&": raw & c2, "^": raw ^ c2, "-": raw - c2}[operator]
	direct = {"|": c1 | c2, "&": c1 & c2, "^": c1 ^ c2, "-": c1 - c2}[operator]
	assert reflected.overlapType(direct) == OverlapType.IDENTICAL


def test_T1_raw_minuend_is_not_the_subtrahend(c1, c2):
	''' T1: `raw - cornish` must be raw MINUS cornish, not the reverse '''
	raw_minus = c1.astObject - c2       # __rsub__: c1 is the minuend
	other_way = c2 - c1
	assert _membership(raw_minus) == (False, True, False, False)
	assert _membership(other_way) == (False, False, True, False)


def test_T1_pyast_regions_have_no_operators(c1, c2):
	''' T1 premise: two RAW regions cannot be combined by Python operators at all '''
	with pytest.raises(TypeError):
		c1.astObject | c2.astObject


@pytest.mark.parametrize("bad", [5, "x", None, 3.5], ids=["int", "str", "None", "float"])
def test_T1_wrong_kind_operand_operators(c1, bad):
	''' T1: non-region operands go through the NotImplemented protocol '''
	for operator in ("|", "&", "^", "-"):
		with pytest.raises(TypeError):
			eval(f"c1 {operator} bad", {"c1": c1, "bad": bad})


@pytest.mark.parametrize("bad", [5, None, "x"], ids=["int", "None", "str"])
def test_T1_wrong_kind_operand_named_methods(c1, bad):
	''' T1: named methods raise TypeError directly, naming the type '''
	for method in ("union", "intersection", "difference", "symmetricDifference"):
		with pytest.raises(TypeError) as excinfo:
			getattr(c1, method)(bad)
		assert type(bad).__name__ in str(excinfo.value)


@pytest.mark.parametrize("frame_operand", [
	lambda: Ast.Frame(2),
	lambda: Ast.SkyFrame("System=ICRS"),
	lambda: CornishFrame(naxes=2),
	lambda: ASTICRSFrame(),
], ids=["Ast.Frame", "Ast.SkyFrame", "ASTFrame", "ASTICRSFrame"])
def test_T1_frame_operands_are_TypeErrors_on_every_path(c1, frame_operand):
	'''
	T1: a frame slips past the permissive bridge unwrap AND satisfies the
	is_sky/naxes frame gate, so only the region-only coercion stops it. It must
	be a TypeError on EVERY construction path — never
	CoordinateSystemsCouldNotBeMapped.
	'''
	frame = frame_operand()
	with pytest.raises(TypeError):
		c1 | frame
	with pytest.raises(TypeError):
		c1.union(frame)
	with pytest.raises(TypeError):
		ASTCompoundRegion(regions=[c1, frame], operation=RegionOperation.OR)


def test_T1_plus_raises_for_every_operand(c1, c2):
	''' T1/§2.5: `+` never returns NotImplemented on the base class '''
	from cornish.region.region import _PLUS_REMOVED_MSG
	polygon = c1.toPolygon(npoints=20)
	for operand in (c2, c2.astObject, 5, "x", None):
		with pytest.raises(TypeError) as excinfo:
			polygon + operand
		assert str(excinfo.value) == _PLUS_REMOVED_MSG
	with pytest.raises(TypeError):
		5 + polygon
	with pytest.raises(TypeError):
		"x" + polygon


# ===========================================================================
# T2 — identities and compound uncertainty (§2.7 / D26)
# ===========================================================================

def test_T2_de_morgan(c1, c2):
	''' T2: ~(a | b) == (~a & ~b) and ~(a & b) == (~a | ~b) '''
	assert (~(c1 | c2)).overlapType((~c1) & (~c2)) == OverlapType.IDENTICAL
	assert (~(c1 & c2)).overlapType((~c1) | (~c2)) == OverlapType.IDENTICAL


def test_T2_symmetric_difference_welds_to_native_xor(c1, c2):
	''' T2: `a ^ b` is the same region as a raw native-XOR compound '''
	native = Ast.CmpRegion(c1.astObject, c2.astObject, int(RegionOperation.XOR))
	assert (c1 ^ c2).overlapType(native) == OverlapType.IDENTICAL


def test_T2_difference_equals_and_with_negation(c1, c2):
	''' T2: `a - b` == `a & ~b` '''
	assert (c1 - c2).overlapType(c1 & (~c2)) == OverlapType.IDENTICAL


def test_T2_compound_uncertainty_is_component_borne(c1, c2):
	''' T2/D26: a compound reports no uncertainty, and takes no such parameter '''
	for compound in (c1 | c2,
	                 c1.union(c2),
	                 ASTCompoundRegion(regions=[c1, c2], operation=RegionOperation.OR)):
		assert compound.uncertainty is None
	with pytest.raises(TypeError):
		ASTCompoundRegion(regions=[c1, c2], operation=RegionOperation.OR, uncertainty=1 * u.arcsec)


def test_T2_component_uncertainties_govern_compound_overlap():
	''' T2/X6: the components' uncertainties, not the compound's, decide overlap '''
	def pair(uncertainty):
		a = ASTCircle(center=[10, 20], radius=1.0, uncertainty=uncertainty)
		b = ASTCircle(center=[11, 20], radius=1.0, uncertainty=uncertainty)
		shifted_a = ASTCircle(center=[10.005, 20], radius=1.0, uncertainty=uncertainty)
		shifted_b = ASTCircle(center=[11.005, 20], radius=1.0, uncertainty=uncertainty)
		return (a | b), (shifted_a | shifted_b)

	coarse, coarse_shifted = pair(0.05 * u.deg)
	assert coarse.overlapType(coarse_shifted) == OverlapType.IDENTICAL
	default, default_shifted = pair(None)
	assert default.overlapType(default_shifted) == OverlapType.PARTIAL


# ===========================================================================
# T3 — mutation and liveness (§4.1 / D28)
# ===========================================================================

@pytest.mark.parametrize("operator", ["|", "&", "^", "-"])
def test_T3_operators_never_mutate_operands(c1, c2, operator):
	''' T3: both operands survive every operator untouched '''
	pristine_1, pristine_2 = c1.copy(), c2.copy()
	{"|": c1 | c2, "&": c1 & c2, "^": c1 ^ c2, "-": c1 - c2}[operator]
	assert c1.overlapType(pristine_1) == OverlapType.IDENTICAL
	assert c2.overlapType(pristine_2) == OverlapType.IDENTICAL
	assert c1.isNegated is False and c2.isNegated is False


def test_T3_invert_leaves_the_source_alone(c1):
	''' T3: ~a is a negated COPY '''
	negated = ~c1
	assert negated.isNegated is True
	assert c1.isNegated is False
	assert type(negated) is ASTCircle # not a compound


def test_T3_invert_carries_python_side_state(c1):
	'''
	T3: ~a preserves .wcs and _uncertainty through ASTObject.copy(). A future
	copy() that stops carrying Python-side attributes fails HERE.
	'''
	marker = object()
	c1.wcs = marker
	negated = ~c1
	assert negated.wcs is marker
	assert negated._uncertainty is c1._uncertainty


def test_T3_components_and_flattened_are_copies(c1, c2):
	''' T3/X9: mutating a returned child must never reach the parent '''
	compound = c1 | c2
	before = _membership(compound)
	compound.components[0].negate()
	compound.flattened()[0].negate()
	assert _membership(compound) == before

	# ...while raw decompose() children ARE live: a pyast change to copy
	# semantics would make this companion assertion fail, and it is the reason
	# 'components' copies at all
	raw_child = compound.astObject.decompose()[0]
	raw_child.negate()
	assert _membership(compound) != before


# ===========================================================================
# T4 — frame mixing (§2.1)
# ===========================================================================

def test_T4_galactic_and_icrs_combine_honestly(c1):
	''' T4: two sky systems map honestly; membership is right in both systems' terms '''
	galactic = ASTCircle(frame=ASTSkyFrame(system="GALACTIC"), center=[120, 30], radius=1.0)
	union = c1 | galactic
	assert union.pointInRegion([10, 20]) is True                       # c1's centre, ICRS
	assert union.pointInRegion(SkyCoord(l=120 * u.deg, b=30 * u.deg, frame="galactic")) is True


@pytest.mark.parametrize("reflected", [False, True], ids=["direct", "reflected"])
def test_T4_result_frame_is_the_first_written_operand(c1, reflected):
	''' T4: the result's System is the first-WRITTEN operand's, both orders '''
	galactic = ASTCircle(frame=ASTSkyFrame(system="GALACTIC"), center=[120, 30], radius=1.0)
	if reflected:
		icrs_first = c1.astObject | galactic
		galactic_first = galactic.astObject | c1
	else:
		icrs_first = c1 | galactic
		galactic_first = galactic | c1
	assert icrs_first.frame().system == "ICRS"
	assert galactic_first.frame().system == "GALACTIC"


# ===========================================================================
# T5 — the radian-trap frame gate (§2.2 / D24)
# ===========================================================================

def _all_construction_paths(first, second):
	''' every public way to build a compound from two operands '''
	return [
		("|", lambda: first | second),
		("&", lambda: first & second),
		("^", lambda: first ^ second),
		("-", lambda: first - second),
		("union", lambda: first.union(second)),
		("intersection", lambda: first.intersection(second)),
		("difference", lambda: first.difference(second)),
		("symmetricDifference", lambda: first.symmetricDifference(second)),
		("ctor/OR", lambda: ASTCompoundRegion(regions=[first, second], operation=RegionOperation.OR)),
		("ctor/AND", lambda: ASTCompoundRegion(regions=[first, second], operation=RegionOperation.AND)),
		("ctor/XOR", lambda: ASTCompoundRegion(regions=[first, second], operation=RegionOperation.XOR)),
	]


def test_T5_sky_x_basic_is_refused_on_every_path(c1, basic_circle):
	''' T5: the silent radian-map trap is gated on every construction path '''
	for name, build in _all_construction_paths(c1, basic_circle):
		with pytest.raises(CoordinateSystemsCouldNotBeMapped, match="RADIANS"):
			build()
	# ...and in the other operand order
	for name, build in _all_construction_paths(basic_circle, c1):
		with pytest.raises(CoordinateSystemsCouldNotBeMapped):
			build()


def test_T5_naxes_mismatch_is_refused_on_every_path(basic_circle):
	''' T5: differing axis counts are refused with the same taxonomy '''
	three_d = ASTCircle(frame=CornishFrame(naxes=3), center=[1, 1, 1], radius=1.0)
	for name, build in _all_construction_paths(basic_circle, three_d):
		with pytest.raises(CoordinateSystemsCouldNotBeMapped):
			build()


def test_T5_raw_ast_still_accepts_the_mixed_pair(c1, basic_circle):
	'''
	T5 behavior pin: raw AST still CONSTRUCTS the sky x basic compound and still
	reports overlap 1 (NONE) rather than 0 (NO_FRAME_MAPPING). An upstream fix
	must be noticed here rather than silently making the gate redundant.
	'''
	mixed = Ast.CmpRegion(c1.astObject, basic_circle.astObject, int(RegionOperation.OR))
	assert mixed is not None
	assert c1.astObject.overlap(basic_circle.astObject) == 1


def test_T5_raw_ast_does_catch_genuine_mismatches(c1, basic_circle):
	''' T5 behavior pin: the trap is exclusive to the same-naxes sky/basic mix '''
	three_d_sky_mismatch = Ast.Circle(Ast.Frame(3), 1, [1.0, 1.0, 1.0], [1.0])
	with pytest.raises(Ast.AstError):
		Ast.CmpRegion(c1.astObject, three_d_sky_mismatch, int(RegionOperation.OR))
	assert c1.astObject.overlap(three_d_sky_mismatch) == 0
	with pytest.raises(Ast.AstError):
		Ast.CmpRegion(basic_circle.astObject, three_d_sky_mismatch, int(RegionOperation.OR))
	assert basic_circle.astObject.overlap(three_d_sky_mismatch) == 0


def test_T5_overlapType_answers_the_mixed_pair_without_ast(c1, basic_circle):
	''' T5: overlapType is the honest non-raising API — it must not relay AST's NONE '''
	assert c1.overlapType(basic_circle) == OverlapType.NO_FRAME_MAPPING
	assert basic_circle.overlapType(c1) == OverlapType.NO_FRAME_MAPPING
	with pytest.raises(CoordinateSystemsCouldNotBeMapped):
		c1.overlaps(basic_circle)


def _domain_pair():
	pixel = ASTCircle(frame=CornishFrame(ast_object=Ast.Frame(2, "Domain=PIXEL")),
	                  center=[50, 50], radius=10.0)
	grid = ASTCircle(frame=CornishFrame(ast_object=Ast.Frame(2, "Domain=GRID")),
	                 center=[50, 50], radius=10.0)
	return pixel, grid


def test_T5_gate_passing_but_ast_refused_is_chained():
	'''
	T5/X30: equal-naxes basic frames with different Domains pass the cornish gate
	(Domain is AST's business) but AST refuses them. That escape must arrive as
	CoordinateSystemsCouldNotBeMapped with the Ast error chained — never a raw
	AstError.
	'''
	pixel, grid = _domain_pair()
	for name, build in _all_construction_paths(pixel, grid):
		with pytest.raises(CoordinateSystemsCouldNotBeMapped) as excinfo:
			build()
		assert isinstance(excinfo.value.__cause__, Ast.AstError), name
	# raw AST is honest about this pair
	assert pixel.astObject.overlap(grid.astObject) == 0


@pytest.mark.parametrize("method", ["union", "intersection"], ids=["union", "intersection"])
def test_T5_ast_refusal_escaping_at_a_later_fold_level(method):
	'''
	T5: three operands with the odd one LAST — the fold pairs the first two, then
	the mismatched survivor rides up and fails at the SECOND fold level, later
	than the gate loop that let it through. Still chained, never a raw AstError.
	'''
	pixel, grid = _domain_pair()
	more_pixel = ASTCircle(frame=CornishFrame(ast_object=Ast.Frame(2, "Domain=PIXEL")),
	                       center=[60, 50], radius=10.0)
	with pytest.raises(CoordinateSystemsCouldNotBeMapped) as excinfo:
		getattr(pixel, method)(more_pixel, grid)
	assert isinstance(excinfo.value.__cause__, Ast.AstError)


# ===========================================================================
# T8 — n-ary operations and tree shape (§3)
# ===========================================================================

def _tree_depth(ast_region) -> int:
	'''
	Maximum decompose() depth, LEAF-INCLUSIVE: the root counts as depth 1 and
	leaves are counted, so a balanced fold of n operands measures
	ceil(log2 n) + 1 for EVERY n, powers of two included. (The bare ceil(log2 n)
	is the edge count under a different convention; mixing the two is a known
	way to write a wrong assertion here.)
	'''
	deepest = 0
	stack = [(ast_region, 1)]
	while stack:
		node, depth = stack.pop()
		deepest = max(deepest, depth)
		a, b, *_ = node.decompose()
		if b is not None:
			stack.append((a, depth + 1))
			stack.append((b, depth + 1))
	return deepest


def _grid_of_boxes(n:int):
	''' n non-overlapping 0.4 deg boxes along a declination row '''
	return [ASTBox.fromCentreAndCorner(frame=ASTICRSFrame(),
	                                   centre=[10 + 1.5 * i, 20], corner=[10 + 1.5 * i + 0.2, 20.2])
	        for i in range(n)]


def test_T8_64_box_union_is_balanced():
	''' T8: a 64-operand union folds to depth ceil(log2 64) + 1 == 7 '''
	boxes = _grid_of_boxes(64)
	union = boxes[0].union(*boxes[1:])
	assert _tree_depth(union.astObject) == 7
	assert len(union.flattened()) == 64


def test_T8_membership_matches_any_of_the_operands():
	''' T8: the union's membership is exactly "any operand contains the point" '''
	boxes = _grid_of_boxes(16)
	union = boxes[0].union(*boxes[1:])
	for ra in np.linspace(9.0, 35.0, 60):
		for dec in (19.9, 20.0, 20.1):
			point = [float(ra), float(dec)]
			assert union.pointInRegion(point) == any(b.pointInRegion(point) for b in boxes)


def test_T8_union_matches_the_equivalent_moc():
	''' T8: the compound and a MOC built from the same pieces agree '''
	boxes = _grid_of_boxes(8)
	union = boxes[0].union(*boxes[1:])
	moc = ASTMoc.fromRegions(boxes, max_order=12)
	for ra in np.linspace(9.0, 22.0, 40):
		point = [float(ra), 20.0]
		assert union.pointInRegion(point) == moc.pointInRegion(point)


@pytest.mark.parametrize("n", [3, 5, 7], ids=["3", "5", "7"])
def test_T8_odd_operand_counts_fold_correctly(n):
	''' T8: the odd leftover rides up a level without being dropped or duplicated '''
	boxes = _grid_of_boxes(n)
	union = boxes[0].union(*boxes[1:])
	assert len(union.flattened()) == n
	for box in boxes:
		assert union.pointInRegion(list(np.asarray(box.points[0], dtype=float)))


@pytest.mark.parametrize("method", ["union", "intersection"])
def test_T8_zero_operands_raise_naming_the_method(c1, method):
	''' T8: each method's zero-operand message names ITSELF '''
	with pytest.raises(ValueError) as excinfo:
		getattr(c1, method)()
	assert f"{method}() needs at least one other region" in str(excinfo.value)


def test_T8_ctor_uses_the_balanced_fold(c1):
	''' T8/Δ4: the constructor's tree shape is balanced too, not a left fold '''
	boxes = _grid_of_boxes(8)
	compound = ASTCompoundRegion(regions=boxes, operation=RegionOperation.OR)
	assert _tree_depth(compound.astObject) == 4 # ceil(log2 8) + 1


def test_T8_ctor_xor_matches_the_operator_and_native_forms(c1, c2):
	''' T8/Δ5: the ctor's XOR arm composes exactly as `a ^ b` does '''
	from_ctor = ASTCompoundRegion(regions=[c1, c2], operation=RegionOperation.XOR)
	native = Ast.CmpRegion(c1.astObject, c2.astObject, int(RegionOperation.XOR))
	assert from_ctor.overlapType(c1 ^ c2) == OverlapType.IDENTICAL
	assert from_ctor.overlapType(native) == OverlapType.IDENTICAL
	# AST recognizes the (a & ~b) | (~a & b) pattern and stores it as a native
	# XOR node, so .operation reads XOR for cornish-built symmetric differences
	# too (verified 2026-09-03; SPEC-03A §2.7's "reads OR" claim was wrong and is
	# corrected in the spec)
	assert from_ctor.operation is RegionOperation.XOR
	assert (c1 ^ c2).operation is RegionOperation.XOR


def test_T8_ctor_xor_refuses_more_than_two_regions(c1, c2, far):
	''' T8/Δ5: n-ary XOR chains are refused rather than silently left-folded '''
	with pytest.raises(ValueError) as excinfo:
		ASTCompoundRegion(regions=[c1, c2, far], operation=RegionOperation.XOR)
	assert "exactly 2 regions" in str(excinfo.value)


def test_T8_wcs_unanimous_or_nothing(c1, c2):
	''' T8/§2.3: one origin carries; conflicting origins give None, never a guess '''
	W, W2 = object(), object()
	c1.wcs = W
	c2.wcs = None
	assert (c1 | c2).wcs is W
	assert ASTCompoundRegion(regions=[c1, c2], operation=RegionOperation.OR).wcs is W

	c2.wcs = W2
	assert (c1 | c2).wcs is None
	assert ASTCompoundRegion(regions=[c1, c2], operation=RegionOperation.OR).wcs is None

	c2.wcs = W # the same object, not merely equal
	assert (c1 | c2).wcs is W


def test_T8_wcs_is_evaluated_once_over_the_flat_list(far):
	'''
	T8/§2.3: the n >= 3 counterexample. Re-deriving the rule pairwise inside the
	fold gives (W, W2) -> None, then (None, W) -> W, silently resurrecting one of
	two conflicting origins.
	'''
	W, W2 = object(), object()
	a, b, c = _grid_of_boxes(3)
	a.wcs, b.wcs, c.wcs = W, W2, W
	assert a.union(b, c).wcs is None
	assert ASTCompoundRegion(regions=[a, b, c], operation=RegionOperation.OR).wcs is None

	a.wcs, b.wcs, c.wcs = None, W, None
	assert a.union(b, c).wcs is W
	assert ASTCompoundRegion(regions=[a, b, c], operation=RegionOperation.OR).wcs is W


def test_T8_wcs_mixed_raw_and_cornish_operands(c1, c2):
	''' T8/§2.3: raw operands contribute nothing — and never an AttributeError '''
	W = object()
	c1.wcs = W
	assert (c1 | c2.astObject).wcs is W
	assert (c2.astObject | c1).wcs is W
	assert ASTCompoundRegion(regions=[c2.astObject, c1], operation=RegionOperation.OR).wcs is W
	assert c1.union(c2.astObject).wcs is W


def test_T8_intermediate_compounds_carry_no_wcs():
	''' T8/§2.3: only the FINAL compound is given the collected .wcs '''
	W = object()
	a, b, c, d = _grid_of_boxes(4)
	for region in (a, b, c, d):
		region.wcs = W
	union = a.union(b, c, d)
	assert union.wcs is W
	assert all(child.wcs is None for child in union.components)


# ===========================================================================
# T11 — performance regression on the shape-sensitive operations (§3.3)
# ===========================================================================

def _best_of(callable_, repeats:int = 3) -> float:
	import time
	best = float("inf")
	for _ in range(repeats):
		start = time.perf_counter()
		callable_()
		best = min(best, time.perf_counter() - start)
	return best


def test_T11_balanced_fold_beats_a_left_fold():
	'''
	T11: the two SHAPE-SENSITIVE operations — construction and overlap — must be
	comfortably faster on the balanced tree. §3.3 measured 5x and 18x at n = 100;
	the 0.5 threshold is generous slack for a loaded CI machine. pointinregion is
	deliberately NOT the metric: it is shape-INsensitive (O(leaf count)), so it
	cannot distinguish the trees at all — asserted below as a pin of that finding.
	'''
	boxes = _grid_of_boxes(128)
	probe = ASTCircle(center=[10, 20], radius=0.05)

	def build_balanced():
		return boxes[0].union(*boxes[1:])

	def build_left():
		result = boxes[0].astObject
		for box in boxes[1:]:
			result = Ast.CmpRegion(result, box.astObject, int(RegionOperation.OR))
		return result

	balanced_build = _best_of(build_balanced)
	left_build = _best_of(build_left)
	assert balanced_build < 0.5 * left_build, (balanced_build, left_build)

	balanced, left = build_balanced().astObject, build_left()
	balanced_overlap = _best_of(lambda: balanced.overlap(probe.astObject))
	left_overlap = _best_of(lambda: left.overlap(probe.astObject))
	assert balanced_overlap < 0.5 * left_overlap, (balanced_overlap, left_overlap)

	# shape-insensitivity of pointinregion, pinned
	point = [10.0, 20.0]
	balanced_point = _best_of(lambda: balanced.pointinregion(np.radians(point)), repeats=5)
	left_point = _best_of(lambda: left.pointinregion(np.radians(point)), repeats=5)
	assert 1 / 3 < (balanced_point / left_point) < 3, (balanced_point, left_point)


# ===========================================================================
# T12 — dilation and scaling (§2.6 / D5)
# ===========================================================================

def test_T12_dilatedBy_sky(c1):
	''' T12: Quantity and bare-degree dilation, positive and negative '''
	assert c1.dilatedBy(0.5).radius.to_value(u.deg) == pytest.approx(1.5)
	assert c1.dilatedBy(30 * u.arcmin).radius.to_value(u.deg) == pytest.approx(1.5)
	assert c1.dilatedBy(-0.5).radius.to_value(u.deg) == pytest.approx(0.5)
	assert c1.radius.to_value(u.deg) == pytest.approx(1.0) # untouched


def test_T12_dilatedBy_basic(basic_circle):
	''' T12: bare numbers on a basic frame; a Quantity there is a ValueError '''
	assert basic_circle.dilatedBy(5).radius == pytest.approx(15.0)
	with pytest.raises(ValueError):
		basic_circle.dilatedBy(5 * u.deg)


def test_T12_dilatedBy_invalid(c1):
	''' T12: shrinking past zero is a ValueError naming both radii; a bool is a TypeError '''
	with pytest.raises(ValueError) as excinfo:
		c1.dilatedBy(-2.0)
	message = str(excinfo.value)
	assert "radius 1 deg" in message and "radius of -1 deg" in message # both radii named
	with pytest.raises(ValueError):
		c1.dilatedBy(-1.5 * u.deg)
	with pytest.raises(TypeError):
		c1.dilatedBy(True)
	with pytest.raises(TypeError):
		c1.dilatedBy(c1)


def test_T12_scaledBy(c1):
	''' T12: scaling contracts '''
	assert c1.scaledBy(2).radius.to_value(u.deg) == pytest.approx(2.0)
	with pytest.raises(TypeError):
		c1.scaledBy(True)
	with pytest.raises(TypeError):
		c1.scaledBy(2 * u.deg)
	for bad in (0, -1, float("inf"), float("nan")):
		with pytest.raises(ValueError):
			c1.scaledBy(bad)


def test_T12_operators_delegate(c1):
	''' T12: `+`, `*` and `/` are sugar for the named methods '''
	assert (c1 + 0.5).overlapType(c1.dilatedBy(0.5)) == OverlapType.IDENTICAL
	assert (c1 * 2).overlapType(c1.scaledBy(2)) == OverlapType.IDENTICAL
	assert (c1 / 2).overlapType(c1.scaledBy(0.5)) == OverlapType.IDENTICAL


def test_T12_plus_message_scope(c1, c2):
	'''
	T12: only the `+`/region branch adopts the removal text — a user who typed
	`circle * region` must not be lectured about `+`.
	'''
	from cornish.region.region import _PLUS_REMOVED_MSG
	with pytest.raises(TypeError) as plus:
		c1 + c2
	assert str(plus.value) == _PLUS_REMOVED_MSG
	with pytest.raises(TypeError) as times:
		c1 * c2
	assert "scale a circle" in str(times.value)
	assert "'+'" not in str(times.value)


def test_T12_no_radd_dilation(c1):
	''' T12: `5 + circle` raises — sugar is not owed commutativity (deliberate) '''
	from cornish.region.region import _PLUS_REMOVED_MSG
	with pytest.raises(TypeError) as excinfo:
		5 + c1
	assert str(excinfo.value) == _PLUS_REMOVED_MSG


@pytest.mark.parametrize("uncertainty,expect_unc_block", [
	("omitted", True),   # DEFAULT_UNCERTAINTY sentinel -> 1 arcsec materialized
	(None, False),       # explicit opt-out stays an opt-out
	(0.02 * u.deg, True),
], ids=["omitted", "explicit-None", "explicit-value"])
def test_T12_uncertainty_fidelity_all_three_states(uncertainty, expect_unc_block):
	''' T12/D15: dilation must not silently change which of the three states applies '''
	from cornish.region.region import DEFAULT_UNCERTAINTY
	if uncertainty == "omitted":
		circle = ASTCircle(center=[10, 20], radius=1.0)
		expected = DEFAULT_UNCERTAINTY
	else:
		circle = ASTCircle(center=[10, 20], radius=1.0, uncertainty=uncertainty)
		expected = uncertainty

	for grown in (circle.dilatedBy(0.5), circle.scaledBy(2.0)):
		# identity, not equality: the DEFAULT_UNCERTAINTY sentinel must pass
		# through as itself, or "omitted" silently becomes "explicitly 1 arcsec"
		assert grown._uncertainty is expected
		# and the AST object must actually carry (or not carry) an Unc block —
		# a faithful Python-side value over a wrong AST object is the SPEC-04
		# bug this pins against
		assert ("Unc =" in str(grown.astObject)) is expect_unc_block


def test_T12_wcs_carried_by_both_methods(c1):
	''' T12/§2.3: single-parent operations preserve provenance '''
	W = object()
	c1.wcs = W
	assert c1.dilatedBy(0.1).wcs is W
	assert c1.scaledBy(1.1).wcs is W
	assert (c1 + 0.1).wcs is W
	assert (c1 * 1.1).wcs is W
	assert (c1 / 1.1).wcs is W
