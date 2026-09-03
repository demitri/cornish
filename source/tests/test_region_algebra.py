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
