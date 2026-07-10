
'''
SPEC-10 tests: exception hierarchy, enums, the overlap truth table, and the
unraised-exception AST-walk gate.
'''

import ast as python_ast
import glob
import os
import warnings

import numpy as np
import pytest
import astropy.units as u
import starlink.Ast as Ast

import cornish.exc as exc
from cornish import ASTCircle, ASTICRSFrame
from cornish import ASTFrame as CornishFrame
from cornish.enums import (OverlapType, RegionOperation, MeshType, PixelComparison,
                           Marker, SkySystem, SpecSystem, TimeSystem)


# ---------------------------------------------------------------------------
# exception hierarchy
# ---------------------------------------------------------------------------

ALL_CUSTOM_EXCEPTIONS = [
	exc.FrameNotFoundException, exc.NotA2DRegion, exc.NotASkyRegion,
	exc.CoordinateSystemsCouldNotBeMapped, exc.NoWCSFound, exc.IncompleteHeader,
	exc.UnsupportedASTClass, exc.SerializationNotPossible, exc.RegionNotConnected,
]

@pytest.mark.parametrize("exception_class", ALL_CUSTOM_EXCEPTIONS)
def test_exception_parented(exception_class):
	''' every custom exception is importable and inherits from CornishError '''
	assert issubclass(exception_class, exc.CornishError)
	assert issubclass(exception_class, Exception)

def test_nowcsfoumd_deprecated_alias():
	assert exc.NoWCSFoumd is exc.NoWCSFound


# ---------------------------------------------------------------------------
# enums welded to the AST source of truth
# ---------------------------------------------------------------------------

def test_region_operation_welds():
	assert RegionOperation.AND == Ast.AND
	assert RegionOperation.OR == Ast.OR
	assert RegionOperation.XOR == Ast.XOR

def test_pixel_comparison_welds():
	assert PixelComparison.LT == Ast.LT
	assert PixelComparison.LE == Ast.LE
	assert PixelComparison.EQ == Ast.EQ
	assert PixelComparison.GE == Ast.GE
	assert PixelComparison.GT == Ast.GT
	assert PixelComparison.NE == Ast.NE

@pytest.mark.parametrize("member", list(SkySystem))
def test_sky_system_canonical(member):
	''' every SkySystem value is canonical: set() then get() returns it unchanged '''
	frame = Ast.SkyFrame(f"System={member.value}")
	assert frame.get("System") == member.value

@pytest.mark.parametrize("member", list(SpecSystem))
def test_spec_system_canonical(member):
	frame = Ast.SpecFrame(f"System={member.value}")
	assert frame.get("System") == member.value

@pytest.mark.parametrize("member", list(TimeSystem))
def test_time_system_canonical(member):
	frame = Ast.TimeFrame(f"System={member.value}")
	assert frame.get("System") == member.value

def test_mesh_type_values():
	''' AST's getregionmesh convention: 0 = interior ("surface"), 1 = boundary '''
	assert MeshType.SURFACE == 0
	assert MeshType.BOUNDARY == 1

def test_marker_values_distinct():
	''' the old dict silently lost the small-circle style to a duplicate key '''
	assert len({m.value for m in Marker}) == len(list(Marker))
	assert Marker.SMALL_CIRCLE == 1 and Marker.CIRCLE == 4


# ---------------------------------------------------------------------------
# overlap truth table (two circles arranged six ways + the unmappable case)
# ---------------------------------------------------------------------------

class TestOverlapTruthTable:

	@pytest.fixture
	def base(self):
		return ASTCircle(center=[30.0, 45.0], radius=2.0)

	def test_none(self, base):
		disjoint = ASTCircle(center=[100.0, -45.0], radius=2.0)
		assert base.overlapType(disjoint) == OverlapType.NONE
		assert base.overlaps(disjoint) is False

	def test_inside(self, base):
		bigger = ASTCircle(center=[30.0, 45.0], radius=5.0)
		assert base.overlapType(bigger) == OverlapType.INSIDE
		assert base.isFullyWithin(bigger) is True
		assert base.overlaps(bigger) is True

	def test_contains(self, base):
		smaller = ASTCircle(center=[30.0, 45.0], radius=0.5)
		assert base.overlapType(smaller) == OverlapType.CONTAINS
		assert base.fullyEncloses(smaller) is True

	def test_partial(self, base):
		shifted = ASTCircle(center=[32.0, 45.0], radius=2.0)
		assert base.overlapType(shifted) == OverlapType.PARTIAL
		assert base.overlaps(shifted) is True

	def test_identical(self, base):
		same = ASTCircle(center=[30.0, 45.0], radius=2.0)
		assert base.overlapType(same) == OverlapType.IDENTICAL
		assert base.isIdenticalTo(same) is True

	def test_negation(self, base):
		negated = ASTCircle(center=[30.0, 45.0], radius=2.0)
		negated.negate()
		assert base.overlapType(negated) == OverlapType.NEGATION
		assert base.isNegationOf(negated) is True
		assert base.overlaps(negated) is False

	def test_no_frame_mapping(self, base):
		# a domain-pinned pixel frame cannot be aligned with a sky frame
		# (a plain default-domain Frame CAN be, so the domain matters here)
		grid_frame = CornishFrame(naxes=2)
		grid_frame.domain = "GRID"
		pixel_circle = ASTCircle(frame=grid_frame, center=[50, 50], radius=10.0)
		assert base.overlapType(pixel_circle) == OverlapType.NO_FRAME_MAPPING
		with pytest.raises(exc.CoordinateSystemsCouldNotBeMapped):
			base.overlaps(pixel_circle)


# ---------------------------------------------------------------------------
# deprecated constants layer
# ---------------------------------------------------------------------------

def test_deprecated_constant_warns_and_maps():
	import cornish.constants as constants
	with pytest.warns(DeprecationWarning):
		assert constants.SYSTEM_ICRS == "ICRS"
	with pytest.warns(DeprecationWarning):
		# the old value "BEOPCH" was a typo AST rejects; the alias now yields the corrected value
		assert constants.SYSTEM_BESSELIAN == "BEPOCH"
	with pytest.warns(DeprecationWarning):
		assert constants.AST_BOUNDARY_MESH == 1

def test_unknown_constant_attributeerror():
	import cornish.constants as constants
	with pytest.raises(AttributeError):
		constants.NO_SUCH_CONSTANT

def test_package_import_clean_of_deprecations():
	''' importing cornish itself must not touch the deprecated alias layer '''
	import importlib, sys
	with warnings.catch_warnings():
		warnings.simplefilter("error", DeprecationWarning)
		saved = {name: module for name, module in sys.modules.items() if name.startswith("cornish")}
		for name in saved:
			del sys.modules[name]
		try:
			importlib.import_module("cornish")
		finally:
			sys.modules.update(saved)


# ---------------------------------------------------------------------------
# formerly-unraised exception paths + the AST-walk gate
# ---------------------------------------------------------------------------

def test_unraised_exception_gate():
	'''
	The rule-of-two gate (SPEC-10 §1): walk every module for expression
	statements whose value is a call to a name ending in Exception/Error —
	the constructed-but-never-raised pattern (three real instances found so far).

	Blind spot: this cannot catch an exception *assigned* to a variable and
	never raised; that residue remains a human-review duty.
	'''
	package_dir = os.path.join(os.path.dirname(__file__), "..", "cornish")
	offenders = []
	for path in glob.glob(os.path.join(package_dir, "**", "*.py"), recursive=True):
		with open(path) as source_file:
			tree = python_ast.parse(source_file.read(), filename=path)
		for node in python_ast.walk(tree):
			if not (isinstance(node, python_ast.Expr) and isinstance(node.value, python_ast.Call)):
				continue
			func = node.value.func
			name = func.id if isinstance(func, python_ast.Name) else \
				(func.attr if isinstance(func, python_ast.Attribute) else None)
			if name and (name.endswith("Exception") or name.endswith("Error")):
				offenders.append(f"{os.path.relpath(path, package_dir)}:{node.lineno}: {name}(...) constructed but never raised")
	assert not offenders, "\n".join(offenders)


def test_no_bare_exception_raised():
	''' SPEC-10 exit criterion (HANDOFF row 1b): nothing raises bare Exception in cornish/ '''
	package_dir = os.path.join(os.path.dirname(__file__), "..", "cornish")
	offenders = []
	for path in glob.glob(os.path.join(package_dir, "**", "*.py"), recursive=True):
		with open(path) as source_file:
			tree = python_ast.parse(source_file.read(), filename=path)
		for node in python_ast.walk(tree):
			if isinstance(node, python_ast.Raise) and isinstance(node.exc, python_ast.Call):
				func = node.exc.func
				if isinstance(func, python_ast.Name) and func.id == "Exception":
					offenders.append(f"{os.path.relpath(path, package_dir)}:{node.lineno}")
	assert not offenders, "bare `raise Exception(...)` found:\n" + "\n".join(offenders)


def test_boundary_mesh_ast_error_chained():
	''' SPEC-10 §1: the one AST-error catch site chains instead of printing '''
	import inspect
	from cornish.region.region import ASTRegion
	source = inspect.getsource(ASTRegion.boundaryPointMesh)
	assert "print(" not in source
	assert "from e" in source


def test_constants_star_import():
	''' codex round-1: `from cornish.constants import *` must still bind the deprecated names (warning) '''
	namespace = {}
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", DeprecationWarning)
		exec("from cornish.constants import *", namespace)
	assert namespace["SYSTEM_ICRS"] == "ICRS"
	assert namespace["AST_BOUNDARY_MESH"] == 1
	assert namespace["CENTER_CORNER"] == 0


def test_mesh_size_restored():
	'''
	codex round-1: a temporary MeshSize (npoints) must never leak — the
	restore is in a `finally`, so this pin covers the success path and the
	structure covers the failure path (a genuine getregionmesh failure is not
	constructible from a valid region).
	'''
	from cornish import ASTCircle
	circle = ASTCircle(center=[30, 45], radius=2.0)
	original = circle.meshSize
	mesh = circle.boundaryPointMesh(npoints=57)
	assert circle.meshSize == original
	circle.interiorPointMesh(npoints=57)
	assert circle.meshSize == original


# ---------------------------------------------------------------------------
# pyast kwarg-swallowing (the N1/oper=/frame= bug class) — behavior + gate
# ---------------------------------------------------------------------------

class TestCompoundRegionOperationTruthTable:
	'''
	Review finding (P1): pyast silently ignores the `oper=` keyword on
	Ast.CmpRegion — every compound region cornish ever built was an OR.
	These tests pin the now-positional delivery for all three operations.
	'''

	@pytest.fixture
	def overlapping_circles(self):
		c1 = ASTCircle(center=[30.0, 45.0], radius=2.0)
		c2 = ASTCircle(center=[32.0, 45.0], radius=2.0)
		return c1, c2

	def check(self, compound, in_both, in_first_only):
		''' membership probes at a point inside both circles / inside only the first '''
		assert compound.pointInRegion([31.0, 45.0]) is in_both
		assert compound.pointInRegion([28.5, 45.0]) is in_first_only

	def test_and(self, overlapping_circles):
		from cornish import ASTCompoundRegion
		from cornish import _pyast_bridge as bridge
		compound = ASTCompoundRegion(regions=overlapping_circles, operation=RegionOperation.AND)
		self.check(compound, in_both=True, in_first_only=False)
		assert int(bridge.dump_value(compound.astObject, "Operator")) == Ast.AND

	def test_or(self, overlapping_circles):
		from cornish import ASTCompoundRegion
		from cornish import _pyast_bridge as bridge
		compound = ASTCompoundRegion(regions=overlapping_circles, operation=Ast.OR)
		self.check(compound, in_both=True, in_first_only=True)
		assert int(bridge.dump_value(compound.astObject, "Operator")) == Ast.OR

	def test_xor(self, overlapping_circles):
		from cornish import ASTCompoundRegion
		from cornish import _pyast_bridge as bridge
		compound = ASTCompoundRegion(regions=overlapping_circles, operation=RegionOperation.XOR)
		self.check(compound, in_both=False, in_first_only=True)
		assert int(bridge.dump_value(compound.astObject, "Operator")) == Ast.XOR

	def test_region_add_operator_is_genuinely_and(self, overlapping_circles):
		'''
		the public `region1 + region2` (base-class ASTRegion.__add__) requests
		AND — and now actually gets it. (Circles override + for dilation and
		refuse region operands per D5, so exercise it via polygons.)
		'''
		p1 = overlapping_circles[0].toPolygon(npoints=60)
		p2 = overlapping_circles[1].toPolygon(npoints=60)
		compound = p1 + p2
		self.check(compound, in_both=True, in_first_only=False)

	def test_operation_validation(self, overlapping_circles):
		from cornish import ASTCompoundRegion
		with pytest.raises(ValueError):
			ASTCompoundRegion(regions=overlapping_circles)          # operation required
		with pytest.raises(TypeError):
			ASTCompoundRegion(regions=overlapping_circles, operation="AND")
		with pytest.raises(ValueError):
			ASTCompoundRegion(regions=overlapping_circles, operation=99)
		with pytest.raises(TypeError):
			ASTCompoundRegion(regions=[overlapping_circles[0], "not a region"], operation=Ast.OR)


def test_frameset_base_frame_constructor():
	''' review finding sibling: Ast.FrameSet's `frame=` kwarg was swallowed too — this ctor path was dead '''
	from cornish import ASTFrameSet, ASTICRSFrame
	fs = ASTFrameSet(base_frame=ASTICRSFrame())
	assert fs.currentFrame.system == "ICRS"
	fs2 = ASTFrameSet(base_frame=Ast.SkyFrame("System=Galactic"))
	assert fs2.currentFrame.system == "GALACTIC"


def test_no_kwargs_into_pyast_constructors():
	'''
	The rule-of-two gate for the pyast kwarg-swallowing bug class: pyast
	accepts and silently IGNORES unknown keyword arguments on many C-level
	constructors (verified victims: region `unc=`, CmpRegion `oper=`,
	FrameSet `frame=`). Every `Ast.<Class>(...)` call in cornish must pass
	its arguments positionally; a (callee, kwarg) pair may be allowlisted
	here only with a comment citing empirical verification that the kwarg
	is honored.

	Blind spots: cannot see calls built dynamically (getattr) or kwargs
	splatted from dicts, and relies on the project-wide `import starlink.Ast
	as Ast` convention — if pyast is ever imported under another name or its
	classes imported directly, extend this gate to resolve imports. None of
	these shapes exist today.
	'''
	import ast as python_ast
	import glob, os
	allowlist = set()  # e.g. ("SomeClass", "some_kwarg") — none verified safe today
	package_dir = os.path.join(os.path.dirname(__file__), "..", "cornish")
	offenders = []
	for path in glob.glob(os.path.join(package_dir, "**", "*.py"), recursive=True):
		with open(path) as source_file:
			tree = python_ast.parse(source_file.read(), filename=path)
		for node in python_ast.walk(tree):
			if not isinstance(node, python_ast.Call) or not node.keywords:
				continue
			func = node.func
			# match Ast.<UpperName>(...) — pyast classes are capitalized
			if isinstance(func, python_ast.Attribute) and isinstance(func.value, python_ast.Name) \
					and func.value.id == "Ast" and func.attr[:1].isupper():
				for keyword in node.keywords:
					if (func.attr, keyword.arg) not in allowlist:
						offenders.append(f"{os.path.relpath(path, package_dir)}:{node.lineno}: "
						                 f"Ast.{func.attr}({keyword.arg}=...) — pyast may silently ignore this kwarg; pass positionally")
	assert not offenders, "\n".join(offenders)


# ---------------------------------------------------------------------------
# ast_object validation gate (rule of two: moc.py and compound_region.py both
# had gaps here — mechanize the invariant across every wrapper)
# ---------------------------------------------------------------------------

def _wrapper_classes():
	from cornish import (ASTBox, ASTCircle, ASTPolygon, ASTCompoundRegion, ASTMoc,
	                     ASTFrame, ASTFrameSet, ASTFITSChannel)
	from cornish.mapping.frame.sky_frame import ASTSkyFrame
	return [ASTBox, ASTCircle, ASTPolygon, ASTCompoundRegion, ASTMoc,
	        ASTFrame, ASTSkyFrame, ASTFrameSet, ASTFITSChannel]

@pytest.mark.parametrize("wrapper_class", _wrapper_classes(), ids=lambda c: c.__name__)
@pytest.mark.parametrize("bad_object", ["garbage", 42], ids=["str", "int"])
def test_wrong_type_ast_object_rejected(wrapper_class, bad_object):
	'''
	Every ASTObject wrapper must reject a wrong-type ast_object with TypeError
	at construction — not accept it silently and fail confusingly later
	(the compound_region gap found in review).
	'''
	with pytest.raises(TypeError):
		wrapper_class(ast_object=bad_object)

@pytest.mark.parametrize("wrapper_class", _wrapper_classes(), ids=lambda c: c.__name__)
def test_wrong_ast_class_ast_object_rejected(wrapper_class):
	''' a genuine AST object of the WRONG class must also be rejected with TypeError '''
	with pytest.raises(TypeError):
		wrapper_class(ast_object=Ast.UnitMap(2))
