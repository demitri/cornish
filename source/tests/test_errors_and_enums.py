
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
