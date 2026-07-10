
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
	from cornish.mapping import ASTMapping
	return [ASTBox, ASTCircle, ASTPolygon, ASTCompoundRegion, ASTMoc,
	        ASTFrame, ASTSkyFrame, ASTFrameSet, ASTFITSChannel, ASTMapping]

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
	# a UnitMap is the wrong class for everything except ASTMapping itself
	# (a UnitMap IS an Ast.Mapping); use a FitsChan there
	from cornish.mapping import ASTMapping
	wrong_object = Ast.FitsChan() if wrapper_class is ASTMapping else Ast.UnitMap(2)
	with pytest.raises(TypeError):
		wrapper_class(ast_object=wrong_object)


def test_skyframe_empty_string_params_not_ignored():
	''' codex follow-up: falsy-but-provided parameters must conflict or fail, never be silently ignored '''
	from cornish.mapping.frame.sky_frame import ASTSkyFrame
	with pytest.raises(ValueError):
		ASTSkyFrame(ast_object=Ast.SkyFrame(), system="") # conflicting even when falsy
	with pytest.raises(ValueError):
		ASTSkyFrame(system="") # empty string is not a valid system, not an omitted one
	with pytest.raises(ValueError):
		ASTSkyFrame(equinox="") # chained from AST's parse error, not a raw Ast.DTERR
	with pytest.raises(ValueError):
		ASTSkyFrame(epoch="")


def test_attribute_setters_never_leak_ast_errors():
	'''
	Rule of two, mechanized: AST attribute setters (Equinox, then System) were
	found leaking raw Ast.AstError subclasses twice — all user-value setters
	now route through ASTObject._setAttribute, which chains into ValueError.
	'''
	frame = CornishFrame(naxes=2)
	with pytest.raises(ValueError) as excinfo:
		frame.system = "garbage"
	assert isinstance(excinfo.value.__cause__, Ast.AstError)
	from cornish.mapping.frame.sky_frame import ASTSkyFrame
	with pytest.raises(ValueError):
		ASTSkyFrame().equinox = "garbage"


def test_epoch_rejects_nonfinite():
	''' AST silently accepts Epoch=nan/inf; the cornish setter is the only gate '''
	from cornish.mapping.frame.sky_frame import ASTSkyFrame
	frame = ASTSkyFrame()
	for bad in ("nan", float("nan"), float("inf")):
		with pytest.raises(ValueError):
			frame.epoch = bad
	with pytest.raises(TypeError):
		frame.epoch = True
	# NumPy scalars are numbers (codex follow-up)
	frame.epoch = np.int64(1975)
	assert frame.epoch == pytest.approx(1975.0)
	frame.epoch = np.float32(1980.0)
	assert frame.epoch == pytest.approx(1980.0)
	# exotic Reals normalize through float() before reaching AST
	from fractions import Fraction
	frame.epoch = Fraction(3961, 2)
	assert frame.epoch == pytest.approx(1980.5)


def test_frame_both_ast_object_and_naxes_zero():
	''' truthiness fix pin: naxes=0 alongside ast_object is still a conflict '''
	with pytest.raises(ValueError):
		CornishFrame(ast_object=Ast.Frame(2), naxes=0)


def test_set_attribute_comma_injection_safe():
	'''
	sonnet round: astSet parses comma-separated settings lists, so an
	unquoted value containing ", Attr=Val" silently truncated the value AND
	set the other attribute. The seam quotes values; pin both halves.
	'''
	frame = CornishFrame(naxes=2)
	frame.domain = "AAA"
	frame.title = "Hello, Domain=BBB"
	assert frame.title == "Hello, Domain=BBB" # value preserved verbatim
	assert frame.domain == "AAA"              # other attribute untouched
	frame.setLabelForAxis(axis=1, label="RA, deg")
	assert frame.label(axis=1) == "RA, deg"


def test_mesh_size_rejects_too_small():
	''' sonnet round note: meshSize silently clamped values < 5 while claiming to require them '''
	circle = ASTCircle(center=[30, 45], radius=2.0)
	with pytest.raises(ValueError):
		circle.meshSize = 3
	with pytest.raises(TypeError):
		circle.meshSize = True
	circle.meshSize = 6
	assert circle.meshSize == 6
	# beyond MAX_MESH_SIZE, mesh spacing is finer than any boundary uncertainty
	# (and AST silently rewrites values beyond C int range besides)
	from cornish.region.region import MAX_MESH_SIZE
	with pytest.raises(ValueError):
		circle.meshSize = MAX_MESH_SIZE + 1
	with pytest.raises(ValueError):
		circle.meshSize = 2**31
	# both mesh methods route npoints through the same validated setter
	with pytest.raises(ValueError):
		circle.interiorPointMesh(npoints=3)
	with pytest.raises(ValueError):
		circle.boundaryPointMesh(npoints=3)
	with pytest.raises(TypeError):
		circle.interiorPointMesh(npoints=True)


def test_integer_parameters_accept_numpy_ints():
	'''
	Rule of two (third instance): integer parameters follow the moc.py
	operator.index pattern everywhere — NumPy ints accepted, floats and
	bools rejected, no silent truncation.
	'''
	circle = ASTCircle(center=[30, 45], radius=2.0)
	circle.meshSize = np.int64(50)
	assert circle.meshSize == 50
	mesh = circle.boundaryPointMesh(npoints=np.int64(40))
	assert len(mesh) > 0
	with pytest.raises(TypeError):
		circle.meshSize = 50.0
	# toSTCS digits
	assert circle.toSTCS(digits=np.int64(8)) == circle.toSTCS(digits=8)
	# frame axis accessors
	frame = CornishFrame(naxes=2)
	assert frame.unit(axis=np.int64(1)) == frame.unit(axis=1)
	# out-of-range axes are a loud ValueError, never a raw Ast.AXIIN leak
	for bad_axis in (0, -1, 3):
		with pytest.raises(ValueError):
			frame.unit(axis=bad_axis)
	# polygon downsize maxvert
	polygon = circle.toPolygon(npoints=50)
	import astropy.units as _u
	downsized = polygon.downsize(maxerr=1 * _u.arcsec, maxvert=np.int64(10))
	assert len(downsized.points) <= 10
	with pytest.raises(TypeError):
		polygon.downsize(maxerr=1 * _u.arcsec, maxvert=10.5)
