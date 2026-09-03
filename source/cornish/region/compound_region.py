
import logging
import warnings
from typing import Iterable, Union, List, Tuple

import astropy.units as u
import starlink.Ast as Ast

from .box import ASTBox
from .circle import ASTCircle
from .region import (ASTRegion, _unwrap_region, _require_compatible_frames,
                     _collect_wcs, _AST_COMBINE_FAILED_MSG)
from .polygon import ASTPolygon
from ..enums import RegionOperation
from ..exc import CoordinateSystemsCouldNotBeMapped
from .. import _pyast_bridge as bridge
from .._validation import as_integer
from .._registry import register_region, region_from_ast

__all__ = ["ASTCompoundRegion"]

logger = logging.getLogger("cornish")

def _symmetric_difference_ast(first: Ast.Region, second: Ast.Region) -> Ast.CmpRegion:
	'''
	``(first & ~second) | (~first & second)``, entirely at the pyast level.

	Operands must already be unwrapped and frame-gated by the caller. Raw pyast
	regions have no Python operators, so the composition is hand-built; the
	negations act on COPIES (``decompose()`` hands out live pointers, and the
	caller's operands must not be mutated), and every ``Ast.CmpRegion`` argument
	is positional because pyast silently ignores the ``oper=`` keyword.

	cornish builds symmetric differences this way rather than as a native
	``Ast.XOR`` node because AST serializes a native XOR to STC-S as
	``Difference``, which IVOA STC 1.33 defines as the ASYMMETRIC difference
	(SPEC-03A §2.7 / D25). AST's own internal storage of an XOR *is* this
	composition, so nothing is lost but the top-level ``Operator = 3`` label.
	'''
	not_first = first.copy();   not_first.negate()
	not_second = second.copy(); not_second.negate()
	left = Ast.CmpRegion(first, not_second, int(RegionOperation.AND))
	right = Ast.CmpRegion(not_first, second, int(RegionOperation.AND))
	return Ast.CmpRegion(left, right, int(RegionOperation.OR))


def _balanced_fold(operands, operation: RegionOperation) -> Ast.CmpRegion:
	'''
	Reduce an already-gated flat operand list to one ``Ast.CmpRegion``, pairwise,
	in ceil(log2 n) levels (SPEC-03A §3.2).

	The shape matters: every ``Ast.CmpRegion`` deep-copies its operands, so a
	left fold re-copies the growing tree (O(n^2) copy volume). Measured at
	n = 1000: 1.34 s to build left-folded versus 37 ms balanced, and
	``overlap``/``getregiondisc`` blow a two-minute budget on the left fold
	versus ~8 s balanced. ``pointinregion`` is shape-INsensitive (it tests every
	leaf regardless of nesting), which is why it is not the metric a regression
	test should watch.

	No ``.wcs``, uncertainty or gate logic lives here: the frame gate runs once
	per input against the first operand up front (compound intermediates inherit
	their children's frames, so re-gating per level is redundant), and the
	``.wcs`` rule is computed once by the caller over the flat input list.

	:param operands: cornish and/or raw regions; never XOR-composed here (the
		caller routes XOR through :func:`_symmetric_difference_ast`, since this
		fold would build native ``Ast.XOR`` nodes)
	'''
	def combine(first, second):
		# unwrapping defensively at every call: `regs` may hold cornish wrappers,
		# raw operands, and — at later levels — raw intermediates alongside a
		# still-wrapped odd survivor riding up.
		# kwarg-swallowing trap: `oper=` is silently ignored — positional only.
		return Ast.CmpRegion(bridge._unwrap(first), bridge._unwrap(second), int(operation))

	regs = list(operands)
	while len(regs) > 1:
		nxt = [combine(regs[i], regs[i + 1]) for i in range(0, len(regs) - 1, 2)]
		if len(regs) % 2:
			nxt.append(regs[-1]) # the odd leftover rides up a level
		regs = nxt
	# the n == 1 edge is unreachable from the frozen call sites (the constructor
	# requires >= 2 regions and union()/intersection() require >= 1 other); the
	# unwrap exists so this helper is total
	return bridge._unwrap(regs[0])


@register_region("CmpRegion")
class ASTCompoundRegion(ASTRegion):
	'''
	A region that is created as the composite of multiple regions.

	Regions are composited two at a time in the order they are supplied,
	e.g. if regions=[r1, r2, r3, r4]
	the result would be

	    region = compound( compound( compound(r1, r2), r3 ) r4 )

	all using the same operator as specified.

	:param regions: a list of regions to compound
	:param operation: one of :class:`cornish.enums.RegionOperation` (or the equivalent ``starlink.Ast.AND``/``OR``/``XOR`` constants)
	'''
	def __init__(self, ast_object=None, regions:Iterable[Union[ASTRegion, Ast.Region]]=None, operation:int=None):
		if ast_object is not None:
			# `is not None`, not truthiness: an empty regions list or a zero-valued
			# operation must still be rejected as a conflicting parameter
			if regions is not None or operation is not None:
				raise ValueError("If the ast_object is specified, no other parameter is accepted.")
			if not isinstance(ast_object, Ast.CmpRegion):
				raise TypeError(f"The 'ast_object' provided was not of type starlink.Ast.CmpRegion (got '{type(ast_object)}').")
			super().__init__(ast_object=ast_object)
			return

		if regions is None:
			raise ValueError("The 'regions' parameter must contain at least two regions.")
		regions = list(regions) # accept any iterable; the loop below consumes this copy
		if len(regions) < 2:
			raise ValueError("The 'regions' parameter must contain at least two regions.")
		for r in regions:
			_unwrap_region(r) # TypeError naming the type; a frame operand must not
			                  # reach the frame gate, which it would satisfy (§2.2)

		if operation is None:
			raise ValueError("An 'operation' must be specified (one of RegionOperation.AND/OR/XOR).")
		operation = RegionOperation(as_integer(operation, "operation")) # ValueError for out-of-range ints

		# The SPEC-03A §2.2 frame gate runs UNCONDITIONALLY, before the operation
		# split: a gate scoped to the non-XOR arm would leave XOR ungated, and
		# this constructor is the most ordinary construction path in the codebase.
		first = _unwrap_region(regions[0])
		for r in regions[1:]:
			_require_compatible_frames(first, _unwrap_region(r))

		if operation == RegionOperation.XOR and len(regions) != 2:
			raise ValueError(f"operation=XOR requires exactly 2 regions (got "
			                 f"{len(regions)}); n-ary symmetric differences are "
			                 f"not supported — compose explicitly")

		wcs = _collect_wcs(regions) # §2.3, ONCE over the input list

		try: # §2.2 boundary catch, spanning both arms
			if operation == RegionOperation.XOR:
				# never a native Ast.XOR node (§2.7 / D25) — the explicit composition
				super().__init__(ast_object=_symmetric_difference_ast(first, _unwrap_region(regions[1])))
			else:
				super().__init__(ast_object=_balanced_fold(regions, operation))
		except Ast.AstError as e:
			raise CoordinateSystemsCouldNotBeMapped(_AST_COMBINE_FAILED_MSG) from e

		self.wcs = wcs

	@property
	def operation(self) -> RegionOperation:
		'''
		How this compound's two children combine: AND, OR, or XOR.

		Read from the object's text dump: pyast exposes no attribute for it
		(``get("Oper")`` raises ``BADAT``), and the dump's label is ``Operator``,
		not the C attribute name ``Oper``.

		Note that XOR is not exclusive to foreign compounds: AST's own
		``CmpRegion`` constructor recognizes the ``(a & ~b) | (~a & b)``
		composition that cornish builds for a symmetric difference and stores it
		as a native XOR node, so ``(a ^ b).operation`` reads XOR too (verified
		2026-09-03 — SPEC-03A §2.7 predicted OR here and was corrected).
		'''
		return RegionOperation(int(bridge.dump_value(self.astObject, "Operator")))

	@property
	def uncertainty(self):
		'''
		Always ``None``: a compound's boundary uncertainty is COMPONENT-BORNE (D26).

		``Ast.CmpRegion`` accepts no uncertainty argument (a fourth positional is
		parsed as an attribute-settings string) and a compound's dump carries no
		top-level ``Unc``. What governs a compound's overlap semantics is the
		components' own uncertainties, which are copied in with the operands:
		verified — compounds of circles carrying 0.05 deg uncertainties, offset
		by 0.005 deg, report IDENTICAL, where the same geometry with AST-default
		uncertainties reports PARTIAL.

		So :class:`ASTCompoundRegion` takes no ``uncertainty`` parameter, and the
		three-way D15 default does not apply to compounds. Set the uncertainty
		you want on the operands before combining them.
		'''
		return None

	@uncertainty.setter
	def uncertainty(self, unc):
		raise NotImplementedError(
			"A compound region's uncertainty is component-borne (D26): set 'uncertainty' "
			"on the component regions before combining them."
		)

	@property
	def points(self):
		raise NotImplementedError(
			"Compound regions do not return points (Ast.CmpRegion objects do not carry them — "
			"AST returns an empty array, which would be a silent nothing); walk the tree with "
			"the 'components' property (or 'flattened()') and read .points from each leaf."
		)

	@property
	def area(self):
		'''
		Not implemented for compound regions.

		No approximation is returned behind a property: use
		``region.toMoc(max_order=...).area``, which is exact to the MOC's
		resolution and says what that resolution is.
		'''
		raise NotImplementedError(
			"exact areas of compound regions are not implemented; "
			"region.toMoc(max_order=...).area gives the area to the MOC's resolution"
		)

	@property
	def components(self) -> Tuple[ASTRegion, ASTRegion]:
		'''
		The two child regions (copies, registry-wrapped).

		AST composes regions pairwise, so a compound always has exactly two
		children; the children of an n-ary cornish union are themselves
		:class:`ASTCompoundRegion` objects (use :meth:`flattened` for the
		leaves).

		The children are COPIES: ``decompose()`` returns live pointers into this
		compound (verified), and handing those out would let
		``compound.components[0].negate()`` silently corrupt the parent.

		Foreign native-XOR compounds decompose to AST's internal composition —
		(¬A ∧ B) and (A ∧ ¬B) — not the two original operands (the originals are
		not recoverable through the sanctioned ``decompose`` API).
		cornish-built symmetric differences have the same structure, so the two
		cases introspect identically.

		This property does not cache: each access re-decomposes and re-copies,
		so a child is never stale after the parent is negated or re-attributed.
		'''
		a, b, *_ = self.astObject.decompose()
		return (region_from_ast(a.copy()), region_from_ast(b.copy()))

	def componentRegions(self) -> List[ASTRegion]:
		'''
		Deprecated alias for :attr:`components` (returns a list rather than a tuple).

		.. deprecated::
			Use :attr:`components`.
		'''
		warnings.warn("ASTCompoundRegion.componentRegions() is deprecated; use the "
		              "'components' property instead.", DeprecationWarning, stacklevel=2)
		return list(self.components)

	def flattened(self) -> Tuple[ASTRegion, ...]:
		'''
		All non-compound leaf regions of this tree, left-to-right, as copies.

		The walk is iterative (an explicit stack): a deep foreign left-fold
		compound must not hit Python's recursion limit. The compound-vs-leaf test
		is class-free — ``decompose()`` on a LEAF region returns ``None`` as its
		second element (verified across Circle/Box/Polygon/Moc) — which also
		keeps the region-class-ladder gate intact.

		Leaves carry their effective negation state as stored (e.g. the leaves of
		a difference or symmetric difference include negated copies). Flattening
		loses the operator structure by design — it answers "which pieces", not
		"which formula"; use :attr:`components` to walk structure.
		'''
		leaves = []
		stack = [self.astObject]
		while stack:
			node = stack.pop()
			a, b, *_ = node.decompose()
			if b is None:
				# a leaf: copy at emission (every decompose() return is live)
				leaves.append(region_from_ast(node.copy()))
			else:
				stack.append(b)
				stack.append(a) # a on top => left-to-right emission
		return tuple(leaves)

	def toPolygon(self, npoints=200, maxerr:u.Quantity=1.0*u.arcsec) -> ASTPolygon:
		'''
		Not yet implemented for compound regions (SPEC-03 repairs this class).

		The previous implementation was triply broken (missing import, a bound
		method passed as a frame, and — fundamentally — a "bounding box" built
		from the component regions' DEFINING points, e.g. a circle's centre and
		one circumference point, yielding a degenerate sliver). No output it
		ever produced was correct, so it raises loudly until SPEC-03 rebuilds
		it. For compound geometry today, convert to a MOC instead:
		``region.toMoc()`` has well-defined area/containment/bounding-circle support
		(to the resolution of its ``max_order``).
		'''
		raise NotImplementedError(
			"ASTCompoundRegion.toPolygon is not yet implemented (its previous algorithm was "
			"incorrect); use toMoc() for compound-geometry coverage, or take componentRegions() "
			"and combine their polygons yourself. SPEC-03 rebuilds this method."
		)
