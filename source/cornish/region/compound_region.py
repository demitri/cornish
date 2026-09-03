
import logging
import warnings
from typing import Iterable, Union, List, Tuple

import astropy.units as u
import starlink.Ast as Ast

from .box import ASTBox
from .circle import ASTCircle
from .region import ASTRegion
from .polygon import ASTPolygon
from ..enums import RegionOperation
from .._validation import as_integer
from .._registry import register_region, region_from_ast

__all__ = ["ASTCompoundRegion"]

logger = logging.getLogger("cornish")

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
			if not isinstance(r, (ASTRegion, Ast.Region)):
				raise TypeError("The regions provided must be of type ASTRegion or starlink.Ast.Region.")

		if operation is None:
			raise ValueError("An 'operation' must be specified (one of RegionOperation.AND/OR/XOR).")
		operation = RegionOperation(as_integer(operation, "operation")) # ValueError for out-of-range ints

		# Maintain a list of regions that make up the compound region.
		# How will this work with different types of chained operations?
		# Postpone: too much effort for an initial version...
		#self._regions = list()

		r1 = None
		r2 = None
		compound_region = None

		while len(regions) > 0:
			if r1 is None:
				r1 = regions.pop(0) # get first item
				r2 = regions.pop(0)
			else:
				r1 = compound_region
				r2 = regions.pop(0)

			if isinstance(r1, ASTRegion):
				r1 = r1.astObject
			if isinstance(r2, ASTRegion):
				r2 = r2.astObject

			# pyast SILENTLY IGNORES the `oper=` keyword on Ast.CmpRegion (the same
			# argument-parsing bug class as region constructors' `unc=`): a kwarg
			# form parses fine and always builds an OR region. The operation MUST
			# be passed positionally — verified 2026-07-10; do not "clean this up".
			compound_region = Ast.CmpRegion( r1, r2, operation ) # todo: 'series' parameter?
			if compound_region is None:
				# not observed in practice (Ast.CmpRegion raises on error), but if it
				# ever happens a half-built object must not escape
				raise ValueError(f"AST could not combine the provided regions (Ast.CmpRegion returned None for {r1.Class} and {r2.Class}).")

		super().__init__(ast_object=compound_region)

	@property
	def points(self):
		raise NotImplementedError(
			"Compound regions do not return points (Ast.CmpRegion objects do not carry them); "
			"call componentRegions() and read .points from each component instead."
		)

	@property
	def area(self):
		'''
		The area of the compound region on the sphere. [Not yet implemented.]
		'''
		raise NotImplementedError()

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
