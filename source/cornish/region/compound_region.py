
import logging
import warnings
from typing import Iterable, Union, List, Tuple

import numpy as np
import astropy.units as u
import starlink.Ast as Ast

from .box import ASTBox
from .circle import ASTCircle
from .region import (ASTRegion, _unwrap_region, _require_compatible_frames,
                     _collect_wcs, _AST_COMBINE_FAILED_MSG, _degenerate_disc,
                     _disc_escape_points)
from .polygon import ASTPolygon
from ..enums import RegionOperation
from ..exc import CoordinateSystemsCouldNotBeMapped, RegionNotConnected
from .. import _pyast_bridge as bridge
from .._validation import as_integer
from .._registry import register_region, region_from_ast

__all__ = ["ASTCompoundRegion"]

logger = logging.getLogger("cornish")

#: A boundary gap is a SEAM (a jump from one component's arc to another's) when
#: it exceeds this multiple of the local median gap. Measured seams were 54x and
#: 927x the median in the probe geometry, while boundary curvature never
#: approaches 5x, so the separation is wide.
SEAM_FACTOR = 5

#: Width of the centred window for the LOCAL median gap. Local, not global:
#: MeshSize applies per component, so a union of a 0.1 deg and a 2 deg circle
#: has two legitimate spacing regimes and a global median would misclassify the
#: coarse regime's ordinary steps as seams.
_LOCAL_MEDIAN_WINDOW = 11


def _separations(points: np.ndarray, is_sky: bool) -> np.ndarray:
	'''
	Geodesic separation between consecutive points, in FRAME units.

	Angular separation (haversine) on a sky frame, Euclidean otherwise. Planar
	distance on sky coordinates is NOT acceptable: a polar component's mesh
	longitudes legitimately span the full 0..360 deg range — a dec 89.5 deg
	circle's mesh spans 0 to 356.4 deg with naive planar "steps" up to 356 deg
	where the true angular steps are a uniform 113 arcsec — so planar gaps would
	shred a polar arc into spurious seams.

	:param points: (n, 2) array in frame units (radians on a sky frame)
	'''
	first, second = points[:-1], points[1:]
	if not is_sky:
		return np.hypot(second[:, 0] - first[:, 0], second[:, 1] - first[:, 1])
	lon1, lat1 = first[:, 0], first[:, 1]
	lon2, lat2 = second[:, 0], second[:, 1]
	sin_half_lat = np.sin(0.5 * (lat2 - lat1))
	sin_half_lon = np.sin(0.5 * (lon2 - lon1))
	haversine = sin_half_lat ** 2 + np.cos(lat1) * np.cos(lat2) * sin_half_lon ** 2
	return 2.0 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))


def _local_medians(gaps: np.ndarray) -> np.ndarray:
	''' the median gap over a centred, end-clamped window, for every gap '''
	half = _LOCAL_MEDIAN_WINDOW // 2
	count = len(gaps)
	return np.array([np.median(gaps[max(0, i - half):min(count, i + half + 1)])
	                 for i in range(count)])


def _edge_median(gaps: np.ndarray, at_start: bool, fallback: float) -> float:
	'''
	The median of an arc's leading or trailing window of gaps.

	A one-point arc has no gaps and therefore no local scale of its own; a
	zero would make its chaining tolerance zero and strand it. ``fallback`` (the
	whole mesh's median spacing) stands in for exactly that case.
	'''
	if len(gaps) == 0:
		return fallback
	window = gaps[:_LOCAL_MEDIAN_WINDOW] if at_start else gaps[-_LOCAL_MEDIAN_WINDOW:]
	return float(np.median(window))


def _split_into_arcs(points: np.ndarray, is_sky: bool) -> list:
	'''
	Split a compound's boundary mesh into per-arc index ranges at its seams.

	Compound boundary meshes on AST 9.3 are ordered WITHIN arcs — consecutive
	steps match nearest-neighbour spacing — with a small number of large jumps
	where AST switches from one boundary arc to the next. That is what makes
	exact-vertex polygonization tractable without rasterization.
	'''
	gaps = _separations(points, is_sky)
	if len(gaps) == 0:
		return [points]
	medians = _local_medians(gaps)
	# a zero local median (a degenerate run of coincident points) must not make
	# every following gap a "seam"; compare against the largest sane scale seen
	positive = medians[medians > 0]
	floor = float(np.min(positive)) if len(positive) else 0.0
	thresholds = SEAM_FACTOR * np.where(medians > 0, medians, floor)
	break_after = np.nonzero(gaps > thresholds)[0]

	arcs, start = [], 0
	for index in break_after:
		arcs.append(points[start:index + 1])
		start = index + 1
	arcs.append(points[start:])
	return [arc for arc in arcs if len(arc) > 0]


def _chain_arcs(arcs: list, is_sky: bool) -> list:
	'''
	Chain boundary arcs into closed loops (§4.3 step 3).

	Greedy, with a deterministic lowest-index tie-break: append the unused arc
	whose nearest endpoint (either orientation — an arc may need reversing) lies
	within the chaining tolerance of the current loop end; when the loop end
	comes back within tolerance of the loop's own start, close it and begin a new
	loop with the next unused arc.

	The tolerance is local, like the seam test: ``SEAM_FACTOR`` times the larger
	of the two local median spacings at the two ends being joined. Arc endpoints
	at a union's waist sit within a step or two of each other, while arcs
	belonging to a different loop sit a seam away.

	An open chain is NEVER silently closed across a gap.

	:raises ValueError: no unused arc, and not the loop's own start, lies within
		tolerance of the current loop end
	'''
	def separation(point_a, point_b) -> float:
		return float(_separations(np.array([point_a, point_b]), is_sky)[0])

	arc_gaps = [_separations(arc, is_sky) for arc in arcs]
	all_gaps = np.concatenate([gaps for gaps in arc_gaps if len(gaps)]) if arc_gaps else np.array([])
	mesh_scale = float(np.median(all_gaps)) if len(all_gaps) else 0.0
	head_scale = [_edge_median(gaps, True, mesh_scale) for gaps in arc_gaps]
	tail_scale = [_edge_median(gaps, False, mesh_scale) for gaps in arc_gaps]

	unused = set(range(len(arcs)))
	loops = []
	while unused:
		index = min(unused)
		unused.discard(index)
		loop = [arcs[index]]
		loop_start_scale = head_scale[index]
		current_end = arcs[index][-1]
		current_scale = tail_scale[index]

		while True:
			# Closing the loop and continuing it are the SAME kind of choice, so
			# they compete on distance: whichever endpoint is nearest wins, ties
			# breaking on the lowest arc index. Testing closure first instead
			# closes the loop prematurely whenever the tolerance (5x the local
			# spacing) happens to exceed the gap back to the loop's start —
			# which it does at coarse npoints, mis-chaining an ordinary union.
			closing_tolerance = SEAM_FACTOR * max(current_scale, loop_start_scale)
			closing_distance = separation(current_end, loop[0][0])
			best = None
			best_distance = closing_distance if closing_distance <= closing_tolerance else float("inf")

			for candidate in sorted(unused):
				tolerance = SEAM_FACTOR * max(current_scale, head_scale[candidate],
				                              tail_scale[candidate])
				to_start = separation(current_end, arcs[candidate][0])
				to_end = separation(current_end, arcs[candidate][-1])
				nearest = min(to_start, to_end)
				if nearest <= tolerance and nearest < best_distance:
					best = (candidate, to_end < to_start)
					best_distance = nearest

			if best is None:
				if closing_distance <= closing_tolerance:
					break # the loop closes onto its own start
				# an open chain is never silently closed across a gap
				raise ValueError(
					"this compound region's boundary arcs could not be chained into closed loops: "
					"the sampling may be too coarse (increase npoints), or these components meet "
					"in a way this algorithm cannot resolve — use toMoc() for such geometry. An "
					"open chain is never silently closed across a gap.")
			candidate, reverse = best
			unused.discard(candidate)
			arc = arcs[candidate][::-1] if reverse else arcs[candidate]
			loop.append(arc)
			current_end = arc[-1]
			current_scale = head_scale[candidate] if reverse else tail_scale[candidate]

		loops.append(np.concatenate(loop))
	return loops


def _min_separation(point: np.ndarray, points: np.ndarray, is_sky: bool) -> float:
	''' the smallest geodesic separation between ``point`` and any of ``points`` '''
	if is_sky:
		lon, lat = points[:, 0], points[:, 1]
		sin_half_lat = np.sin(0.5 * (lat - point[1]))
		sin_half_lon = np.sin(0.5 * (lon - point[0]))
		haversine = sin_half_lat ** 2 + np.cos(lat) * np.cos(point[1]) * sin_half_lon ** 2
		return float(np.min(2.0 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))))
	return float(np.min(np.hypot(points[:, 0] - point[0], points[:, 1] - point[1])))


def _has_repeated_vertex(loop: np.ndarray, is_sky: bool) -> bool:
	'''
	True when a chained loop revisits a vertex — two loops pinched at crossing
	points, e.g. the lunes of a symmetric difference. A polygon cannot represent
	that.

	Coincidence is judged per vertex against its own local spacing, not a global
	scalar, for the same reason the seam test is local.
	'''
	gaps = _separations(loop, is_sky)
	if len(gaps) == 0:
		return False
	medians = _local_medians(gaps)
	tolerances = np.concatenate([medians, medians[-1:]]) # one per vertex
	count = len(loop)

	# pairwise separations, vectorized: a python double loop over a 400-point
	# mesh is thousands of AST-free but still needless calls
	if is_sky:
		lon, lat = loop[:, 0], loop[:, 1]
		sin_half_lat = np.sin(0.5 * (lat[:, None] - lat[None, :]))
		sin_half_lon = np.sin(0.5 * (lon[:, None] - lon[None, :]))
		haversine = sin_half_lat ** 2 + np.cos(lat)[:, None] * np.cos(lat)[None, :] * sin_half_lon ** 2
		distances = 2.0 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))
	else:
		distances = np.hypot(loop[:, 0][:, None] - loop[:, 0][None, :],
		                     loop[:, 1][:, None] - loop[:, 1][None, :])

	pair_tolerance = np.minimum(tolerances[:, None], tolerances[None, :])
	coincident = distances <= pair_tolerance
	# consecutive vertices are legitimately close, and the loop's own closure
	# (first vertex against last) is expected: mask the three bands out
	rows, columns = np.indices(distances.shape)
	coincident &= (columns - rows) >= 2
	coincident[0, count - 1] = False
	return bool(coincident.any())


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
		A polygon whose vertices lie exactly on this compound's boundary.

		The compound's boundary mesh is split at its seams into per-component
		arcs, the arcs are chained into closed loops, and the single resulting
		loop becomes the polygon. Every step that could produce a plausible but
		WRONG polygon raises instead:

		* an empty compound (``a & b`` for disjoint operands, ``a - a``,
		  ``a ^ a``) raises :class:`ValueError`. This matters more than it
		  sounds: one of the two empty classes returns a finite bounding disc of
		  about the operand's own scale and a full, ordinary-looking boundary
		  mesh, so an unguarded algorithm would polygonize ``a ^ a`` into a
		  plausible circle with wrong membership everywhere;
		* a boundary with more than one loop — disjoint pieces, or a hole from
		  ``a - b`` with ``b`` inside ``a`` — raises
		  :class:`~cornish.exc.RegionNotConnected`, as does a loop pinched at a
		  crossing point (the lunes of a symmetric difference). Use
		  :meth:`toMoc` for coverage a polygon cannot represent;
		* sampling too coarse to resolve the seams raises :class:`ValueError`,
		  never a silently cut waist or a silently closed gap.

		NOTE: ``MeshSize`` applies PER COMPONENT, so the vertex count scales with
		the number of leaves — an n-leaf union at ``npoints=200`` yields roughly
		n x 200 boundary points.

		Stated resolution limit: components separated by less than the chaining
		tolerance cannot be distinguished from connected ones at the chosen
		``npoints``; finite boundary sampling cannot resolve sub-tolerance
		topology, and the scale shrinks as 1/npoints.

		:param npoints: approximate boundary points sampled PER COMPONENT
		:param maxerr: accepted and ignored (interface uniformity with the
			abstract signature) — the vertices lie exactly on the AST boundary,
			so ``npoints`` alone governs fidelity
		'''
		# --- step 0: emptiness, in three layers (§4.3) ---
		centre, radius = self.astObject.getregiondisc()
		if _degenerate_disc(centre, radius):
			raise ValueError(
				"this compound region is empty (e.g. an intersection of disjoint regions); "
				"an empty region has no boundary polygon")
		# layers 2 and 3: the guarded SURFACE-mesh probe and the belt that its
		# own interior point must test as a member. The probe is fetched exactly
		# once here and REUSED as step 6's containment-fix point.
		try:
			probe = self._interiorProbePoint()
		except ValueError as e:
			raise ValueError(
				"this compound region is empty (e.g. a region differenced with itself); "
				"an empty region has no boundary polygon") from e

		# --- step 1: the boundary mesh, and its geometry in frame units ---
		mesh_degrees = np.asarray(self.boundaryPointMesh(npoints=npoints), dtype=float)
		if len(mesh_degrees) < 3:
			raise ValueError(
				f"this compound region's boundary mesh has only {len(mesh_degrees)} point(s); "
				f"it is degenerate or empty and has no boundary polygon")
		# the bridge speaks AST's (naxes, npoints) layout; the geometry below wants
		# one point per row
		mesh_frame = np.asarray(bridge.to_frame_units(mesh_degrees, self.astObject), dtype=float).T
		is_sky = bridge.is_sky(self.astObject)

		try:
			# --- steps 2 and 3: arcs, then loops ---
			arcs = _split_into_arcs(mesh_frame, is_sky)
			loops = _chain_arcs(arcs, is_sky)

			# --- step 4: reject what a polygon cannot represent ---
			# ...but only ASSERT a topological fact when the sample can support
			# one: every tolerance here is a local median over an 11-gap window,
			# so an arc with fewer gaps than that has no local scale and a
			# multi-loop verdict derived from it is not trustworthy. Saying
			# "increase npoints" is honest; claiming the region is disconnected
			# when it is not would be a confident wrong answer.
			under_resolved = any(len(arc) <= _LOCAL_MEDIAN_WINDOW for arc in arcs)
			if len(loops) != 1:
				if under_resolved:
					raise ValueError(
						f"this compound region's boundary resolved into {len(loops)} loops, but at "
						f"npoints={npoints} some arcs are too short to judge that reliably; "
						f"increase npoints")
				raise RegionNotConnected(
					f"this compound region's boundary has {len(loops)} separate loops (disjoint "
					f"pieces or holes); a polygon cannot represent it — use toMoc()")
			loop = loops[0]
			if _has_repeated_vertex(loop, is_sky):
				if under_resolved:
					raise ValueError(
						f"this compound region's boundary loop appears pinched, but at "
						f"npoints={npoints} some arcs are too short to judge that reliably; "
						f"increase npoints")
				raise RegionNotConnected(
					"this compound region's boundary loop revisits a vertex (two pieces pinched "
					"at a crossing point); a polygon cannot represent it — use toMoc()")

			# --- step 5: the missed-seam backstop ---
			# at coarse npoints the step-2 ratio test can miss a seam, and a missed
			# seam becomes a polygon edge cutting straight across the region (e.g.
			# across a union's waist). Mesh spacing is uniform WITHIN a component's
			# arc, so no legitimate edge of the assembled loop can be a large
			# multiple of its neighbours.
			edges = _separations(loop, is_sky)
			if len(edges):
				medians = _local_medians(edges)
				positive = medians[medians > 0]
				floor = float(np.min(positive)) if len(positive) else 0.0
				limits = SEAM_FACTOR * np.where(medians > 0, medians, floor)
				if np.any(edges > limits):
					raise ValueError(
						"boundary sampling too coarse to polygonize this compound reliably; "
						"increase npoints")
		except Ast.AstError as e:
			# no raw pyast internal error leaves this method (§4.3 step 0)
			raise ValueError("this compound region's boundary could not be polygonized") from e

		# --- step 6: build, with the winding fixed against the interior probe ---
		loop_degrees = np.asarray(bridge.from_frame_units(loop.T, self.astObject), dtype=float)
		polygon = self._orientedPolygon(loop_degrees, probe=probe)

		# --- step 7: the post-condition, because the steps above are heuristics ---
		self._verifyPolygonAgrees(polygon, npoints, loop, is_sky)
		return polygon

	def _verifyPolygonAgrees(self, polygon: ASTPolygon, npoints:int,
	                         boundary: np.ndarray, is_sky: bool) -> None:
		'''
		Refuse a polygon that does not agree with this compound on both sides.

		Steps 2-5 are heuristics over a finite sample, and a seam missed at coarse
		``npoints`` produces a polygon that is plausible and wrong — at
		``npoints=8`` on an ordinary two-circle union, one that disagrees with the
		compound almost everywhere. The heuristics' output is therefore CHECKED
		rather than trusted: every interior mesh point of this compound must lie
		inside the polygon, and every interior mesh point of its complement must
		lie outside.

		Stated blind spot: agreement is checked at a finite sample, so a
		disagreement confined to a region smaller than the mesh spacing can pass.
		That is the same resolution limit ``npoints`` sets everywhere in this
		method, not an additional one.

		:raises ValueError: the polygon and the compound disagree
		'''
		# Points sitting ON the boundary are not evidence of anything: AST's
		# interior meshes legitimately place points in the cusps where two
		# components meet, and at those the polygon's own edge decides membership
		# by floating-point luck. Only disagreements further than one seam-scale
		# from the polygon's vertices count.
		gaps = _separations(boundary, is_sky)
		near_boundary = SEAM_FACTOR * (float(np.median(gaps)) if len(gaps) else 0.0)

		# ...and if that band is as wide as the region itself, EVERY point is
		# "on the boundary", the check below can see nothing, and the polygon is
		# unverifiable rather than verified. That is precisely the coarse-sampling
		# case that produces a plausible wrong polygon, so it raises here.
		_, disc_radius = self.astObject.getregiondisc()
		if near_boundary >= disc_radius:
			raise ValueError(
				f"boundary sampling too coarse to polygonize this compound reliably: at "
				f"npoints={npoints} the boundary is sampled no more finely than the region's own "
				f"size, so the result cannot be verified — increase npoints")

		def disagreements(points, expected_inside: bool) -> tuple:
			''' (real disagreements, disagreements excused as boundary-band) '''
			real = excused = 0
			for point in points:
				if polygon.pointInRegion(list(point)) is expected_inside:
					continue
				frame_point = np.asarray(bridge.to_frame_units(point, self.astObject, squeeze=True),
				                         dtype=float)
				if _min_separation(frame_point, boundary, is_sky) <= near_boundary:
					excused += 1 # on the boundary at this resolution; says nothing
				else:
					real += 1
			return real, excused

		interior = np.asarray(self.interiorPointMesh(npoints=npoints), dtype=float)
		missed, missed_excused = disagreements(interior, True)

		complement = ~self
		try:
			exterior = np.asarray(complement.interiorPointMesh(npoints=npoints), dtype=float)
		except ValueError as e:
			if not isinstance(e.__cause__, Ast.AstError):
				# only AST's own "no mesh here" answer is expected and handled;
				# anything else (a bad npoints, a non-finite coordinate) is a real
				# failure and must not be absorbed by the fallback below
				raise
			# AST cannot mesh an unbounded complement (a compound on a basic
			# frame): the ring of points outside the compound's own bounding disc
			# is always available and probes the same thing
			centre, radius = self.astObject.getregiondisc()
			frame_points, _ = _disc_escape_points(self.astObject, centre, radius)
			exterior = np.asarray(bridge.from_frame_units(np.asarray(frame_points).T, self.astObject),
			                      dtype=float)
		included, included_excused = disagreements(exterior, False)

		if missed or included:
			raise ValueError(
				f"the polygon built from this compound region's boundary does not agree with "
				f"the region ({missed} of {len(interior)} interior points fall outside it, "
				f"{included} of {len(exterior)} exterior points fall inside it); the boundary "
				f"sampling is too coarse to polygonize this compound — increase npoints")

		# The boundary band above EXCUSES disagreements, so excused disagreements
		# must not dominate — and the test is PER SIDE, because the two sides
		# answer different questions and either alone can be uninformative. A
		# polygon that is wrong about the entire interior but happens to be right
		# about the exterior sample pools to just under half and would slip past a
		# combined threshold; each side has to be meaningful on its own, or
		# "verified" is the silent-skip this check exists to prevent.
		unverified = []
		if len(interior) and missed_excused * 2 > len(interior):
			unverified.append(f"{missed_excused} of {len(interior)} interior")
		if len(exterior) and included_excused * 2 > len(exterior):
			unverified.append(f"{included_excused} of {len(exterior)} exterior")
		if unverified:
			raise ValueError(
				f"the polygon built from this compound region's boundary could not be verified: "
				f"{' and '.join(unverified)} probe points disagree with the region but lie within "
				f"one seam-width of the boundary, where the sampling cannot decide — "
				f"increase npoints")
