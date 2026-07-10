
import logging
from typing import Iterable, Union, List

import numpy as np

import astropy
import astropy.units as u
import starlink.Ast as Ast

from .box import ASTBox
from .circle import ASTCircle
from .region import ASTRegion
from .polygon import ASTPolygon
from ..enums import RegionOperation

__all__ = ["ASTCompoundRegion"]

logger = logging.getLogger("cornish")

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

		import operator
		if operation is None:
			raise ValueError("An 'operation' must be specified (one of RegionOperation.AND/OR/XOR).")
		if isinstance(operation, bool):
			raise TypeError("'operation' must be a RegionOperation (or Ast.AND/OR/XOR integer), not a bool.")
		try:
			operation = operator.index(operation) # accepts NumPy ints; rejects floats
		except TypeError:
			raise TypeError(f"'operation' must be a RegionOperation (or Ast.AND/OR/XOR integer); got '{type(operation).__name__}'.") from None
		operation = RegionOperation(operation) # ValueError for out-of-range ints

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

	def componentRegions(self) -> List[ASTRegion]:
		'''
		Returns a list of region objects that comprise this component region.
		'''
		(map1, map2, series, invert1, invert2) = self.astObject.decompose()
		regions = list()
		for r in [map1, map2]:
			if r.isapolygon():
				regions.append(ASTPolygon(ast_object=r))
			elif r.isabox():
				regions.append(ASTBox(ast_object=r))
			elif r.isacircle():
				regions.append(ASTCircle(ast_object=r))
			elif r.isacmpregion():
				regions.append(ASTCompoundRegion(ast_object=r))
			#elif r.isaellipse():
			#	regions.append(ASTEllipse(r))
			#elif r.isanullregion():
			#	regions.append( ... ?)
			else:
				raise NotImplementedError(f"A region type was found that is not yet handled: {type(r)}")
		return regions

	def toPolygon(self) -> ASTPolygon: #, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec) -> ASTPolygon:
		'''
		Return a single polygon that bounds the total of the regions that make up this compound region.
		'''

		# TODO: check that frame match; deal with otherwise!!
		logging.warning("ASTCompoundRegion currently doesn't check that component regions are in the same frame!")

		region1, region2 = self.componentRegions()
		if region1.overlaps(region2):
			raise NotImplementedError("")
		else:

			# this doesn't have to be ra/dec; it's the first and second axis
			r1_ra  = region1.points.T[0]
			r1_dec = region1.points.T[1]
			r2_ra  = region2.points.T[0]
			r2_dec = region2.points.T[1]

			points = np.array([np.concatenate((r1_ra,r2_ra)), np.concatenate((r1_dec,r2_dec))])

			# .. todo:: convert points from one frame to the other

			# create bounding polygon
			polygon = ASTPolygon(frame=region1.frame, points=np.array([
			    [min(points[0]), min(points[1])],
			    [min(points[0]), max(points[1])],
			    [max(points[0]), max(points[1])],
			    [max(points[0]), min(points[1])]
			]))

			if not polygon.pointInRegion(region1.points[0]):
				polygon.negate()

			return polygon
