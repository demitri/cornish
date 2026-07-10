
from __future__ import annotations # remove in Python 3.10

import logging
from typing import Iterable, Union

import numpy as np
import astropy
import astropy.units as u
import starlink.Ast as Ast

from .region import ASTRegion
from ..exc import NotASkyRegion
from .._validation import as_integer

__all__ = ['ASTMoc']

logger = logging.getLogger("cornish")

class ASTMoc(ASTRegion):
	'''
	An IVOA Multi-Order Coverage (MOC) map: a HEALPix-based description of an
	arbitrary patch (or patches) of sky. A MOC describes its own cell set exactly,
	but represents an arbitrary input region only to the resolution of its
	smallest cells (set by ``max_order``).

	A MOC is a full :class:`ASTRegion`, so containment tests, overlap logic,
	bounding circles, meshes, and plotting all work on it directly. Because any
	region — including compound regions — can be added to a MOC, converting to a
	MOC is the recommended way to flatten complicated compound geometry
	(e.g. the union of many survey field footprints) into a single well-behaved object
	that can also be serialized to the standard IVOA string and FITS forms
	understood by Aladin, MOCpy, VizieR, pgsphere (``smoc``), etc.

	The ``max_order`` parameter sets the HEALPix resolution of the map; each
	increment halves the cell scale. Order 10 corresponds to ~3.4 arcmin cells,
	order 15 to ~6.4 arcsec, order 18 to ~0.8 arcsec.

	Accepted signatures:

	.. code-block:: python

	    moc = ASTMoc()                      # empty MOC, AST default resolution
	    moc = ASTMoc(max_order=12)          # empty MOC at the given resolution
	    moc = ASTMoc(ast_object=ast_object) # wrap an existing starlink.Ast.Moc

	:param ast_object: an existing :class:`starlink.Ast.Moc` object to wrap
	:param max_order: HEALPix order of the finest cells used, 0..27
	'''
	def __init__(self, ast_object:Ast.Moc=None, max_order:int=None):
		if ast_object is not None:
			if max_order is not None:
				raise ValueError("Cannot specify both 'ast_object' and 'max_order'.")
			if isinstance(ast_object, Ast.Moc):
				super().__init__(ast_object=ast_object)
				return
			else:
				raise TypeError(f"The 'ast_object' provided was not of type starlink.Ast.Moc (got '{type(ast_object)}').")

		if max_order is None:
			super().__init__(ast_object=Ast.Moc())
		else:
			max_order = as_integer(max_order, "max_order")
			if not (0 <= max_order <= 27):
				raise ValueError(f"'max_order' must be in the range 0..27 (got {max_order}).")
			super().__init__(ast_object=Ast.Moc(f"MaxOrder={max_order}"))

	@classmethod
	def fromRegion(cls, region:Union[ASTRegion, Ast.Region], max_order:int=10) -> ASTMoc:
		'''
		Create a MOC covering the provided region.

		:param region: the region to cover; must be defined in a sky frame
		:param max_order: HEALPix order of the finest cells used (see class docstring)
		'''
		moc = cls(max_order=max_order)
		moc.add(region)
		return moc

	@classmethod
	def fromRegions(cls, regions:Iterable[Union[ASTRegion, Ast.Region]], max_order:int=10) -> ASTMoc:
		'''
		Create a MOC covering the union of all provided regions,
		e.g. the total footprint of a collection of images.

		:param regions: regions to cover; each must be defined in a sky frame
		:param max_order: HEALPix order of the finest cells used (see class docstring)
		'''
		moc = cls(max_order=max_order)
		for region in regions:
			moc.add(region)
		return moc

	@classmethod
	def fromString(cls, moc_string:str) -> ASTMoc:
		'''
		Create a MOC from its IVOA string serialization (ASCII form, e.g.
		``"6/2463 7/9847-9850"``, or the equivalent JSON form).

		:param moc_string: MOC in IVOA ASCII or JSON encoding
		'''
		if not isinstance(moc_string, str):
			raise TypeError(f"ASTMoc.fromString expects a string (got '{type(moc_string).__name__}').")
		moc = cls()
		try:
			moc.astObject.addmocstring(moc_string)
		except Ast.AstError as e:
			# pyast raises library-specific errors (e.g. Ast.INMOC) for malformed input;
			# surface them as a ValueError with the original error chained
			raise ValueError(f"The provided string could not be parsed as a MOC: '{moc_string}'") from e
		return moc

	def add(self, region:Union[ASTRegion, Ast.Region], operation:int=Ast.OR):
		'''
		Combine the given region with this MOC (by default, the union — the region's
		coverage is added to this map's).

		:param region: region to combine; must be defined in a sky frame
		:param operation: one of ``starlink.Ast.OR`` (default), ``AND``, ``XOR``
		'''
		if isinstance(region, ASTRegion):
			ast_region = region.astObject
		elif isinstance(region, Ast.Region):
			ast_region = region
		else:
			raise TypeError(f"ASTMoc.add: expected an ASTRegion or starlink.Ast.Region (got '{type(region)}').")
		try:
			self.astObject.addregion(ast_region, operation)
		except Ast.NOSKY as e:
			raise NotASkyRegion(f"Only regions defined in a sky frame can be added to a MOC (got a region in a '{ast_region.getregionframe().Class}' frame).") from e

	@property
	def maxOrder(self) -> Union[int, None]:
		'''
		The HEALPix order of the finest cells in this MOC.

		Returns ``None`` for an empty MOC created without an explicit ``max_order``
		(the order is not determined until a region is added; AST reports this
		internally with a sentinel value of -1).
		'''
		max_order = int(self.astObject.MaxOrder)
		return None if max_order == -1 else max_order

	@property
	def resolution(self) -> astropy.units.Quantity:
		''' The best resolution of this MOC (the scale of its smallest cells). '''
		return float(self.astObject.MaxRes) * u.arcsec

	@property
	def area(self) -> astropy.units.Quantity:
		''' The sky area covered by this MOC. '''
		# MocArea is in square arcminutes
		return (float(self.astObject.MocArea) * u.arcmin * u.arcmin).to(u.deg * u.deg)

	def toString(self, json:bool=False) -> str:
		'''
		Return the IVOA string serialization of this MOC.

		:param json: if True, return the JSON encoding rather than the ASCII encoding
		'''
		# note: pyast's getmocstring rejects keyword arguments; pass positionally
		return self.astObject.getmocstring(json)

	def toPolygon(self, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec) -> "ASTPolygon":
		'''
		Not yet implemented for MOCs.

		A MOC boundary is a staircase of HEALPix cell edges and the points of its
		boundary mesh are not ordered, so they cannot be used directly as polygon
		vertices (doing so was verified to produce an invalid polygon). A correct
		implementation needs to trace the boundary contour in order first.
		Until then, prefer using the MOC itself: it already supports containment,
		overlap, area, and bounding-circle operations.
		'''
		raise NotImplementedError(
			"ASTMoc.toPolygon is not yet implemented: the MOC boundary mesh is unordered "
			"and would produce an invalid polygon. Use the ASTMoc directly for containment, "
			"overlap, area, and boundingCircle operations."
		)
