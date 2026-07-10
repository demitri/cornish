
'''
Notes.

From the AST documentation: "The Box class does not define any new attributes beyond those
which are applicable to all Regions.", i.e. for AST, there is nothing special about a box
beyond a specific means to define it, i.e. a corner and center or two corner points.
'''

from __future__ import annotations # remove in Python 3.10

import logging
from typing import Union, Iterable

import numpy as np
import starlink
import starlink.Ast as Ast
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord

from .region import ASTRegion, DEFAULT_UNCERTAINTY
from ..mapping import ASTMapping
from ..mapping import ASTFrame, ASTFrameSet
from .. import _pyast_bridge as bridge

__all__ = ["ASTBox"]

CENTRE_CORNER = 0
CORNER_CORNER = 1

logger = logging.getLogger("cornish")

class ASTBox(ASTRegion):
	'''
	ASTBox is an ASTRegion that represents a box with sides parallel to the axes of an ASTFrame.

	Boxes are created via the factory class methods:

	.. code-block:: python

		b = ASTBox.fromCentreAndCorner(frame, centre=..., corner=...)
		b = ASTBox.fromCorners(frame, corners=(c1, c2))
		b = ASTBox(ast_box)   # where ast_box is an existing Ast.Box object

	Points can be given in any of the bridge-accepted forms — SkyCoord, Quantity
	pairs, or two-element containers of bare values (read as degrees on sky
	frames, native units otherwise), e.g.

	.. code-block:: python

		(1,2)
		[1,2]
		np.array([1,2])

	The 'frame' parameter may be a cornish or ``starlink.Ast`` frame, frame set
	(the current frame governs), or region.

	A Box is similar to an Interval, the only real difference being that the Interval
	class allows some axis limits to be unspecified. Note, a Box will only look like a box
	if the Frame geometry is approximately flat. For instance, a Box centered close to a pole
	in a SkyFrame will look more like a fan than a box (the Polygon class can be used to
	create a box-like region close to a pole).

	:param ast_object: an existing object of type :class:`starlink.Ast.Box`
	'''
	def __init__(self, ast_object:starlink.Ast.Box=None):

		if not isinstance(ast_object, Ast.Box):
			raise TypeError(f"An ASTBox can only be initialized with a starlink.Ast.Box object (got '{type(ast_object)}'); use the fromCentreAndCorner/fromCorners factory methods to create new boxes.")

		super().__init__(ast_object=ast_object)

	@classmethod
	def fromCentreAndCorner(cls, frame:Union[Ast.Frame, ASTFrame], centre:Iterable=None, corner:Iterable=None, center:Iterable=None,
	                        uncertainty:Union[None, float, u.Quantity, ASTRegion, Ast.Region]=DEFAULT_UNCERTAINTY) -> ASTBox:
		'''
		Create a new ASTBox object defined by the provided corner and centre points.

		:param frame: the frame the provided points lie in; accepts cornish or ``starlink.Ast`` frames, frame sets (current frame governs), or regions
		:param centre: the coordinate of the point at the centre of the box in the frame provided
		:param corner: the coordinate of the point at any corner of the box in the frame provided
		:param center: synonym for 'centre', ignored if 'centre' is defined
		:param uncertainty: uncertainty of the boundary: a Quantity (angular, sky frames only), a number (degrees on sky frames), or a Region; omit for the 1 arcsec default (sky frames), or pass ``None`` explicitly for AST's internal default
		'''
		if center is not None and centre is None:
			centre = center

		if frame is None:
			raise ValueError("A frame the region is defined in must be provided.")
		elif not all([x is not None for x in [centre, corner]]):
			raise ValueError("Both a corner and centre coordinates must be provided to define this box region.")

		f = bridge._unwrap(frame) # frame sets and regions are legal: the current/encapsulated frame governs
		c = bridge.to_frame_units(centre, f, squeeze=True)
		k = bridge.to_frame_units(corner, f, squeeze=True)

		if uncertainty is DEFAULT_UNCERTAINTY and not bridge.is_sky(f):
			uncertainty = None # the angular 1 arcsec default is meaningless on a basic frame; AST's internal default governs
		unc_region = bridge.as_uncertainty_region(uncertainty, f, c)

		# pyast SILENTLY IGNORES the `unc=` keyword on region constructors and
		# rejects a positional None, so the uncertainty region must be appended
		# POSITIONALLY and only when non-None (SPEC-04A §8 / D15). Do not "clean
		# this up" into the kwarg form — it parses fine and silently drops the
		# uncertainty.
		args = (f, CENTRE_CORNER, c, k) + ((unc_region,) if unc_region is not None else ())
		box = ASTBox(ast_object=Ast.Box(*args))
		box._uncertainty = uncertainty
		return box

	@classmethod
	def fromCorners(cls, frame:Union[Ast.Frame, ASTFrame], corners:Iterable[Iterable]=None,
	                uncertainty:Union[None, float, u.Quantity, ASTRegion, Ast.Region]=DEFAULT_UNCERTAINTY) -> ASTBox:
		'''
		Create a new ASTBox object defined by two corner points.

		:param frame: the frame the provided points lie in; accepts cornish or ``starlink.Ast`` frames, frame sets (current frame governs), or regions
		:param corners: a collection (list, tuple, array, etc.) of coordinates of two opposite corners of the box in the frame provided
		:param uncertainty: uncertainty of the boundary: a Quantity (angular, sky frames only), a number (degrees on sky frames), or a Region; omit for the 1 arcsec default (sky frames), or pass ``None`` explicitly for AST's internal default
		'''
		if frame is None:
			raise ValueError("A frame the region is defined in must be provided.")
		if corners is None:
			raise ValueError("Two corner coordinates must be provided to define this box region.")

		c1 = corners[0]
		c2 = corners[1]

		if not all([x is not None for x in [c1, c2]]):
			raise ValueError("Two corner coordinates must be provided to define this box region.")

		f = bridge._unwrap(frame) # frame sets and regions are legal: the current/encapsulated frame governs
		p1 = bridge.to_frame_units(c1, f, squeeze=True)
		p2 = bridge.to_frame_units(c2, f, squeeze=True)

		if uncertainty is DEFAULT_UNCERTAINTY and not bridge.is_sky(f):
			uncertainty = None # the angular 1 arcsec default is meaningless on a basic frame; AST's internal default governs
		unc_region = bridge.as_uncertainty_region(uncertainty, f, p1)

		# pyast SILENTLY IGNORES the `unc=` keyword on region constructors and
		# rejects a positional None, so the uncertainty region must be appended
		# POSITIONALLY and only when non-None (SPEC-04A §8 / D15). Do not "clean
		# this up" into the kwarg form — it parses fine and silently drops the
		# uncertainty.
		args = (f, CORNER_CORNER, p1, p2) + ((unc_region,) if unc_region is not None else ())
		box = ASTBox(ast_object=Ast.Box(*args))
		box._uncertainty = uncertainty
		return box

	@property
	def centre(self) -> np.ndarray:
		'''
		Returns the location of the Box's centre as a coordinate pair, in degrees if a sky frame.

		:returns: a Numpy array of points (axis1, axis2)
		'''
		# 'getregionpoints' returns two points: (the centre, a corner)
		centre_point, corner_point = bridge.from_frame_units(self.astObject.getregionpoints(), self.astObject)
		return centre_point

	@property
	def center(self) -> np.ndarray:
		'''
		Returns "self.centre". This is a British library, after all.
		'''
		return self.centre

	@property
	def corner(self) -> np.ndarray:
		'''
		Returns the location of one of the box's corners as a coordinate pair, in degrees if a sky frame.

		:returns: a Numpy array of points (axis1, axis2)
		'''
		centre_point, corner_point = bridge.from_frame_units(self.astObject.getregionpoints(), self.astObject)
		return corner_point

	def corners(self, mapping=None) -> np.ndarray:
		'''
		Returns all four corners of the box, mapped through the provided frame set.

		The corners are computed in this box's own frame from its centre/corner
		definition, transformed with the frame set's forward mapping, and
		returned in the frame set's current frame units (degrees for sky).

		:param mapping: an :class:`ASTFrameSet` (or ``starlink.Ast.FrameSet``) whose base frame matches this box's frame
		:returns: an array of four points, shape (4, 2)
		'''
		if mapping is None:
			raise ValueError("A frame set must be provided to return a list of corner points.")
		ast_mapping = bridge._unwrap(mapping)
		if not isinstance(ast_mapping, Ast.FrameSet):
			raise TypeError(f"The 'mapping' parameter must be a frame set (an object with base and current frames); got '{type(mapping)}'.")

		base = ast_mapping.getframe(Ast.BASE)       # live pointers, used read-only (V2)
		current = ast_mapping.getframe(Ast.CURRENT)

		# build the four corners in this box's own frame
		centre_point, corner_point = bridge.from_frame_units(self.astObject.getregionpoints(), self.astObject)
		d_x, d_y = corner_point - centre_point
		corner_points = np.array([
			[centre_point[0] - d_x, centre_point[1] - d_y],
			[centre_point[0] - d_x, centre_point[1] + d_y],
			[centre_point[0] + d_x, centre_point[1] + d_y],
			[centre_point[0] + d_x, centre_point[1] - d_y]
		])

		pts = bridge.to_frame_units(corner_points, base)
		out = ast_mapping.tran(pts, True)
		return bridge.from_frame_units(out, current) # norm owns any wrapping

	def toPolygon(self, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec) -> ASTPolygon:
		'''
		Returns a four-vertex ASTPolygon that describes this box in the same frame.

		The parameters 'npoints' and 'maxerr' are ignored.
		'''
		# Note: other region methods use the boundary mesh points technique to get a polygon.
		# A box is a simple enough shape to get a precise polygon from the provided points.

		from .polygon import ASTPolygon # avoid circular import

		centre_point, corner_point = bridge.from_frame_units(self.astObject.getregionpoints(), self.astObject)
		d_x, d_y = centre_point - corner_point

		polygon_points = np.array([
			[centre_point[0] - d_x, centre_point[1] - d_y],
			[centre_point[0] - d_x, centre_point[1] + d_y],
			[centre_point[0] + d_x, centre_point[1] + d_y],
			[centre_point[0] + d_x, centre_point[1] - d_y]
		])

		polygon = ASTPolygon(frame=self.frame(), points=polygon_points)
		if polygon.containsPoint(centre_point) is False:
			polygon.negate()
		polygon.wcs = self.wcs # propagate the originating WCS, when known (SPEC-04 §2)
		return polygon

	@property
	def area(self) -> u.Quantity:
		'''
		The area of the box within its frame (e.g. on a Cartesian plane or sphere).
		'''
		frame = self.frame() # create variable here as frame() creates a copy

		if frame.isSkyFrame:
			# Girard's theorem on the box's four-vertex polygon
			return self.toPolygon().area.to(u.deg * u.deg)

		elif frame.system == 'Cartesian' and frame.domain == "GRID":

			centre_point, corner_point = self.astObject.getregionpoints().T # native units; no conversion
			half_dx, half_dy = centre_point - corner_point
			return abs(2 * half_dx) * abs(2 * half_dy) * u.pixel * u.pixel

		else:
			raise NotImplementedError(f"The area computation and units for a frame with .system='{frame.system}' and/or .domain='{frame.domain}' has not been written.")
