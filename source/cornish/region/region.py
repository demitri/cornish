
from __future__ import annotations # remove in Python 3.10
# Needed for forward references, see:
# https://stackoverflow.com/a/33533514/2712652

import logging
from abc import ABCMeta, abstractproperty, abstractmethod
from typing import Union, Iterable, Tuple

import math

import numpy as np
import astropy
import astropy.units as u
import starlink
import starlink.Ast as Ast

import cornish.region # to avoid circular imports below - better way?
from ..mapping import ASTFrame, ASTFrameSet, ASTMapping
from ..exc import NotA2DRegion, CoordinateSystemsCouldNotBeMapped
from ..enums import MeshType, OverlapType
from .._validation import as_integer
from .. import _pyast_bridge as bridge

#: Module-level default for region constructors' ``uncertainty`` parameter.
#: Identity against this singleton distinguishes "parameter omitted" (the 1
#: arcsec default is materialized and delivered to AST on sky frames; on basic
#: frames an angular default is meaningless and AST's internal default governs)
#: from an explicitly-passed value (SPEC-04A §8 three-way rule).
DEFAULT_UNCERTAINTY = 1 * u.arcsec

#: Ceiling for ``meshSize``. Physically anchored rather than the C-int limit:
#: 10^7 points sample even a full great-circle boundary (the longest possible,
#: 360 deg) at ~0.13 arcsec spacing — finer than the 1 arcsec default boundary
#: uncertainty, below which additional mesh points are meaningless — while a
#: mesh this size already costs ~160 MB and minutes of compute. Larger requests
#: are almost certainly errors (a value beyond C int is silently rewritten by
#: AST besides). Adjust here, deliberately, if a use case ever needs more.
MAX_MESH_SIZE = 10_000_000


__all__ = ['ASTRegion']

logger = logging.getLogger("cornish")

'''
Copied from documentation, to be implemented.

Properties of ASTRegion over those from ASTFrame

  * Adaptive: Should the area adapt to changes in the coordinate system?
  * Negated: Has the original region been negated?
  * Closed: Should the boundary be considered to be inside the region?
  * MeshSize: Number of points used to create a mesh covering the Region
  * FillFactor: Fraction of the Region which is of interest
  * Bounded: Is the Region bounded?

Region-specific methods:

  * astGetRegionBounds: Get the bounds of a Region
  * astGetRegionFrame: Get a copy of the Frame represent by a Region
  * astGetRegionFrameSet: Get a copy of the Frameset encapsulated by a Region
  * astGetRegionMesh: Get a mesh of points covering a Region (points returned are normalized)
  * astGetRegionPoints: Get the positions that define a Region
  * astGetUnc: Obtain uncertainty information from a Region
  * astMapRegion: Transform a Region into a new coordinate system
  * astNegate: Toggle the value of the Negated attribute
  * astOverlap: Determines the nature of the overlap between two Regions
  * astMask<X>: Mask a region of a data grid
  * astSetUnc: Associate a new uncertainty with a Region
  * astShowMesh: Display a mesh of points on the surface of a Region

'''

class ASTRegion(ASTFrame, metaclass=ABCMeta):
	'''
	Represents a region within a coordinate system.
	This is an abstract superclass - there is no supported means to create an ASTRegion object directly
	(see :py:class:`ASTBox`, :py:class:`ASTPolygon`, etc.).

	Accepted signatures for creating an ASTRegion:

	.. code-block:: python

		r = ASTRegion(ast_object)

	:param ast_object:
	:param uncertainty:
	'''
	def __init__(self, ast_object:starlink.Ast.Region=None, uncertainty=None):
		super().__init__(ast_object=ast_object)
		self._uncertainty = uncertainty
		#: The pixel<->world :class:`ASTFrameSet` this region was born from, when
		#: known (e.g. set by ``fromFITSHeader``); ``None`` otherwise (SPEC-04 §2).
		self.wcs = None

	def __add__(self, other_region):
		# TODO: check data type, handle both ASTRegion and the ast_object region?
		from .compound_region import ASTCompoundRegion # forward import to avoid circular import error
		return ASTCompoundRegion(regions=[self, other_region], operation=Ast.AND)

	@classmethod
	def fromFITSHeader(cls, fits_header=None, maxerr:astropy.units.Quantity=1*u.arcsec, maxvert:int=200) -> ASTRegion:
		'''
		Factory method to create a region covering the field described by the
		provided FITS header; the returned object will be as specific as
		possible (currently always an :py:class:`ASTPolygon`).

		This is a dispatcher over the shape-specific factories (SPEC-04 §3.4):
		today it returns :meth:`ASTPolygon.fromFITSHeader`; a future heuristic
		may return e.g. an :class:`ASTCircle` where one describes the field
		exactly. The returned region carries its pixel<->world frame set in
		``.wcs``.

		:param fits_header: a FITS header (astropy, fitsio, an array of cards, a dict, or a string)
		:param maxerr: maximum deviation of the polygon from the true field outline
		:param maxvert: maximum number of vertices in the returned polygon
		'''
		if fits_header is None:
			raise ValueError("This method requires a 'fits_header' to be set.")

		from .polygon import ASTPolygon # avoid circular import
		return ASTPolygon.fromFITSHeader(header=fits_header, maxerr=maxerr, maxvert=maxvert)

	@classmethod
	def fromSTCS(cls, stcs:str) -> ASTRegion:
		'''
		Factory method that creates a region from an STC-S description string.

		STC-S is the IVOA standard string serialization of sky regions, e.g. the
		``s_region`` column of ObsCore tables:

		.. code-block::

		    Circle ICRS 30.0 45.0 2.0
		    Polygon ICRS 10 20 10 21 11 21 11 20
		    Union ICRS ( Circle 30 45 2 Circle 32 44 1 )

		:param stcs: an STC-S region description
		:returns: an :class:`ASTRegion` subclass instance matching the description
		'''
		from ..channel.channel_io import ListSource # avoid circular import
		import cornish.region as cr

		if not isinstance(stcs, str):
			raise TypeError(f"ASTRegion.fromSTCS expects a string (got '{type(stcs).__name__}').")

		stcs_channel = Ast.StcsChan(ListSource(stcs), None)
		try:
			ast_object = stcs_channel.read()
		except Ast.AstError as e:
			# pyast raises library-specific errors (e.g. Ast.BADIN) for malformed input;
			# surface them as a ValueError with the original error chained
			raise ValueError(f"The provided string could not be parsed as an STC-S region: '{stcs}'") from e
		if ast_object is None or not isinstance(ast_object, Ast.Region):
			raise ValueError(f"The provided string could not be read as an STC-S region: '{stcs}'")

		# wrap in the appropriate cornish class
		if isinstance(ast_object, Ast.Polygon):
			return cr.ASTPolygon(ast_object=ast_object)
		elif isinstance(ast_object, Ast.Circle):
			return cr.ASTCircle(ast_object=ast_object)
		elif isinstance(ast_object, Ast.Box):
			return cr.ASTBox(ast_object=ast_object)
		elif isinstance(ast_object, Ast.CmpRegion):
			return cr.ASTCompoundRegion(ast_object=ast_object)
		else:
			raise NotImplementedError(f"STC-S produced a region class not yet wrapped by cornish: '{ast_object.Class}'")

	def toMoc(self, max_order:int=10) -> "cornish.region.ASTMoc":
		'''
		Return an IVOA Multi-Order Coverage map (:class:`ASTMoc`) covering this region.

		This region must be defined in a sky frame. A MOC is often the most robust
		way to work with complicated geometry (e.g. compound regions): it has a
		well-defined area, bounding circle, and standard serializations.

		:param max_order: HEALPix order of the finest MOC cells (order 10 ≈ 3.4', order 15 ≈ 6.4")
		'''
		from .moc import ASTMoc # avoid circular import
		return ASTMoc.fromRegion(self, max_order=max_order)

	def toSTCS(self, digits:int=16) -> str:
		'''
		Return the IVOA STC-S string description of this region, e.g. ``Circle ICRS 30 45 2``.

		This is the format used by the ``s_region`` column of IVOA ObsCore tables
		and is understood by TAP services, Aladin, and (via ``regions``) DS9.
		Compound regions serialize as, e.g., ``Union ICRS ( Circle ... Circle ... )``.

		Known limitation (AST 9.3): a box serializes to STC-S ``PositionInterval``,
		which the AST STC-S *reader* fails to parse back — a box's STC-S output is
		valid for external services but does not round-trip through
		:meth:`fromSTCS`. Convert with :meth:`toPolygon` first if a round-trip is needed.

		:param digits: number of significant digits written for coordinate values;
			the default (16) makes the geometry round-trip losslessly through
			:meth:`fromSTCS` (STC-S remains lossy for uncertainty and refpos)
		:raises cornish.exc.SerializationNotPossible: if AST cannot represent this region in STC-S
		'''
		from ..channel.channel_io import ListSink # avoid circular import

		digits = as_integer(digits, "digits")
		if not (1 <= digits <= 17):
			raise ValueError(f"'digits' must be between 1 and 17 (got {digits}).")

		# write a copy so the formatting attribute never leaks into this region
		formatted_copy = self.astObject.copy()
		formatted_copy.set(f"Digits={digits}")

		sink = ListSink()
		stcs_channel = Ast.StcsChan(None, sink)
		n_written = stcs_channel.write(formatted_copy)
		if n_written == 0 or len(sink.lines) == 0:
			from ..exc import SerializationNotPossible
			raise SerializationNotPossible(f"This region ({self.__class__.__name__}) could not be serialized to STC-S.")
		return " ".join(line.strip() for line in sink.lines).strip()

	@property
	def points(self) -> np.ndarray:
		'''
		The array of points that define the region. The interpretation of the points is dependent on the type of shape in question.

		* :class:`ASTBox`: returns two points; the center and a box corner.
		* :class:`Circle`: returns two points; the center and a point on the circumference.
		* :class:`CmpRegion`: no points returned; to get points that define a compound region, call this method on each of the component regions via the method "decompose".
		* :class:`Ellipse`: three points: 1) center, 2) end of one axis, 3) end of the other axis
		* :class:`Interval`: two points: 1) lower bounds position, 2) upper bounds position (reversed when interval is an excluded interval)
		* :class:`NullRegion`: no points
		* :class:`PointList`: positions that the list was constructed with
		* :class:`Polygon`: vertex positions used to construct the polygon
		* :class:`Prism`: no points (see CmpRegion)

		+-----------------------+--------------------------------------------------------------------------------+
		| Class                 | Return Value                                                                   |
		+=======================+================================================================================+
		| :class:`ASTBox`       | returns two points; the center and a box corner.                               |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTCircle`    | returns two points; the center and a point on the circumference.               |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTCmpRegion` | no points returned; to get points that define a compound region,               |
		|                       | call this method on each of the component regions via the method "decompose".  |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTEllipse`   | three points: 1) center, 2) end of one axis, 3) end of the other axis          |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTInterval`  | two points: 1) lower bounds position                                           |
		|                       | 2) upper bounds position (reversed when interval is an excluded interval)      |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTNullRegion`| no points                                                                      |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTPointList` | positions that the list was constructed with                                   |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTPolygon`   | vertex positions used to construct the polygon                                 |
		+-----------------------+--------------------------------------------------------------------------------+
		| :class:`ASTPrism`     | no points (see CmpRegion)                                                      |
		+-----------------------+--------------------------------------------------------------------------------+

		NOTE: points returned reflect the current coordinate system and may be different from the initial construction

		:returns: NumPy array of coordinate points in degrees, shape (ncoord,2), e.g. [[ra1,dec1], [ra2, dec2], ..., [ra_n, dec_n]]
		'''

		# getregionpoints returns data as [[x1, x2, ..., xn], [y1, y2, ..., yn]];
		# the bridge normalizes, converts units, and transposes

		region_points = self.astObject.getregionpoints()

		if self.isNegated:
			# reverse order to reflect the definition of the polygon
			# (negation is region semantics, not units — it stays at this call site)
			region_points = np.fliplr(region_points)

		return bridge.from_frame_units(region_points, self.astObject)

	@property
	def uncertainty(self):
		'''
		The uncertainty associated with the boundary of this region, as requested
		at construction (a Quantity, a number, a Region, or ``None``).

		The uncertainty in any point on the boundary of the region is found by
		shifting the uncertainty region so that it is centered at the boundary
		point being considered; the area it covers then represents the
		uncertainty in the boundary position (assumed the same for all points).

		This Python-side value was actually delivered to AST exactly when the
		object was built through the cornish constructors/factories that accept
		an ``uncertainty`` parameter (and it was not ``None``); for objects that
		wrap an externally-created AST region it is ``None`` and AST's internal
		default governs. (Historical note: before the SPEC-04 migration no
		cornish uncertainty ever reached AST — the ``unc=`` keyword is silently
		ignored by pyast.)
		'''
		return self._uncertainty

	@uncertainty.setter
	def uncertainty(self, unc):
		raise NotImplementedError(
			"The uncertainty of an existing region cannot be changed (no astSetUnc "
			"exposure exists in pyast); pass 'uncertainty' to the region's constructor "
			"or factory method instead."
		)

	@property
	def isAdaptive(self):
		'''
		Boolean attribute that indicates whether the area adapt to changes in the coordinate system.
		'''
		return self.astObject.get("Adaptive") == "1"

	@isAdaptive.setter
	def isAdaptive(self, newValue:bool):
		if newValue in [True, 1, "1"]:
			self.astObject.set("Adaptive=1")
		elif newValue in [False, 0, "0"]:
			self.astObject.set("Adaptive=0")
		else:
			raise ValueError("ASTRegion.isAdaptive property must be one of [True, False, 1, 0].")

	@property
	def isNegated(self):
		''' Boolean attribute that indicates whether the original region has been negated. '''
		return self.astObject.get("Negated") == "1"

	@isNegated.setter
	def isNegated(self, newValue:bool):
		if newValue in [True, 1, "1"]:
			self.astObject.set("Negated=1")
		elif newValue in [False, 0, "0"]:
			self.astObject.set("Negated=0")
		else:
			raise ValueError("ASTRegion.isNegated property must be one of [True, False, 1, 0].")

	def negate(self):
		''' Negate the region, i.e. points inside the region will be outside, and vice versa. '''
		self.astObject.negate()

	@property
	def isClosed(self) -> bool:
		'''
		Boolean attribute that indicates whether the boundary be considered to be inside the region.
		'''
		return self.astObject.get("Closed") == "1"

	@isClosed.setter
	def isClosed(self, newValue:bool):
		if newValue in [True, 1, "1"]:
			self.astObject.set("Closed=1")
		elif newValue in [False, 0, "0"]:
			self.astObject.set("Closed=0")
		else:
			raise ValueError("ASTRegion.isClosed property must be one of [True, False, 1, 0].")

	@property
	def isBounded(self) -> bool:
		''' Boolean attribute that indicates whether the region is bounded. '''
		return self.astObject.get("Bounded") == "1"

	@isBounded.setter
	def isBounded(self, newValue:bool):
		if newValue in [True, 1, "1"]:
			self.astObject.set("Bounded=1")
		elif newValue in [False, 0, "0"]:
			self.astObject.set("Bounded=0")
		else:
			raise ValueError("ASTRegion.isBounded property must be one of [True, False, 1, 0].")

	def frame(self) -> ASTFrame:
		'''
		Returns a copy of the frame encapsulated by this region.

		Note that the frame is not directly accessible; both this method and the underlying ``starlink-pyast`` function returns a copy.
		'''
		# this is not a property since a new object is being returned.
		ast_frame = self.astObject.getregionframe() # "A pointer to a deep copy of the Frame represented by the Region."
		return ASTFrame.frameFromAstObject(ast_frame)

	def frameSet(self) -> ASTFrameSet:
		'''
		Returns a copy of the frameset encapsulated by this region.

		From AST docs:

		::

		  The base Frame is the Frame in which the box was originally
		  defined, and the current Frame is the Frame into which the
		  Region is currently mapped (i.e. it will be the same as the
		  Frame returned by astGetRegionFrame).

		'''
		raise NotImplementedError("getregionframeset() has not yet been exposed to the Python interface.")
		return ASTFrameSet(ast_object=self.astObject.getregionframeset())

	@property
	def meshSize(self) -> int:
		''' Number of points used to create a mesh covering the region. '''
		#return int(self.astObject.get("MeshSize"))
		return int(self.astObject.MeshSize)

	@meshSize.setter
	def meshSize(self, newValue:int):
		newValue = as_integer(newValue, "meshSize")
		if newValue < 5:
			# raise rather than silently clamp: a rewritten value is a hidden surprise
			raise ValueError(f"'meshSize' must be at least 5 (got {newValue}).")
		if newValue > MAX_MESH_SIZE:
			raise ValueError(f"'meshSize' must be at most {MAX_MESH_SIZE} (got {newValue}); "
			                 f"see MAX_MESH_SIZE for the physical rationale.")
		self.astObject.set("MeshSize={0}".format(newValue))


	@property
	def fillFactor(self):
		''' <Fraction of the Region which is of interest> '''
		# TODO: properly document, see p. 812 of documentation
		return self.astObject.get("FillFactor")

	@fillFactor.setter
	def fillFactor(self, newValue):
		raise NotImplementedError("Setting 'fillFactor' is not yet implemented.")

	@property
	def bounds(self) -> Tuple:
		'''
		Returns lower and upper coordinate points that bound this region.
		'''

		lower_bounds, upper_bounds = self.astObject.getregionbounds()

		# lower_bounds and upper_bounds are n-dimensional arrays
		# e.g. for a 2D image,
		# [-10,5], [10,20] <- (ra, dec) or pixel bounds

		lower_bounds = bridge.from_frame_units(lower_bounds, self.astObject)
		upper_bounds = bridge.from_frame_units(upper_bounds, self.astObject)

		return (lower_bounds, upper_bounds)

	def boundingBox(self) -> ASTBox:
		'''
		Returns an ASTBox region that bounds this region where the box sides align with RA, dec.
		'''
		from cornish import ASTBox # import here to avoid circular import
		return ASTBox.fromCorners(frame=self.frame(), corners=self.bounds)
		#raise NotImplementedError()
		# use the "bounds" method above to create a bounding box

	def boundingCircle(self) -> ASTCircle:
		'''
		Returns the smallest circle (:py:class:`ASTCircle`) that bounds this region.

		It is up to the caller to know that this is a 2D region (only minimal checks are made).
		:raises cornish.exc.NotA2DRegion: raised when attempting to get a bounding circle for a region that is not 2D
		'''

		if self.naxes != 2:
			raise NotA2DRegion(f"A bounding circle can only be computed on a 2D region; this region has {self.naxes} axes.")

		from .circle import ASTCircle
		centre, radius = self.astObject.getregiondisc() # returns frame units (radians on sky)
		# note that AST signals "no disc" with Ast.BAD (a large *finite* float), not NaN;
		# this guard stays ahead of the bridge call for its domain-specific message
		if (Ast.BAD in centre) or (radius == Ast.BAD) or \
		   not (np.all(np.isfinite(centre)) and math.isfinite(radius)):
			# e.g. an empty region (such as a MOC of the AND of disjoint regions);
			# fail here rather than let bad values propagate
			raise ValueError(f"No bounding circle could be determined for this region (AST returned centre={centre}, radius={radius}); the region may be empty or unbounded.")
		return ASTCircle(frame=self,
		                 center=bridge.from_frame_units(centre, self.astObject),
		                 radius=bridge.from_frame_distance(radius, self.astObject))

	def overlapType(self, region:Union[ASTRegion, starlink.Ast.Region]) -> OverlapType:
		'''
		Return the nature of the overlap between this region and the provided one
		as an :class:`~cornish.enums.OverlapType` (the ``astOverlap`` return code).

		Unlike the boolean predicates below, this includes
		:attr:`OverlapType.NO_FRAME_MAPPING` rather than raising.

		:param region: the region to compare against
		'''
		if region is None:
			raise ValueError("'None' was provided as the second region.")
		if isinstance(region, ASTRegion):
			return_value = self.astObject.overlap(region.astObject)
		elif isinstance(region, starlink.Ast.Region):
			return_value = self.astObject.overlap(region)
		else:
			raise TypeError(f"Unexpected object given for region; expected either ASTRegion or starlink.Ast.Region (got '{type(region)}').")
		return OverlapType(return_value)

	def _checkedOverlapType(self, region) -> OverlapType:
		''' overlapType, raising when the two regions' coordinate systems cannot be mapped. '''
		overlap = self.overlapType(region)
		if overlap == OverlapType.NO_FRAME_MAPPING:
			raise CoordinateSystemsCouldNotBeMapped("The provided region's coordinate system could not be mapped to this region's system.")
		return overlap

	def overlaps(self, region) -> bool:
		'''
		Return ``True`` if this region overlaps with the provided one.
		'''
		return self._checkedOverlapType(region) in \
			(OverlapType.INSIDE, OverlapType.CONTAINS, OverlapType.PARTIAL, OverlapType.IDENTICAL)

	def isIdenticalTo(self, region:ASTRegion) -> bool:
		'''
		Returns 'True' if this region is identical (to within their uncertainties) to the provided region, 'False' otherwise.
		'''
		return self._checkedOverlapType(region) == OverlapType.IDENTICAL

	def isFullyWithin(self, region:ASTRegion) -> bool:
		'''
		Returns 'True' if this region is fully within the provided region.
		'''
		return self._checkedOverlapType(region) == OverlapType.INSIDE

	def fullyEncloses(self, region:ASTRegion) -> bool:
		'''
		Returns 'True' if this region fully encloses the provided region.
		'''
		return self._checkedOverlapType(region) == OverlapType.CONTAINS

	def isNegationOf(self, region) -> bool:
		'''
		Returns 'True' if this region is the exact negation of the provided region.
		'''
		return self._checkedOverlapType(region) == OverlapType.NEGATION

	def maskOnto(self, image=None, mapping=None, fits_coordinates:bool=True, lower_bounds=None, mask_inside=True, mask_value=float("NaN")):
		'''
		Apply this region as a mask on top of the provided image; note: the image values are overwritten!

		:param image: numpy.ndarray of pixel values (or other array of values)
		:param mapping: mapping from this region to the pixel coordinates of the provided image
		:param fits_coordinates: use the pixel coordinates of a FITS file (i.e. origin = [0.5, 0.5] for 2D)
		:param lower_bounds: lower bounds of provided image, only specify if not using FITS coordinates
		:param mask_inside: True: mask the inside of this region; False: mask outside of this region
		:param mask_value: the value to set the masked image pixels to
		:returns: number of pixels in image masked
		'''

		# coded now for numpy arrays, but set ndim,shape for anything else
		ndim = len(image.shape)

		# assert number of axes in image == # of outputs in the mapping
		if ndim != mapping.numberOfOutputCoordinates:
			raise ValueError(f"The number of dimensions in the provided image ({ndim}) does not match the number of output dimensions of the provided mapping ({mapping.numberOfOutputCoordinates}).")

		if fits_coordinates:
			# use the pixel coordinates for FITS files -> origin at [0.5, 0.5]
			lower_bounds = [0.5 for x in range(ndim)]
		else:
			# must use provided lower bounds
			if lower_bounds is None:
				raise ValueError("'lower_bounds' must be provided (or specify 'fits_coordinates=True' to use FITS coordinates.")
#		upper_bounds = list()
#		for idx, n in enumerate(shape):
#			upper_bounds.append(lower_bounds[idx] + n)

		npix_masked = self.astObject.mask(mapping.astObject, mask_inside, lower_bounds, image, mask_value)
		return npix_masked

	def regionWithMapping(self, map=None, frame=None) -> ASTRegion:
		'''
		Returns a new ASTRegion with the coordinate system from the supplied frame.

		Corresponds to the ``astMapRegion`` C function (``starlink.Ast.mapregion``).

		:param map: A mapping that can convert coordinates from the system of the current region to that of the supplied frame.
		:param frame: A frame containing the coordinate system for the new region.
		:returns: new ASTRegion with a new coordinate system
		'''
		if frame is None:
			raise ValueError("A frame must be specified.")
		if map is None:
			map = frame # attempt to use the frame as a mapper (an ASTFrame is a subclass of ASTMapper)

		# Would be nice to be able to create an instance of the same subclass of ASTRegion
		# - how to inspect the object for this information?

		if isinstance(map, starlink.Ast.Mapping):
			ast_map = map
		elif isinstance(map, (ASTMapping, ASTFrameSet)): # frame sets contain mappings
			ast_map = map.astObject
		else:
			raise TypeError("Expected 'map' to be one of these two types: starlink.Ast.Mapping, ASTMapping.")

		if isinstance(frame, starlink.Ast.Frame):
			ast_frame = frame
		elif isinstance(frame, (ASTFrame, ASTFrameSet)):
			ast_frame = frame.astObject
		else:
			raise TypeError("Expected 'frame' to be one of these two types: starlink.Ast.Frame, ASTFrame.")

		new_ast_region = self.astObject.mapregion(ast_map, ast_frame)

		# wrap in the appropriate cornish class
		if isinstance(new_ast_region, starlink.Ast.Polygon):
			new_region = cornish.region.ASTPolygon(ast_object=new_ast_region)
		elif isinstance(new_ast_region, starlink.Ast.Box):
			new_region = cornish.region.ASTBox(ast_object=new_ast_region)
		elif isinstance(new_ast_region, starlink.Ast.Circle):
			new_region = cornish.region.ASTCircle(ast_object=new_ast_region)
		else:
			from ..exc import UnsupportedASTClass
			raise UnsupportedASTClass(f"ASTRegion.regionWithMapping: unhandled region type '{new_ast_region.Class}'.")
		new_region.wcs = self.wcs # propagate the originating WCS, when known (SPEC-04 §2)
		return new_region

	def mapRegionMesh(self, mapping=None, frame=None):
		'''
		Returns a new ASTRegion that is the same as this one but with the specified coordinate system.

		Parameters
		----------
		mapping : `~cornish.mapping.ASTMapping` class
			The mapping to transform positions from the current ASTRegion to those specified by the given frame.
		frame : `~cornish.frame.ASTFrame` class
			Coordinate system frame to convert the current ASTRegion to.

		Returns
		-------
		region : `ASTRegion`
			A new region object covering the same area but in the frame specified in `frame`.

		Raises
		------
		Exception
			An exception is raised for missing parameters.
		'''
		if mapping is None or frame is None:
			raise ValueError("A mapping and frame is required to be passed to 'mapRegionMesh'.")

		# check it's the correct type
		if not isinstance(mapping, ASTMapping):
			raise TypeError("The object passed to 'mapping' needs to be an ASTMapping object.")

		if not isinstance(frame, ASTFrame):
			raise TypeError("The object passed to 'frame' needs to be an ASTFrame object.")

		self.astObject.mapregionmesh( mapping, frame )

	def boundaryPointMesh(self, npoints:int=None) -> np.ndarray:
		'''
		Returns an array of evenly distributed points that cover the boundary of the region.
		For example, if the region is a box, it will generate a list of points that trace the edges of the box.

		The default value of 'npoints' is 200 for 2D regions and 2000 for three or more dimensions.

		:param npoints: the approximate number of points to generate in the mesh
		:returns: list of points in degrees
		'''
		# The starlink.AST object uses the attribute "MeshSize" to determine the number of points to
		# use. This should be specified when building the mesh - the attribute doesn't seem to be used
		# anywhere else. This method shouldn't change the value in case that's not true, but we'll make
		# this one step here.

		#if npoints is not None and not isinstance(npoints, int):
		#	raise Exception("The parameter 'npoints' must be an integer ('{1}' provided).".format(npoints))

		if npoints is None:
			pass # use default meshSize value
		else:
			# use provided value
			#old_mesh_size = self.astObject.get("MeshSize")
			#self.astObject.set("MeshSize={0}".format(npoints))
			old_mesh_size = self.meshSize
			self.meshSize = npoints

		try:
			mesh = self.astObject.getregionmesh(MeshType.BOUNDARY)
			# As of starlink-pyast v.3.15.4, getregionmesh returns normalized points;
			# from_frame_units norms again, which is a harmless idempotent call (M14).
		except Ast.MBBNF as e:
			raise ValueError("A mesh could not be generated: no bounding box could be found for this region's mapping.") from e
		finally:
			if npoints is not None:
				# restore original value even when the mesh generation fails
				self.meshSize = old_mesh_size

		return bridge.from_frame_units(mesh, self.astObject) # a list of point pairs, not two parallel arrays

	def interiorPointMesh(self, npoints:int=None):
		'''
		Returns an array of evenly distributed points that cover the surface of the region.
		For example, if the region is a box, it will generate a list of points that fill the interior area of the box.

		The default value of 'npoints' is 200 for 2D regions and 2000 for three or more dimensions.

		:param npoints: the approximate number of points to generate in the mesh
		:returns: array of points in degrees
		'''
		# See discussion of "MeshSize" in method "boundaryPointMesh".

		if npoints is None:
			pass # use default value
		else:
			# route through the property so its validation (type, bounds)
			# applies identically to boundaryPointMesh and interiorPointMesh
			old_mesh_size = self.meshSize
			self.meshSize = npoints

		# The returned "points" array from getregionmesh() will be a 2-dimensional array with shape (ncoord,npoint),
		# where "ncoord" is the number of axes within the Frame represented by the Region,
		# and "npoint" is the number of points in the mesh (see attribute "MeshSize").
		try:
			mesh = self.astObject.getregionmesh(MeshType.SURFACE) # AST's "surface" mesh is the interior
			# As of starlink-pyast v.3.15.4, getregionmesh returns normalized points;
			# from_frame_units norms again, which is a harmless idempotent call (M14).
		finally:
			if npoints is not None:
				# restore original value even when the mesh generation fails
				self.meshSize = old_mesh_size

		return bridge.from_frame_units(mesh, self.astObject)

	def containsPoint(self, point:Union[Iterable, astropy.coordinates.SkyCoord]=None) -> bool:
		'''
		Returns ``True`` if the provided point lies inside this region, ``False`` otherwise.

		This method is a direct synonym for :meth:`pointInRegion`.
		The name "containsPoint" is more appropriate from an object perspective,
		but the ``pointInRegion`` method is present for consistency with the AST library.

		:param point: a coordinate point in the same frame as this region
		'''
		return self.pointInRegion(point=point)

	def pointInRegion(self, point:Union[Iterable, astropy.coordinates.SkyCoord, np.ndarray]) -> bool:
		'''
		Returns ``True`` if the provided point lies inside this region, ``False`` otherwise.

		The point may be a SkyCoord (any system — converted), a Quantity pair, or
		bare values: degrees when this region is in a sky frame, the frame's
		native units otherwise (e.g. pixels for a pixel-frame region).

		:param point: a coordinate point in the same frame as this region
		'''
		if point is None:
			raise ValueError("A test point must be specified.")

		return self.astObject.pointinregion(bridge.to_frame_units(point, self.astObject, squeeze=True))

	@abstractproperty
	def area(self) -> astropy.units.quantity.Quantity:
		# subclasses must implement
		raise NotImplementedError()

	@abstractmethod
	def toPolygon(self, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec) -> ASTPolygon:
		'''
		Method that guarantees returning a polygon that describes or approximates this region.

		This method provides a common interface to create polygons from different region types.
		Calling this on an ASTPolygon will return itself; calling it on an ASTCircle
		will return a polygon that approximates the circle. The parameters 'npoints' and
		'maxerr' will be used only when appropriate.

		:param npoints: number of points to sample from the region's boundary for the resulting polygon
		:param maxerr:
		'''
		pass

