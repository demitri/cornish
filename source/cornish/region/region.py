
from __future__ import annotations # remove in Python 3.10
# Needed for forward references, see:
# https://stackoverflow.com/a/33533514/2712652

from abc import ABCMeta, abstractproperty
from typing import Union, Iterable

import math
from math import radians as deg2rad
from math import degrees as rad2deg

import numpy as np
import astropy
import astropy.units as u
import starlink
import starlink.Ast as Ast

import cornish.region # to avoid circular imports below - better way?
from ..mapping import ASTFrame, ASTFrameSet, ASTMapping
from ..exc import NotA2DRegion, CoordinateSystemsCouldNotBeMapped

__all__ = ['ASTRegion']

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
  * astGetRegionMesh: Get a mesh of points covering a Region
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
	This is an abstract superclass - there is no means to create an ASTRegion object directly
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
	
	def __add__(self, other_region):
		# TODO: check data type, handle both ASTRegion and the ast_object region?
		from .compound_region import ASTCompoundRegion # forward import to avoid circular import error
		return ASTCompoundRegion(regions=[self, other_region], operation=Ast.AND)
	
	@classmethod
	def fromFITSHeader(cls, fits_header=None, uncertainty:float=4.848e-6):
		'''
		Factory method to create a region from the provided FITS header; the returned object will be as specific as possible (but probably an ASTPolygon).
		
		The frame is determined from the FITS header.
		
		:param fits_header: a FITS header (Astropy, fitsio, an array of cards)
		:param uncertainty: defaults to 4.848e-6, an uncertainty of 1 arcsec
		'''
		if fits_header is None:
			raise ValueError("This method requires a 'fits_header' to be set.")
		
		# import here to avoid circular imports
		from .box import ASTBox
		from .circle import ASTCircle
		from .polygon import ASTPolygon
		from ..channel import ASTFITSChannel
		
		# get frame from header
		fits_channel = ASTFITSChannel(header=fits_header)
		
		# does this channel contain a frame set?
		frame_set = fits_channel.frameSet
		if frame_set is None:
			raise ValueError("The provided FITS header does not describe a region (e.g. not an image, does not contain a WCS that AST can read).")
		
		frame = frame_set.baseFrame
		
		# support n-dimensional regions
		
		# define needed parameters for region creation below
		dimensions = fits_channel.dimensions
		n_dim = len(dimensions)
		cornerPoint = [0.5 for x in range(n_dim)]
		cornerPoint2 = [dimensions[x] + 0.5 for x in range(n_dim)]
		#cornerPoint=[0.5,0.5], # center of lower left pixel
		#cornerPoint2=[dimensions[0]+0.5, dimensions[1]+0.5])

		if n_dim > 2:
			raise NotImplementedError("HDUs describing dimensions greater than 2 not currently supported.")
			
		#if isinstance(frame, ASTFrame):
		#	self.frame = frame
		#elif isinstance(frame, starlink.Ast.Frame):
		#	self.frame = ASTFrame(frame=frame)
		#else:
		#	raise Exception("ASTBox: unexpected frame type specified ('{0}').".format(type(frame)))

		#if all([cornerPoint,centerPoint]) or all([cornerPoint,cornerPoint2]) or dimensions is not None:
		#	if dimensions is not None:
		#		input_form = CORNER_CORNER
		#		p1 = [0.5,0.5] # use 0.5 to specify the center of each pixel
		#		p2 = [dimensions[0]+0.5, dimensions[1]+0.5]
		#	elif centerPoint is None:
		#		input_form = CORNER_CORNER
		#		p1 = [cornerPoint[0], cornerPoint[1]]
		#		p2 = [cornerPoint2[0], cornerPoint2[1]]
		#		dimensions = [math.fabs(cornerPoint[0] - cornerPoint2[0]),
		#					  math.fabs(cornerPoint[1] - cornerPoint2[1])]
		#	else:
		#		input_form = CENTER_CORNER
		#		p1 = [centerPoint[0], centerPoint[1]]
		#		p2 = [cornerPoint[0], cornerPoint[1]]
		#		dimensions = [2.0 * math.fabs(centerPoint[0] - cornerPoint[0]),
		#					  2.0 * math.fabs(centerPoint[1] - cornerPoint[1])]

		# input_form constants (define properly elsewhere?)
		CENTER_CORNER = 0
		CORNER_CORNER = 1

		input_form = CORNER_CORNER
		p1 = [cornerPoint[0], cornerPoint[1]]
		p2 = [cornerPoint2[0], cornerPoint2[1]]
		dimensions = [math.fabs(cornerPoint[0] - cornerPoint2[0]),
					  math.fabs(cornerPoint[1] - cornerPoint2[1])]


		#dimensions = [dimensions[0], dimensions[1]]
		#logger.debug("Setting dims: {0}".format(self.dimensions))

		ast_object = Ast.Box( frame.astObject, input_form, p1, p2, unc=uncertainty )
		
		# create the mapping from pixel to sky (or whatever is there) if available
		mapping = frame_set.astObject.getmapping() # defaults are good
		current_frame = frame_set.astObject.getframe(starlink.Ast.CURRENT)
		
		# create a new region with new mapping
		ast_object = ast_object.mapregion(mapping, current_frame)
		
		if isinstance(ast_object, Ast.Box):
			from .box import ASTBox # avoid circular imports
			return ASTBox(ast_object=ast_object)
		elif isinstance(ast_object, Ast.Circle):
			from .circle import ASTCircle # avoid circular imports
			return ASTCircle(ast_object=ast_object)
		elif isinstance(ast_object, Ast.Polygon):
			return ASTPolygon(ast_object=ast_object)
		else:
			raise Exception(f"Unexpected region type encountered: {type(ast_object)}.")

	@property
	def points(self) -> np.ndarray:
		'''
		The array of points that define the region. The interpretation of the points is dependent on the type of shape in question.
		
		Box: returns two points; the center and a box corner.
		Circle: returns two points; the center and a point on the circumference.
		CmpRegion: no points returned; to get points that define a compound region, call this method on each of the component regions via the method "decompose".
		Ellipse: three points: 1) center, 2) end of one axis, 3) end of the other axis
		Interval: two points: 1) lower bounds position, 2) upper bounds position (reversed when interval is an excluded interval)
		NullRegion: no points
		PointList: positions that the list was constructed with
		Polygon: vertex positions used to construct the polygon
		Prism: no points (see CmpRegion)

		NOTE: points returned reflect the current coordinate system and may be different from the initial construction
		
		:returns: NumPy array of coordinate points in degrees, shape (ncoord,2), e.g. [[ra1,dec1], [ra2, dec2], ..., [ra_n, dec_n]]
		'''
		
		# getregionpoints returns data as [[x1, x2, ..., xn], [y1, y2, ..., yn]]
		# transpose the points before returning
		# also, normalize points to expected bounds
		return np.rad2deg(self.astObject.norm(self.astObject.getregionpoints())).T

	@property
	def uncertainty(self):
		'''
		Uncertainties associated with the boundary of the Box.
					
		The uncertainty in any point on the boundary of the Box is found by
		shifting the supplied "uncertainty" Region so that it is centered at
		the boundary point being considered. The area covered by the shifted
		uncertainty Region then represents the uncertainty in the boundary
		position. The uncertainty is assumed to be the same for all points.
		'''
		return self._uncertainty
			
	@uncertainty.setter
	def uncertainty(self, unc):
		raise Exception("Setting the uncertainty currently doesn't do anything.")
		self._uncertainty = unc

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
			raise Exception("ASTRegion.adaptive property must be one of [True, False, 1, 0].")
	
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
			raise Exception("ASTRegion.isNegated property must be one of [True, False, 1, 0].")
	
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
			raise Exception("ASTRegion.isClosed property must be one of [True, False, 1, 0].")
	
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
			raise Exception("ASTRegion.isBounded property must be one of [True, False, 1, 0].")
	
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
		if isinstance(newValue, int):
			if newValue < 5:
				newValue = 5
			self.astObject.set("MeshSize={0}".format(newValue))
		else:
			raise Exception("ASTRegion.meshSize: an integer value of at least 5 is required.")
	

	@property
	def fillFactor(self):
		''' <Fraction of the Region which is of interest> '''
		# TODO: properly document, see p. 812 of documentation
		return self.astObject.get("FillFactor")
	
	@fillFactor.setter
	def fillFactor(self, newValue):
		raise Exception("TODO: document and implement")

	def bounds(self):
		'''
		Return the upper and lower bounds of a box that contains this regions.
		'''
		
		lower_bounds, upper_bounds = self.astObject.getregionbounds()
		
		# lower_bounds and upper_bounds are n-dimensional arrays
		# e.g. for a 2D image,
		# [-10,5], [10,20] <- ra, dec or pixel bounds
		
		raise Exception("test units?")
		# .. todo:: if SkyFrame convert from ra -> deg
		
		return (lower_bounds, upper_bounds)
		
	def boundingBox(self):
		'''
		
		'''
		raise NotImplementedError()
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
		centre, radius = self.astObject.getregiondisc() # returns radians
		return ASTCircle(frame=self, center=np.rad2deg(centre), radius=rad2deg(radius))

	def overlaps(self, region) -> bool:
		'''
		Return ``True`` if this region overlaps with the provided one.
		'''
		if region is None:
			raise ValueError("'None' was provided as the second region.")
		if isinstance(region, ASTRegion):
			return_value = self.astObject.overlap(region.astObject)
		elif isinstance(region, starlink.Ast.Region):
			return_value = self.astObject.overlap(region)
		else:
			raise ValueError(f"Unexpected object given for region; expected either ASTRegion or starlink.Ast.Region (got '{type(region)}').")
			
		if return_value == 0:
			raise CoordinateSystemsCouldNotBeMapped("The provided region's coordinate system could not be mapped to this region's system.")
		elif return_value == 1:
			return False 			# no overlap
		elif return_value == 2:
			return True 			# this region is completely inside the provded region
		elif return_value == 3:
			return True 			# the provded region is completely inside the first region
		elif return_value == 4:
			return True 			# there is partial overlap
		elif return_value == 5:
			return True				# the resions are identical to within their uncertainties
		elif return_value == 6:
			return False			# the second region is the exact negation of this region to within their uncertainties

	def isIdenticalTo(self, region:ASTRegion) -> bool:
		'''
		Returns 'True' if this region is identical (to within their uncertainties) to the provided region, 'False' otherwise.
		'''
		if region is None:
			raise ValueError("'None' was provided as the second region.")
		if isinstance(region, ASTRegion):
			return_value = self.astObject.overlap(region.astObject)
		elif isinstance(region, starlink.Ast.Region):
			return_value = self.astObject.overlap(region)
		else:
			raise ValueError(f"Unexpected object given for region; expected either ASTRegion or starlink.Ast.Region (got '{type(region)}').")

		if return_value == 0:
			raise CoordinateSystemsCouldNotBeMapped("The provided region's coordinate system could not be mapped to this region's system.")
		else:
			return return_value == 5

	def isFullyWithin(self, region:ASTRegion) -> bool:
		'''
		Returns 'True' if this region is fully within the provided region.
		'''
		if region is None:
			raise ValueError("'None' was provided as the second region.")
		if isinstance(region, ASTRegion):
			return_value = self.astObject.overlap(region.astObject)
		elif isinstance(region, starlink.Ast.Region):
			return_value = self.astObject.overlap(region)
		else:
			raise ValueError(f"Unexpected object given for region; expected either ASTRegion or starlink.Ast.Region (got '{type(region)}').")

		if return_value == 0:
			raise CoordinateSystemsCouldNotBeMapped("The provided region's coordinate system could not be mapped to this region's system.")
		else:
			return return_value == 2

	def fullyEncloses(self, region:ASTRegion) -> bool:
		'''
		Returns 'True' if this region fully encloses the provided region.
		'''
		if region is None:
			raise ValueError("'None' was provided as the second region.")
		if isinstance(region, ASTRegion):
			return_value = self.astObject.overlap(region.astObject)
		elif isinstance(region, starlink.Ast.Region):
			return_value = self.astObject.overlap(region)
		else:
			raise ValueError(f"Unexpected object given for region; expected either ASTRegion or starlink.Ast.Region (got '{type(region)}').")

		if return_value == 0:
			raise CoordinateSystemsCouldNotBeMapped("The provided region's coordinate system could not be mapped to this region's system.")
		else:
			return return_value == 3

	def isNegationOf(self, region):
		'''
		Returns 'True' if this region is the exact negation of the provided region.
		'''
		if region is None:
			raise ValueError("'None' was provided as the second region.")
		if isinstance(region, ASTRegion):
			return_value = self.astObject.overlap(region.astObject)
		elif isinstance(region, starlink.Ast.Region):
			return_value = self.astObject.overlap(region)
		else:
			raise ValueError(f"Unexpected object given for region; expected either ASTRegion or starlink.Ast.Region (got '{type(region)}').")

		if return_value == 0:
			raise CoordinateSystemsCouldNotBeMapped("The provided region's coordinate system could not be mapped to this region's system.")
		else:
			return return_value == 6

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
		shape = image.shape # <-- unused variable!
		
		# assert number of axes in image == # of outputs in the mapping
		if ndim != mapping.number_of_output_coordinates:
			raise Exception(f"The number of dimensions in the provided image ({ndim}) does not match the number of output dimensions of the provided mapping ({mapping.number_of_output_coordinates}).")
		
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
		
		:param frame: A frame containing the coordinate system for the new region.
		:param map: A mapping that can convert coordinates from the system of the current region to that of the supplied frame.
		:returns: new ASTRegion with a new coordinate system
		'''
		if frame is None:
			raise Exception("A frame must be specified.")
		if map is None:
			map = frame # attempt to use the frame as a mapper (an ASTFrame is a subclass of ASTMapper)
		
		# Would be nice to be able to create an instance of the same subclass of ASTRegion
		# - how to inspect the object for this information?
		
		if isinstance(map, starlink.Ast.Mapping):
			ast_map = map
		elif isinstance(map, (ASTMapping, ASTFrameSet)): # frame sets contain mappings
			ast_map = map.astObject
		else:
			raise Exception("Expected 'map' to be one of these two types: starlink.Ast.Mapping, ASTMap.")

		if isinstance(frame, starlink.Ast.Frame):
			ast_frame = frame
		elif isinstance(map, (ASTFrame, ASTFrameSet)):
			ast_frame = frame.astObject
		else:
			raise Exception("Expected 'frame' to be one of these two types: starlink.Ast.Frame, ASTFrame.")
		
		new_ast_region = self.astObject.mapregion(ast_map, ast_frame)
		
		# This is temporary and probably fragile. Find a replacement for this ASAP.
		# get the returned region type to create the correct wrapper
		#
		# -> make a deep copy, replace obj.astObject with new one (check any properties)
		#
		if new_ast_region.Class == 'Polygon' or isinstance(new_ast_region, starlink.Ast.Polygon):
			return cornish.region.ASTPolygon(ast_object=new_ast_region)
		elif new_ast_region.Class == 'Box' or isinstance(new_ast_region, starlink.Ast.Box):
			return cornish.region.ASTBox(ast_object=new_ast_region)
		else:
			raise Exception("ASTRegion.regionWithMapping: unhandled region type (easy fix).")
	
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
			raise Exception("A mapping and frame is required to be passed to 'mapRegion'.")

		# check it's the correct type
		if not isinstance(mapping, ASTMapping):
			raise Exception("The object passed to 'mapping' needs to be an ASTMapping object.")
		
		if not isinstance(frame, ASTFrame):
			raise Exception("The object passed to 'frame' needs to be an ASTFrame object.")
		
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
			mesh = self.astObject.getregionmesh(1) # surface=1, here "surface" means the boundary
			# if a basic frame is used instead of a sky frame, the points need to be normalized on [0,360)
			mesh = self.astObject.norm(mesh)
		except Ast.MBBNF as e:
			print("AST error: Mapping bounding box not found. ({0})".format(e))
			raise e
		
		if npoints is not None:
			# restore original value
			self.meshSize = old_mesh_size #self.astObject.set("MeshSize={0}".format(old_mesh_size))
		
		return np.rad2deg(mesh).T # returns as a list of pairs of points, not two parallel arrays
		
	def interiorPointMesh(self, npoints:int=None):
		'''
		Returns an array of evenly distributed points that cover the surface of the region.
		For example, if the region is a box, it will generate a list of points that fill the interior area of the box.
		
		The default value of 'npoints' is 200 for 2D regions and 2000 for three or more dimensions.
		
		:param npoints: The approximate number of points to generate in the mesh.
		:returns: array of points in degrees
		'''
		# See discussion of "MeshSize" in method "boundaryPointMesh".
		
		if npoints is not None and not isinstance(npoints, int):
			raise Exception(f"The parameter 'npoints' must be an integer ('{type(npoints)}' provided).")

		if npoints is None:
			pass # use default value
		else:
			# use provided value
			old_mesh_size = self.astObject.get("MeshSize")
			self.astObject.set("MeshSize={0}".format(npoints))

		# The returned "points" array from getregionmesh() will be a 2-dimensional array with shape (ncoord,npoint),
		# where "ncoord" is the number of axes within the Frame represented by the Region,
		# and "npoint" is the number of points in the mesh (see attribute "MeshSize").
		mesh = self.astObject.getregionmesh(0) # surface=0, here "surface" means the interior
		mesh = self.astObject.norm(mesh)
		
		if npoints is not None:
			# restore original value
			self.astObject.set("MeshSize={0}".format(old_mesh_size))

		# .. todo:: double check points as degrees vs radians!

		return np.rad2deg(mesh).T
	
	def containsPoint(self, point:Union[Iterable, astropy.coordinates.SkyCoord]=None) -> bool:
		'''
		Returns ``True`` if the provided point lies inside this region, ``False`` otherwise.

		This method is a direct synonym for :meth:`pointInRegion`.
		The name "containsPoint" is more appropriate for the object oriented format,
		but the ``pointInRegion`` method is present for consistency with the AST library.
		'''
		return self.pointInRegion(point=point)
	
	def pointInRegion(self, point:Union[Iterable, astropy.coordinates.SkyCoord,np.ndarray]=None) -> bool:
		'''
		Returns ``True`` if the provided point lies inside this region, ``False`` otherwise.
		
		If no units are specified degrees are assumed.
		'''
		if point is None:
			raise ValueError("A test point must be specified.")
		
		if isinstance(point, astropy.coordinates.SkyCoord):
			point = [point.ra.to(u.rad).value, point.dec.to(u.rad).value]
		else:
			point = np.deg2rad(point)
		
		return self.astObject.pointinregion(point)
		
	@abstractproperty
	def area(self) -> astropy.units.quantity.Quantity:
		# subclasses  must implement
		raise NotImplementedError()



