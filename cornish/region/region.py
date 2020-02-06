
from abc import ABCMeta
from typing import Union

from math import radians as deg2rad
from math import degrees as rad2deg

import numpy as np
import astropy
import astropy.units as u
import starlink
import starlink.Ast as Ast

import cornish.region # to avoid circular imports below - better way?
from ..ast_object import ASTObject
from ..mapping import ASTFrame, ASTFrameSet, ASTMapping
from ..exc import NotA2DRegion

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
	(see ASTBox, ASTPolygon, etc.).
	
	self.astObject is of type starlink.Ast.Region.
	'''
	
	def __init__(self, ast_object=None, uncertainty=None):
	  '''
	  
	  '''
	  super().__init__(ast_object=ast_object)
	  self._uncertainty = None
		
	@property
	def points(self, units:astropy.units.core.Unit=u.deg):
		'''
		The array of points that define the region.
		
		:param unit: the unit of the points requested (astroy.units.deg or astropy.units.rad)
		:returns: Numpy array of coordinate points.
		'''
		
		# getregionpoints returns data as [[x1, x2, ..., xn], [y1, y2, ..., yn]]
		# transpose the points before returning
		points = self.astObject.norm(self.astObject.getregionpoints())
		if units == u.deg:
			return np.rad2deg(points.T)
		elif units == u.rad:
			return points.T
		else:
			raise ValueError("Requested units for points must be either astropy.units.deg or astropy.units.rad.")
		
# 		if units == u.deg:
# 			return np.rad2deg(self.astObject.getregionpoints().T)
# 		elif units == u.rad:
# 			return self.astObject.getregionpoints().T
# 		else:
# 			raise ValueError("Requested units for points must be either astropy.units.deg or astropy.units.rad.")

	@property
	def uncertainty(self):
		'''
		Uncertainties associated with the boundary of the Box.
					
		The uncertainty in any point on the boundary of the Box is found by
		shifting the supplied "uncertainty" Region so that it is centered at
		the boundary point being considered.  The area covered by the shifted
		uncertainty Region then represents the uncertainty in the boundary
		position.  The uncertainty is assumed to be the same for all points.
		'''
		return self._uncertainty
			
	@uncertainty.setter
	def uncertainty(self, unc):
		raise Exception("Setting the uncertainty currently doesn't do anything.")
		self._uncertainty = unc

	@property
	def isAdaptive(self):
		''' Boolean attribute that indicates whether the area adapt to changes in the coordinate system. '''
		return self.astObject.get("Adaptive")
	
	@isAdaptive.setter
	def isAdaptive(self, newValue):
		if newValue in [True, 1]:
			self.astObject.set("Adaptive=1")
		elif newValue in [False, 0]:
			self.astObject.set("Adaptive=0")
		else:
			raise Exception("ASTRegion.adaptive property must be one of [True, False, 1, 0].")
	
	@property
	def isNegated(self, newValue):
		''' Boolean attribute that indicates whether the original region has been negated. '''
		return self.astObject.get("Negated")
	
	@isNegated.setter
	def isNegated(self, newValue):
		if newValue in [True, 1]:
			self.astObject.set("Negated=1")
		elif newValue in [False, 0]:
			self.astObject.set("Negated=0")
		else:
			raise Exception("ASTRegion.isNegated property must be one of [True, False, 1, 0].")
	
	@property
	def isClosed(self):
		''' Boolean attribute that indicates whether the boundary be considered to be inside the region. '''
		return self.astObject.get("Closed")

	@isClosed.setter
	def isClosed(self, newValue):
		if newValue in [True, 1]:
			self.astObject.set("Closed=1")
		elif newValue in [False, 0]:
			self.astObject.set("Closed=0")
		else:
			raise Exception("ASTRegion.isClosed property must be one of [True, False, 1, 0].")
	
	@property
	def isBounded(self):
		''' Boolean attribute that indicates whether the region is bounded. '''
		return self.astObject.get("Bounded")

	@isBounded.setter
	def isBounded(self, newValue):
		if newValue in [True, 1]:
			self.astObject.set("Bounded=1")
		elif newValue in [False, 0]:
			self.astObject.set("Bounded=0")
		else:
			raise Exception("ASTRegion.isBounded property must be one of [True, False, 1, 0].")
		
	@property
	def meshSize(self):
		''' Number of points used to create a mesh covering the region. '''
		return self.astObject.get("MeshSize")
	
	@meshSize.setter
	def meshSize(self, newValue):
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
		shape = image.shape
		
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

		npix_masked = self.astObject.mask(mapping.astObject, mask_inside, lower_bounds, image, value)
		return npix_masked
		
	def boundingBox(self):
		'''
		
		'''
		pass
		# use the "bounds" method above to create a bounding box
		
	

	def regionWithMapping(self, map=None, frame=None):
		'''
		Returns a new ASTRegion with the coordinate system from the supplied frame.
		
		Corresponds to the astMapRegion C function (starlink.Ast.mapregion).
		
		@param frame A frame containing the coordinate system for the new region.
		@param map A mapping that can convert coordinates from the system of the current region to that of the supplied frame.
		@returns New ASTRegion with a new coordinate system
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
			return cornish.region.ASTBox(ast_box=new_ast_region)
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

	def boundaryPointMesh(self, npoints=None):
		'''
		Returns an array of evenly distributed points that cover the boundary of the region.
		For example, if the region is a box, it will generate a list of points that trace the edges of the box.
		
		The default value of 'npoints' is 200 for 2D regions and 2000 for three or more dimensions.
		
		@param npoints The approximate number of points to generate in the mesh.
		@returns List of points.
		'''
		# The starlink.AST object uses the attribute "MeshSize" to determine the number of points to
		# use. This should be specified when building the mesh - the attribute doesn't seem to be used
		# anywhere else. This method shouldn't change the value in case that's not true, but we'll make
		# this one step here.
		
		if npoints is not None and not isinstance(npoints, int):
			raise Exception("The parameter 'npoints' must be an integer ('{1}' provided).".format(npoints))
		
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
		except Ast.MBBNF as e:
			print("AST error: Mapping bounding box not found. ({0})".format(e))
			raise e
		
		if npoints is not None:
			# restore original value
			self.astObject.meshSize = old_mesh_size #self.astObject.set("MeshSize={0}".format(old_mesh_size))
		
		return mesh.T # returns as a list of pairs of points, not two parallel arrays
		
	def interiorPointMesh(self, npoints=None):
		'''
		Returns an array of evenly distributed points that cover the surface of the region.
		For example, if the region is a box, it will generate a list of points that fill the interior area of the box.
		
		The default value of 'npoints' is 200 for 2D regions and 2000 for three or more dimensions.
		
		:param npoints: The approximate number of points to generate in the mesh.
		:returns: List of points.
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

		return mesh.T
	
	def boundingCircle(self):
		'''
		Returns the smallest circle (ASTCircle) that bounds this region.
		
		It is up to the caller to know that this is a 2D region (only minimal checks are made).
		:raises cornish.exc.NotA2DRegion: raised when attempting to get a bounding circle for a region that is not 2D
		'''
		
		if self.naxes != 2:
			raise NotA2DRegion(f"A bounding circle can only be computed on a 2D region; this region has {self.naxes} axes.")
		
		from .circle import ASTCircle
		centre, radius = self.astObject.getregiondisc() # returns radians
		return ASTCircle(frame=self, center_point=np.rad2deg(centre), radius=rad2deg(radius))
	




