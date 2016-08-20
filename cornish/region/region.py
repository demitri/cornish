from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink
import starlink.Ast as Ast
import astropy.units as u

from ..ast_object import ASTObject
from ..mapping import ASTFrame, ASTFrameSet, ASTMapping
from .polygon import ASTPolygon
from .box import ASTBox

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
		
class ASTRegion(ASTFrame):
	'''
	Represents a region within a coordinate system.
	This is an abstract superclass - there is no means to create an ASTRegion object directly
	(see ASTBox, ASTPolygon, etc.).
	
	self.astObject is of type starlink.Ast.Region.
	'''
	#def __init__(self, ast_frame=None):
	#	'''
	#	
	#	'''
	#	self.astObject = super(ASTRegion, self).__init__(ast_frame=ast_frame)
		
	@property
	def points(self):
		'''
		The array of points that define the region.
		
		@returns Numpy array of coordinate points.
		'''
		
		# getregionpoints returns data as [[x1, x2, ..., xn], [y1, y2, ..., yn]]
		# transpose the points before returning
		return self.astObject.getregionpoints().T

	@property
	def isAdaptive(self):
		''' Boolean attribute that indicates whether the area adapt to changes in the coordinate system. '''
		return self.astObject.get("Adaptive")
	
	@adaptive.setter
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
	
	@adaptive.setter
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
	
	@fillFactor.setter(self, newValue):
		raise Exception("TODO: document and implement")

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
		ast_type = skybox.astObject.__repr__().split("\n")[0].split()[2] # e.g. '< Begin Polygon'
		if ast_type == "Polygon":
			return ASTPolygon(polygon=new_ast_region)
		elif ast_type= "Box":
			return ASTBox()
			
		#return ASTRegion(ast_frame=new_ast_region)
	
	
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
		
		@param npoints The approximate number of points to generate in the mesh.
		@returns List of points.
		'''
		# See discussion of "MeshSize" in method "boundaryPointMesh".
		
		if npoints is not None and not isinstance(npoints, int):
			raise Exception("The parameter 'npoints' must be an integer ('{1}' provided).".format(npoints))

		if npoints is None:
			pass # use default value
		else:
			# use provided value
			old_mesh_size = self.astObject.get("MeshSize")
			self.astObject.set("MeshSize={0}".format(npoints))

		mesh = self.astObject.getregionmesh(0) # surface=0, here "surface" means the interior

		if npoints is not None:
			# restore original value
			self.astObject.set("MeshSize={0}".format(old_mesh_size))

		return mesh.T
	
	
	# Attributes to implement: Adaptive, Negated, Closed, FillFactor, Bounded
	# See p. 811 of documentation


