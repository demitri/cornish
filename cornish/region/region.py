from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink
import starlink.Ast as Ast

from ..ast_object import ASTObject
from ..mapping import ASTFrame, ASTMapping

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
		elif isinstance(map, ASTMapping):
			ast_map = map.astObject
		else:
			raise Exception("Expected 'map' to be one of these two types: (starlink.Ast.Mapping, ASTMap).")

		if isinstance(frame, starlink.Ast.Frame):
			ast_frame = frame
		elif isinstance(map, ASTFrame):
			ast_frame = frame.astObject
		else:
			raise Exception("Expected 'frame' to be one of these two types: (starlink.Ast.Frame, ASTFrame).")
		
		new_ast_region = self.astObject.mapregion(ast_map, ast_frame)

		return ASTFrame(ast_frame=new_ast_region)
	
