from __future__ import (absolute_import, division, print_function, unicode_literals)

from ..ast_object import ASTObject
from ..mapping import ASTFrame

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
	
	'''
	def __init__(self, frame=None):
		# this might not make sense?
		if frame is None:
			raise Exception("ASTRegion requires a frame object to be specified.")
		self.frame = frame


