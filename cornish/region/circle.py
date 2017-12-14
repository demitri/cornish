from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink.Ast as Ast

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame

__all__ = ["ASTCircle"]

CENTER_EDGE = 0
CENTER_RADIUS = 1

class ASTCircle(ASTRegion):
	'''
	ASTCircle is an ASTRegion that represents a circle.
	
	There are two accepted signatures for creating an ASTCircle:
	
	c = ASTCircle(frame, center, edge)
	c = ASTCircle(frame, center, radius)
	
		
	'''
	def __init__(self, ast_circle=None, frame=None, centerPoint=None, edgePoint=None, radius=None):
		
		self._uncertainty = 4.848e-6 # defaults to 1 arcsec
		
		if ast_box is not None:
			if isinstance(ast_circle, starlink.Ast.Circle):
				if not any([centerPoint, frame, radius, edgePoint]):
					self.astObject = ast_circle
					return
				else:
					raise Exception("ASTCircle: cannot specify both an ast_circle and any other parameter.")
			else:
				raise Exception("ASTCircle: The ast_circle provided was not of type starlink.Ast.Circle.")

		# input forms:
		#	CENTER_EDGE   (0) : circle specified by center point and any point on the circumference
		#	CENTER_RADIUS (1) : circle specified by center point and radius
		input_form = None
		
		# check valid combination of parameters
		# -------------------------------------
		if frame is None:
			raise Exception("ASTCircle: A frame must be specified when creating an ASTBox.")
		else:
			if isinstance(frame, ASTFrame):
				self.frame = frame
			elif isinstance(frame, starlink.Ast.frame):
				self.frame = ASTFrame(frame=frame)
		
		if all([centerPoint, edgePoint]) or all([centerPoint, radius]):
			if edgePoint:
				input_form = CENTER_EDGE
			else:
				input_form = CENTER_RADIUS
		else:
			raise Exception("ASTCircle: Either 'centerPoint' and 'edgePoint' OR 'centerPoint' " + \
							"and 'radius' must be specified when creating an ASTCircle.")
			
		if input_form == CENTER_EDGE:
			p1 = [centerPoint[0], centerPoint[1]]
			p2 = [edgePoint[0], edgePoint[1]]
		else:
			p1 = [centerPoint[0], centerPoint[1]]
			p2 = [radius]
		
		self.astObject = Ast.Circle( self.frame.astObject, input_form, p1, p2, unc=self.uncertainty)
			
	@property
	def radius(self):
		'''
		The radius of this circle region.
		
		@returns The radius as a geodesic distance in the associated coordinate system.
		'''
		center = None
		radius = None
		starlink.Ast.CirclePars(self.astObject, center, radius, None)
		
		# TODO: possibly cache this
		
		return radius

	@property
	def center(self):
		'''
		The center of this circle region.
		
		@returns A list of points [x,y] that describe the center of the circle.
		'''
		center = None
		radius = None
		starlink.Ast.CirclePars(self.astObject, center, radius, None)
		
		# TODO: possibly cache this
		
		return center
			
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
			
