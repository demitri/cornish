from __future__ import (absolute_import, division, print_function, unicode_literals)

from math import radians as deg2rad
from math import degrees as rad2deg
import numpy as np
import starlink
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
		'''
		Parameters
		----------
		centerPoint : numpy.ndarray, list, tuple
			Two elements that describe the center point of the circle in the provided frame
	
		edgePoint : numpy.ndarray, list, tuple
			Two elements that describe a point on the circumference of the circle in the provided frame
	
		radius : float
			The radius of the circle to be created.
		'''
		self._uncertainty = 4.848e-6 # defaults to 1 arcsec
		
		if ast_circle is not None:
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
		#	CENTER_RADIUS (1) : circle specified by  point and radius
		input_form = None
		
		# check valid combination of parameters
		# -------------------------------------
		if frame is None:
			raise Exception("ASTCircle: A frame must be specified when creating an ASTCircle.")
		else:
			if isinstance(frame, ASTFrame):
				self.frame = frame
			elif isinstance(frame, Ast.Frame):
				self.frame = ASTFrame(frame=frame)
			else:
				raise Exception("ASTCircle: unexpected frame type specified ('{0}').".format(type(frame)))
		
		# convert np.array types to lists so that the value can be used in 'any' and 'all' comparisons.
		if isinstance(centerPoint, np.ndarray):
			centerPoint = centerPoint.tolist()
		if isinstance(edgePoint, np.ndarray):
			edgePoint = edgePoint.tolist()
		
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
	def radius(self):
		'''
		The radius of this circle region.
		
		@returns The radius as a geodesic distance in the associated coordinate system in degrees.
		'''
		center = None
		radius = None
		
		( center, radius, some_point_on_circumference ) = self.astObject.circlepars()
		
		return rad2deg(radius)
		
		# TODO: possibly cache this

	@property
	def center(self):
		'''
		The center of this circle region.
		
		Returns
		-------
		@returns A list of points [x,y] that describe the center of the circle in degrees.
		'''
		center = None
		radius = None
		
		( center, radius, some_point_on_circumference ) = self.astObject.circlepars()
		
		# TODO: possibly cache this
		
		return np.rad2deg(center)
			
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
			
