
from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink
import starlink.Ast as Ast
import astropy.units as u
import numpy as np

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame

__all__ = ["ASTCircle"]

class ASTPolygon(ASTRegion):
	
	def __init__(self, points=None, frame=None, polygon=None):
		'''
		ASTPolygon is an ASTRegion that represents a polygon, a collection of vertices on a sphere in a 2D plane.

		Points may be provided as a list of coordinate points, e.g.
			[(x1, y1), (x2, y2), ... , (xn, yn)]
		or as two parallel arrays, e.g.
			[[x1, x2, x3, ..., xn], [y1, y2, y3, ..., yn]]
		
		self.astObject is of type starlink.Ast.Polygon.
		
		@param points Points that describe the polygon, may be a list of pairs of points or two parallel arrays of axis points.
		@param frame The frame the provided points lie in, accepts either ASTFrame or starlink.Ast.frame objects.
		@param polygon Create a new ASTPolygon from an existing (or more likely returned) starlink.Ast.Polygon object.
		@param Returns a new ASTPolygon object.
		'''
		
		if polygon is not None:
			if isinstance(ast_box, starlink.Ast.Polygon):
				if not any([frame, points]):
					self.astObject = polygon
					return
				else:
					raise Exception("ASTPolygon: cannot specify both an ast_box and any other parameter.")
			else:
				raise Exception("ASTPolygon: The ast_box provided was not of type starlink.Ast.Polygon.")
		
		if isinstance(frame, starlink.Ast.Frame):
			ast_frame = frame
		elif isinstance(frame, ASTFrame):
			ast_frame = frame.astObject
		else:
			raise Exception("ASTPolygon: The supplied 'frame' object must either be a starlink.Ast.Frame or ASTFrame object.")
		
		if points is None:
			raise Exception("A list of points must be provided to create a polygon. This doesn't seem like an unreasonable request.")
		
		# The problem with accepting both forms is that the case of two points is ambiguous:
		# [[x1,x2], [y1, y2]]
		# [(x1,y1), (x2, y2}]
		# I'm going to argue that two points does not a polygon make.
		if len(points) == 2 and len(points[0]) == 2:
			raise Exception("There are only two points in this polygon, making the point ordering ambiguous. But is it really a polygon?")
		
		# Internally, the starlink.Ast.Polygon constructor takes the parallel array form of points.
		# starlink.Ast.Polygon( ast_frame, points, unc=None, options=None )

		parallel_arrays = not (len(points[0]) == 2)
		
		if parallel_arrays:
			self.astObject = Ast.Polygon(ast_frame, points)
		else:
			if isinstance(points, np.ndarray):
				self.astObject = Ast.Polygon(ast_frame, points.T)
			else:
				# Could be a list or lists or tuples - reshape into parallel array form
				dim1 = np.zeros(len(points))
				dim2 = np.zeros(len(points))
				for idx, (x, y) in points:
					dim1[idx] = x
					dim2[idx] = y
				
				self.astObject = Ast.Polygon(ast_frame, np.array([dim1, dim2]))
		
	def downsize(self, maxerr=None, maxvert=0):
		'''
		Returns a new ASTPolygon that contains a subset of the vertices of this polygon.
		
		The subset is chosen so that the returned polygon is a good approximation of this polygon,
		within the limits specified. The density of points in the new polygon is greater
		where the curvature of the boundary is greater.
		
		The 'maxerr' parameter set the maximum allowed discrepancy between the original and
		new polygons as a geodesic distance within the polygon's coordinate frame. Setting this to zero
		returns a new polygon with the number of vertices set in "maxvert".
		
		The 'maxvert' parameter set the maximum number of vertices the new polygon can have. If this is
		less than 3, the number of vertices in the returned polygon will be the minimum needed
		to achieve the maximum discrepancy specified by "maxerr". The unardoned value is in radians,
		but accepts Astropy unit objects.
		
		@param maxerr Maximum allowed discrepancy in radians between the original and new polygons as a geodesic distance within the polygon's coordinate frame.
		@param maxvert Maximum allowed number of vertices in the returned Polygon.
		@returns A new ASTPolygon.
		'''
		
		# should find some reasonable default values
		if None in [maxerr, maxvert]:
			raise Exception("ASTPolygon.downsize: Both 'maxerr' and 'maxvert' must be specified.")
		
		ast_polygon = self.astObject.downsize(maxerr, maxvert)
		return ASTPolygon(polygon=ast_polygon)



