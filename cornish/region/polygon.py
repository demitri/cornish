
import math

import starlink
import starlink.Ast as Ast
import astropy.units as u
import astropy
import numpy as np

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame, ASTSkyFrame

__all__ = ["ASTCircle"]

class ASTPolygon(ASTRegion):
	
	def __init__(self, points=None, frame=None, ast_polygon=None):
		'''
		ASTPolygon is an ASTRegion that represents a polygon, a collection of vertices on a sphere in a 2D plane.

		Points may be provided as a list of coordinate points, e.g.
			[(x1, y1), (x2, y2), ... , (xn, yn)]
		or as two parallel arrays, e.g.
			[[x1, x2, x3, ..., xn], [y1, y2, y3, ..., yn]]
		
		self.astObject is of type :class:`starlink.Ast.Polygon`.
		
		:param points: Points that describe the polygon, may be a list of pairs of points or two parallel arrays of axis points.
		:param frame: The frame the provided points lie in, accepts either ASTFrame or starlink.Ast.frame objects.
		:param ast_polygon: Create a new ASTPolygon from an existing (or more likely returned) starlink.Ast.Polygon object.
		:returns: Returns a new ASTPolygon object.
		'''
		
		if ast_polygon is not None:
			if isinstance(ast_polygon, starlink.Ast.Polygon):
				if not any([frame, points]):
					self.astObject = ast_polygon
					return
				else:
					raise Exception("ASTPolygon: cannot specify both an ast_polygon and any other parameter.")
			else:
				raise Exception("ASTPolygon: The ast_polygon provided was not of type starlink.Ast.Polygon.")
		
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
	
	@staticmethod
	def fromPoints(radec_pairs:np.ndarray=None, ra=None, dec=None, system=None, skyframe:Ast.SkyFrame=None, expand_by=20*u.pix): # astropy.coordinates.BaseRADecFrame
		'''
		Create an ASTPolygon from an array of points.
		
		:param ra: list of RA points, must be in degrees (or :class:`astropy.units.Quantity` objects)
		:param dec: list of declination points, must be in degrees (or :class:`astropy.units.Quantity` objects)
		:returns: new ASTPolygon object
		'''
		# author: David Berry
		#
		#  This method uses astConvex to find the shortest polygon enclosing a
		#  set of positions on the sky. The astConvex method determines the
		#  required polygon by examining an array of pixel values, so we first
		#  need to create a suitable pixel array. An (M,M) integer array is first
		#  created and initialised to hold zero at every pixel. A tangent plane
		#  projection is then determined that maps the smallest circle containing
		#  the specified (RA,Dec) positions onto the grid of (M,M) pixels. This
		#  projection is then used to convert each (RA,Dec) position into a pixel
		#  position and a value of 1 is poked into the array at each such pixel
		#  position. The astConvex method is then used to determine the shortest
		#  polygon that encloses all pixels that have value 1 in the array.
		#
		#  This polygon is then modified by moving each vertex 20 pixels radially
		#  away from the centre of the bounding disc used to define the extent of
		#  the pixel grid.
		#
		#  Finally, the resulting polygon is mapping from pixel coordinates to
		#  (RA,Dec).
		
		#  Set the required positional accuracy for the polygon vertices, given
		#  as an arc-distance in radians. The following value corresponds to 10
		#  arc-seconds. The size of the array (M) is selected to give pixels
		#  that have this size. Alternatively, specify a non-zero value for M
		#  explicitly, in which case the pixel size will be determined from M.
		ACC = 4.85E-5
		M = 0
		
		#  A SkyFrame describing the (RA,Dec) values.
		#skyfrm = Ast.SkyFrame( "System=FK5,Equinox=J2000,Epoch=1982.0" )
		
		#  The RA values (radians).
# 		ra_list = [ 0.1646434, 0.1798973, 0.1925398, 0.2024329, 0.2053291,
# 		            0.1796907, 0.1761278, 0.1701603, 0.1762123, 0.1689954,
# 		            0.1725925, 0.1819018, 0.1865827, 0.19369, 0.1766037 ]
# 		
# 		#  The Dec values (radians).
# 		dec_list = [ 0.6967545, 0.706133, 0.7176528, 0.729342, 0.740609,
# 		             0.724532, 0.7318467, 0.7273944, 0.7225725, 0.7120513,
# 		             0.7087136, 0.7211723, 0.7199059, 0.7268493, 0.7119532 ]

		# .. todo:: handle various input types (np.ndarray, Quantity)
		ra_list = ra
		dec_list = dec
		
		if isinstance(skyframe, ASTSkyFrame):
			skyframe = skyframe.astObject
		
		#  Create a PointList holding the (RA,Dec) positions.
		plist = Ast.PointList( skyframe, [ra_list, dec_list] )
		
		#  Get the centre and radius of the circle that bounds the points (in
		#  radians).
		(centre,radius) = plist.getregiondisc()
		
		#  Determine the number of pixels (M) along each size of the grid that
		#  will produce pixels equal is size of ACC. If a non-zero value for M
		#  has already been set, use it.
		if M == 0 :
		   M = int( 1 + 2.0*radius/ACC )
		#logger.debug(f"Using grid size {M}")
		
		#  Create a minimal set of FITS-WCS headers that describe a TAN
		#  projection that projects the above circle into a square of M.M
		#  pixels. The reference point is the centre of the circle and is put
		#  at the centre of the square grid. Put the headers into a FitsChan.
		fc = Ast.FitsChan()
		fc["NAXIS1"] = M
		fc["NAXIS2"] = M
		fc["CRPIX1"] = 0.5*( 1 + M )
		fc["CRPIX2"] = 0.5*( 1 + M )
		fc["CRVAL1"] = centre[0]*Ast.DR2D
		fc["CRVAL2"] = centre[1]*Ast.DR2D
		fc["CDELT1"] = 2.0*radius*Ast.DR2D/( M - 1 )
		fc["CDELT2"] = 2.0*radius*Ast.DR2D/( M - 1 )
		fc["CTYPE1"] = 'RA---TAN'
		fc["CTYPE2"] = 'DEC--TAN'
		
		#  Re-wind the FitsChan and read the FrameSet corresponding to the above
		#  FITS headers.
		fc.clear("Card")
		wcs = fc.read()
		
		#  Use this FrameSet to transform all the (RA,Dec) positions into pixel
		#  coordinates within the grid.
		( x_list, y_list ) = wcs.tran( [ra_list, dec_list], False )
		
		#  Create an integer numpy array of the same shape, filled with zeros.
		ar = np.zeros( shape=(M,M), dtype=int )
		
		#  Poke a value 1 into the above array at each pixel position, checking
		#  each such position is inside the array.
		for (x,y) in zip( x_list, y_list ):
		   ix = int( round( x ) )
		   iy = int( round( y ) )
		   if ix >= 1 and ix <= M and iy >= 1 and iy <= M:
		      ar[ iy - 1, ix - 1 ] = 1
		
		#  Create a Polygon representing the convex hull that encloses the
		#  positions. This Polygon is defined in pixel coordinaates within the
		#  grid defined by the above FITS headers.
		pix_poly = Ast.convex( 1, Ast.EQ, ar, [1,1], [M,M], False )
		
		
		if expand_by.to_value(u.pix) > 0:
			#  Now expand the above polygon a bit. First get the vertex positions
			#  from the Polygon.
			(x_list, y_list ) = pix_poly.getregionpoints()
			
			# Transform the centre position from sky to pixel coordinates.
			( x_cen, y_cen ) = wcs.tran( [[centre[0]], [centre[1]]], False )
			
			#  For each vertex, extend it's radial vector by 20 pixels. Create lists
			#  of extended x and y vertex positions. [Expanding about the centroid of
			#  the original vertices may give better results than expanding about the
			#  centre of the bounding disc in some cases].
			x_new = []
			y_new = []
			for (x,y) in zip( x_list, y_list ):
			   dx = x - x_cen[0]
			   dy = y - y_cen[0]
			   old_radius = math.sqrt( dx*dx + dy*dy )
			   new_radius = old_radius + 20
			   factor = new_radius/old_radius
			   dx *= factor
			   dy *= factor
			   x_new.append( dx + x_cen[0] )
			   y_new.append( dy + y_cen[0] )
		
			#  Create a new polygon in pixel coordinates using the extended vertex positions.
			big_poly = Ast.Polygon( wcs.getframe( Ast.BASE ), [ x_new, y_new ] )
			
			# Transform the Polygon into (RA,Dec).
			new_ast_polygon = big_poly.mapregion( wcs, skyframe )
		
		else:
			# Transform the Polygon into (RA,Dec)
			new_ast_polygon = pix_poly.mapregion( wcs, skyframe )
		
		return ASTPolygon(ast_polygon=new_ast_polygon)
		
	
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
		return ASTPolygon(ast_polygon=ast_polygon)



