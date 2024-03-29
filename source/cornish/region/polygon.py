
from __future__ import annotations # remove in Python 3.10

import os
import re
import ast
import math
import logging
from typing import Union, Iterable

import starlink.Ast as Ast
import astropy.units as u
import astropy
from astropy.coordinates.builtin_frames import ICRS as AstropyICRS
import numpy as np

from .box import ASTBox
from .region import ASTRegion
from ..mapping import ASTFrame, ASTSkyFrame

__all__ = ["ASTPolygon"]

logger = logging.getLogger('cornish')

class ASTPolygon(ASTRegion):
	'''
	ASTPolygon is an ASTRegion that represents a polygon, a collection of vertices on a sphere in a 2D plane.

	Accepted signatures for creating an ASTPolygon:

	.. code-block:: python

		p = ASTPolygon(frame, points)
		p = ASTPolygon(fits_header, points)  # get the frame from the FITS header provided
		p = ASTPolygon(ast_object)           # where ast_object is a starlink.Ast.Polygon object

	Points may be provided as a list of coordinate points, e.g.

	.. code-block:: python

		[(x1, y1), (x2, y2), ... , (xn, yn)]

	or as two parallel arrays, e.g.

	.. code-block:: python

		[[x1, x2, x3, ..., xn], [y1, y2, y3, ..., yn]]

	A string format that can be parsed as above is also accepted, e.g.:

	.. code-block::

	    "((131.758,5.366),(131.759,3.766),(132.561,3.767),(133.363,3.766),(133.364,5.366),(132.577,5.367))"

	:param ast_object: create a new ASTPolygon from an existing :class:`starlink.Ast.Polygon` object
	:param frame: the frame the provided points lie in, accepts either :class:`ASTFrame` or :class:`starlink.Ast.Frame` objects
	:param points: points in degrees that describe the polygon, may be a list of pairs of points or two parallel arrays of axis points
	:returns: Returns a new ``ASTPolygon`` object.
	'''
	def __init__(self, ast_object:Ast.Polygon=None,
				 frame:Union[ASTFrame, Ast.Frame, ASTRegion, Ast.Region]=None,
				 points=None,
				 fits_header=None):

		if ast_object:
			if any([frame, points, fits_header]):
				raise ValueError("ASTPolygon: Cannot specify 'ast_object' along with any other parameter.")
			# test object
			if isinstance(ast_object, Ast.Polygon):
				super().__init__(ast_object=ast_object)
				return
			else:
				raise Exception("ASTPolygon: The 'ast_object' provided was not of type starlink.Ast.Polygon.")

		if points is None:
			raise Exception("A list of points must be provided to create a polygon. This doesn't seem like an unreasonable request.")

		# Get the frame from the FITS header
		if fits_header:
			if frame is not None:
				raise ValueError("ASTPolygon: Provide the frame via the 'frame' parameter or the FITS header, but not both.")

			frame_set = ASTFrameSet.fromFITSHeader(fits_header=fits_header).baseFrame # raises FrameNotFoundException

		if isinstance(frame, Ast.Region):
			# a region returns 'True' for being a frame, so catch this
			# before testing for Ast.Frame below
			ast_frame = frame.getregionframe()
		elif isinstance(frame, ASTRegion):
			ast_frame = frame.astObject.getregionframe()
		elif isinstance(frame, Ast.Frame):
			ast_frame = frame
		elif isinstance(frame, ASTFrame):
			ast_frame = frame.astObject
		else:
			raise Exception(f"ASTPolygon: The supplied 'frame' object must either be a starlink.Ast.Frame or ASTFrame object (got '{type(frame)}').")

		# parse points if provided as string
		if isinstance(points, str):
			# acceptable forms:
			#    "[(ra1, dec1), (ra2, dec2), ..., (ran, decn)]"
			#    "((ra1, dec1), (ra2, dec2), ..., (ran, decn))"
			#
			# We're going to use eval here; make sure input only
			# contains numbers (inc. exponentials, e.g. "12.34e-2"), spaces, and "()[]".
			s = re.sub('[^\d\[\]\(\) ,\.eE+-]', '', points)
			try:
				points = np.array(ast.literal_eval(s), dtype=float)
			except ValueError:
				raise Exception(f"Could not parse provided string into an array of coordinate points: '{points}'")

		if ast_frame.isaskyframe():
			points = np.deg2rad(points)

		# The problem with accepting both forms is that the case of two points is ambiguous:
		# [[x1,x2], [y1, y2]]
		# [(x1,y1), (x2, y2}]
		# I'm going to argue that two points does not a polygon make.
		if len(points) == 2 and len(points[0]) == 2:
			raise Exception("There are only two points in this polygon, making the point ordering ambiguous. But is it really a polygon?")

		# Internally, the starlink.Ast.Polygon constructor takes the parallel array form of points.
		# starlink.Ast.Polygon( ast_frame, points, unc=None, options=None )

		parallel_arrays = not len(points[0]) == 2

		if parallel_arrays:
			self.astObject = Ast.Polygon(ast_frame, points)
		else:
			if isinstance(points, np.ndarray):
				self.astObject = Ast.Polygon(ast_frame, points.T)
			else:
				dim1, dim2 = points.T
				self.astObject = Ast.Polygon(ast_frame, np.array([dim1, dim2]))

	@staticmethod
	def fromFITSFilepath(path:Union[str,os.PathLike]=None, hdu:int=1):
		'''
		Create a polygon bounding the region of a 2D image.

		:param path: the path to the FITS file
		:param hdu: the HDU to open, first HDU = 1
		'''
		with astropy.io.fits.open(path) as hdu_list:
			polygon = ASTPolygon.fromFITSHeader(hdu_list[hdu-1].header)
		return polygon

	@staticmethod
	def fromFITSHeader(header=None, uncertainty=4.848e-6):
		'''
		Creates an ASTPolygon in a sky frame from a FITS header. Header of HDU must be a 2D image and contain WCS information.

		:param header: a FITS header
		:param uncertainty: TODO: parameter not yet used
		'''
		if header is None:
			raise ValueError("A FITS header must be provided.")

		from ..channel import ASTFITSChannel # avoids circular import

		# The code below is adapted from code originally provided by David Berry.
		fitsChannel = ASTFITSChannel(header=header)

		# create an ASTFrameSet that contains two frames (pixel grid, WCS) and the mapping between them
		wcsFrameSet = fitsChannel.frameSet

		# Create a Box describing the extent of the image in pixel coordinates.
		#
		# From David Berry:
		#       "Because of the FITS-WCS standard, the base Frame in a FrameSet read
		#       from a FITS header will always represent FITS pixel coordinates, which
		#       are defined by the FITS-WCS standard so that the bottom left (i.e.
		#       first) pixel in a 2D array is centered at (1,1). That means that
		#       (0.5,0.5) is the bottom left corner of the bottom left pixel, and
		#       (dim1+0.5,dim2+0.5) is the top right corner of the top right pixel.
		#       This results in the Box covering the whole image area."
		#
		dims = fitsChannel.dimensions
		corner1 = [0.5,0.5] # center of lower left pixel
		corner2 = [dims[0]+0.5, dims[1]+0.5]
		pixelbox = ASTBox.fromCorners(frame=wcsFrameSet.baseFrame, corners=(corner1, corner2))
						  #cornerPoint=[0.5,0.5], # center of lower left pixel
						  #cornerPoint2=[dims[0]+0.5, dims[1]+0.5])

		#  Map this box into (RA,Dec)
		#
		skybox = pixelbox.regionWithMapping(map=wcsFrameSet, frame=wcsFrameSet) # -> ASTBox

		#  Get the (RA,Dec) at a large number of points evenly distributed around
		#  the polygon. The number of points created is controlled by the
		#  MeshSize attribute of the polygon.
		#
		mesh_deg = skybox.boundaryPointMesh() # np.array of points in degrees

		#  Create a polygon using the vertices in the mesh. This polygon is
		#  defined in a basic Frame (flat geometry) - not a SkyFrame (spherical
		#  geometry). If we used a SkyFrame, then all the mesh points along each
		#  edge of the box would fall exactly on a geodesic (i.e. a great circle),
		#  and so the subsequent call to the downsize function would remove them all
		#  (except the corners). Using a basic Frame means that the downsize function
		#  will use geodesics that are Cartesian straight lines. So points that
		#  deviate by more than the required error form a Cartesian straight line
		#  will be retained by the downsize function.
		#
		degFlatFrame = ASTFrame(naxes=2)
		degFlatFrame.setUnitForAxis(axis=1, unit="deg")
		degFlatFrame.setUnitForAxis(axis=2, unit="deg")

		flatpoly = ASTPolygon(frame=degFlatFrame, points=np.deg2rad(mesh_deg))

		#  Remove mesh points where the polygon is close to a Cartesian straight
		#  line, and retain them where it deviates from a straight line, in order
		#  to achieve an max error of 1 arc-second (4.8E-6 rads).
		#
		downsizedpoly = flatpoly.downsize(maxerr=4.848e-6, maxvert=200) # -> ASTPolygon

		#logger.debug(f"{flatpoly=}")

		# "downsizedpoly" is a polygon in a frame with axes in degrees, but is not a sky frame.

		sky_frame_polygon = ASTPolygon(frame=wcsFrameSet, points=downsizedpoly.points)

		#  Create a polygon with the same vertices but defined in a SkyFrame rather than a flat Frame.
		#downsized_points = downsizedpoly.astObject.norm(downsizedpoly.astObject.getregionpoints())
#		sky_frame_polygon = ASTPolygon(frame=wcsFrameSet,
#			                           points=np.rad2deg(downsizedpoly.astObject.getregionpoints()))

		#  The order in which the vertices are supplied to the polygon constructor above
		#  defines which side of the polygon boundary is the inside and which is the
		#  outside. Since we do not know the order of the points returned by the
		#  getregionpoints or boundarypointmesh methods, we check now that the
		#  central pixel in the FITS image is "inside" the polygon, and negate
		#  the polygon if it is not.
		(a, b) = wcsFrameSet.astObject.tran( [[dims[0]/2], [dims[1]/2]] )
		center = wcsFrameSet.astObject.norm(np.array([a,b]))

		#logger.debug(f"before: {sky_frame_polygon.isBounded=}, {center=}, {a=}, {b=}")

		#if not sky_frame_polygon.astObject.pointinregion( [a[0], b[0]] ):
		if not sky_frame_polygon.astObject.pointinregion( center ):
			# sky_frame_polygon.astObject.negate()
			# The region is all points *outside* of the FITS area.
			# In this case the region is unbounded. One option is to call .negate(),
			# but this only flips a boolean flag in the object and leaves the region unbounded.
			# Instead, just recreate the region by reversing the order of the points.
			sky_frame_polygon = ASTPolygon(frame=wcsFrameSet, points=np.flip(downsizedpoly.points, axis=0))

		#logger.debug(f"after: {sky_frame_polygon.isBounded=}")

		# Return a polygon with the same vertices but defined in a SkyFrame
		# rather than a flat Frame.
		return sky_frame_polygon #ASTPolygon(frame=wcsFrameSet, points=downsizedPolygon.astObject.getregionpoints())

	@staticmethod
	def fromPointsOnSkyFrame(frame:ASTSkyFrame=None, points:np.ndarray=None, expand_by:astropy.units.quantity.Quantity=20*u.pix): # astropy.coordinates.BaseRADecFrame
		'''
		Create an ``ASTPolygon`` specifically in a sky frame from an array of points.

		Points can be provided in degrees either as an array or coordinate pairs, e.g.

		.. code-block:: python

		    np.array([[1,2], [3,4], [5,6]])

		or as parallel arrays of ra,dec:

		.. code-block:: python

		    np.array([[1,3,5], [2,4,6]])

		:param points: coordinate points, either as a list of coordinate pairs or two parallel ra,dec arrays
		:param frame: the frame the points lie in, specified as an ``ASTSkyFrame`` object
		:param expand_by: number of pixels to extend the polygon beyond the provided points
		:returns: new ``ASTPolygon`` object
		'''
		if points is None:
			raise ValueError("Coordinate points must be provided to create a polygon.")

		# .. todo:: handle various input types (e.g. Quantity) (also see below)

		# code below requires two parallel arrays of ra, dec in radians
		points = np.deg2rad(points)
		if len(points[0]) == 2:
			# data is array of coordinate pairs, convert to parallel arrays
			ra_list, dec_list = points.T
		else:
			ra_list, dec_list = points

		# author: David Berry (any errors are Demitri Muna's)
		#
		#  This method uses astConvex to find the shortest polygon enclosing a
		#  set of positions on the sky. The astConvex method determines the
		#  required polygon by examining an array of pixel values, so we first
		#  need to create a suitable pixel array. An (M,M) integer array is first
		#  created and initialised to hold zero axt every pixel. A tangent plane
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

		# .. todo:: handle various input types (e.g. Quantity)
		# if isinstance(points, np.ndarray):
		# 	if len(points.shape) != 2 or points.shape[1] != 2:
		# 		raise ValueError("The shape of the array provided should be (n,2).")
		# 	ra_list, dec_list = np.deg2rad(points.T)
		# elif isinstance(frame, (ASTSkyFrame, Ast.SkyFrame)):
		# 	# if it's a sky frame of some kind, we will expect degrees
		# 	ra_list  = np.deg2rad(ra)
		# 	dec_list = np.deg2rad(dec)

		# convert frame parameter to an Ast.Frame object
		if isinstance(frame, ASTFrame):
			frame = frame.astObject
		elif isinstance(frame, Ast.Frame):
			pass
		else:
			raise ValueError(f"The 'frame' parameter must be either an Ast.SkyFrame or ASTSkyFrame object; got {type(frame)}")

		#  Create a PointList holding the (RA,Dec) positions.
		plist = Ast.PointList( frame, [ra_list, dec_list] )

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
		#  positions. This Polygon is defined in pixel coordinates within the
		#  grid defined by the above FITS headers.
		pix_poly = Ast.convex( 1, Ast.EQ, ar, [1,1], [M,M], False )


		if expand_by.to_value(u.pix) > 0:

			n_pixels = expand_by.to_value(u.pix)

			#  Now expand the above polygon a bit. First get the vertex positions
			#  from the Polygon.
			(x_list, y_list ) = pix_poly.getregionpoints()

			# Transform the centre position from sky to pixel coordinates.
			( x_cen, y_cen ) = wcs.tran( [[centre[0]], [centre[1]]], False )

			#  For each vertex, extend its radial vector by n pixels. Create lists
			#  of extended x and y vertex positions. [Expanding about the centroid of
			#  the original vertices may give better results than expanding about the
			#  centre of the bounding disc in some cases].
			x_new = []
			y_new = []
			for (x,y) in zip( x_list, y_list ):
			   dx = x - x_cen[0]
			   dy = y - y_cen[0]
			   old_radius = math.sqrt( dx*dx + dy*dy )
			   new_radius = old_radius + n_pixels
			   factor = new_radius/old_radius
			   dx *= factor
			   dy *= factor
			   x_new.append( dx + x_cen[0] )
			   y_new.append( dy + y_cen[0] )

			#  Create a new polygon in pixel coordinates using the extended vertex positions.
			big_poly = Ast.Polygon( wcs.getframe( Ast.BASE ), [ x_new, y_new ] )

			# Transform the Polygon into (RA,Dec).
			new_ast_polygon = big_poly.mapregion( wcs, frame )

		else:
			# Transform the Polygon into (RA,Dec)
			new_ast_polygon = pix_poly.mapregion( wcs, frame )

		polygon = ASTPolygon(ast_object=new_ast_polygon)

		# check if we need to negate polygon
		if not polygon.containsPoint(np.rad2deg(centre)):
			polygon.negate()

		return polygon


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
		to achieve the maximum discrepancy specified by "maxerr". The unadorned value is in radians,
		but accepts Astropy unit objects.

		:param maxerr: maximum allowed discrepancy in radians between the original and new polygons as a geodesic distance within the polygon's coordinate frame
		:param maxvert: maximum allowed number of vertices in the returned polygon
		:returns: a new ASTPolygon.
		'''

		# should find some reasonable default values
		if None in [maxerr, maxvert]:
			raise Exception("ASTPolygon.downsize: Both 'maxerr' and 'maxvert' must be specified.")

		ast_polygon = self.astObject.downsize(maxerr, maxvert)
		return ASTPolygon(ast_object=ast_polygon)

	@property
	def area(self) -> astropy.units.quantity.Quantity:
		'''
		Returns the area of the polygon as an :class:`astropy.units.quantity.Quantity`. [Not yet implemented for non-sky frames.]
		'''
		# See: https://stackoverflow.com/questions/1340223/calculating-area-enclosed-by-arbitrary-polygon-on-earths-surface
		# Multiple useful answers on that page.
		# see also: https://github.com/anutkk/sphericalgeometry/blob/master/sphericalgeometry/highlevel.py#L145
		#           https://github.com/spacetelescope/spherical_geometry/blob/fbdc54aa5e5953c5b22723c0982a5f0b45ab2d39/spherical_geometry/polygon.py#L525

		frame = self.frame() # create variable here as frame() creates a copy

		if frame.isSkyFrame:

			# This is Girard's theorem.

			angles = list() # unit: radians
			points = self.points
			n = len(points)

			for idx, v in enumerate(points):
				if idx == 0:
					p1 = points[-1]
				else:
					p1 = points[idx-1]
				if idx == n-1:
					p2 = points[0]
				else:
					p2 = points[idx+1]

				angle = frame.angle(vertex=v, points=(p1,p2)) # -> Quantity
				angles.append(angle.to(u.rad).value)

				print(p1,v, p2, angle)

			#print(angles)
			sum_of_polygon_angles = sum(angles) # radians
			#area = math.pi/180 * (sum_of_polygon_angles - (n-2) * 180)
			#area = math.pi/180 * (sum_of_polygon_angles - (n-2))
			#area = np.deg2rad(sum_of_polygon_angles) - (n-2) * math.pi
			area_sr = (sum_of_polygon_angles - (n-2) * math.pi) * u.sr
			return area_sr

		else:

			# Ref: https://mathworld.wolfram.com/PolygonArea.html
			raise NotImplementedError("The area calculation for a polygon in a non-sky frame has not been immplemented.")

	def toPolygon(self, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec) -> ASTPolygon:
		'''
		Common interface to return a polygon from a region; here 'self' is returned.

		The parameters 'npoints' and 'maxerr' are ignored.
		'''
		return self
