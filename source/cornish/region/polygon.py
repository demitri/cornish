
from __future__ import annotations # remove in Python 3.10

import os
import math
import logging
import numbers
import warnings
from typing import Union, Iterable

import starlink.Ast as Ast
import astropy.units as u
import astropy
import numpy as np

from .box import ASTBox
from .region import ASTRegion
from ..mapping import ASTFrame, ASTSkyFrame, ASTFrameSet
from .. import _pyast_bridge as bridge
from .._validation import as_integer

__all__ = ["ASTPolygon"]

logger = logging.getLogger('cornish')

class ASTPolygon(ASTRegion):
	'''
	ASTPolygon is an ASTRegion that represents a polygon, a collection of vertices on a sphere in a 2D plane.

	Accepted signatures for creating an ASTPolygon:

	.. code-block:: python

		p = ASTPolygon(frame, points)
		p = ASTPolygon(fits_header=fits_header, points=points)  # get the frame from the FITS header provided
		p = ASTPolygon(ast_object)           # where ast_object is a starlink.Ast.Polygon object

	Points may be provided in any of the bridge-accepted forms — SkyCoords,
	Quantities, bare arrays (read as degrees on sky frames, native units
	otherwise) — as a list of coordinate points, e.g.

	.. code-block:: python

		[(x1, y1), (x2, y2), ... , (xn, yn)]

	or as two parallel arrays, e.g.

	.. code-block:: python

		[[x1, x2, x3, ..., xn], [y1, y2, y3, ..., yn]]

	A string format that can be parsed as above is also accepted, e.g.:

	.. code-block::

	    "((131.758,5.366),(131.759,3.766),(132.561,3.767),(133.363,3.766),(133.364,5.366),(132.577,5.367))"

	:param ast_object: create a new ASTPolygon from an existing :class:`starlink.Ast.Polygon` object
	:param frame: the frame the provided points lie in — a cornish or ``starlink.Ast`` frame, frame set (current frame governs), or region
	:param points: points that describe the polygon, may be a list of pairs of points or two parallel arrays of axis points
	:param fits_header: a FITS header whose WCS frame set defines the polygon's frame (the points are then sky coordinates in degrees)
	:returns: Returns a new ``ASTPolygon`` object.
	'''
	def __init__(self, ast_object:Ast.Polygon=None,
				 frame:Union[ASTFrame, Ast.Frame, ASTRegion, Ast.Region]=None,
				 points=None,
				 fits_header=None):

		if ast_object is not None:
			if any([frame is not None, points is not None, fits_header is not None]):
				raise ValueError("Cannot specify 'ast_object' along with any other parameter.")
			# test object
			if isinstance(ast_object, Ast.Polygon):
				super().__init__(ast_object=ast_object)
				return
			else:
				raise TypeError("The 'ast_object' provided was not of type starlink.Ast.Polygon.")

		if points is None:
			raise ValueError("A list of points must be provided to create a polygon. This doesn't seem like an unreasonable request.")

		# Get the frame from the FITS header
		wcs_frameset = None
		if fits_header is not None:
			if frame is not None:
				raise ValueError("Provide the frame via the 'frame' parameter or the FITS header, but not both.")
			wcs_frameset = ASTFrameSet.fromFITSHeader(fits_header=fits_header) # raises FrameNotFoundException
			frame = wcs_frameset # the current (sky) frame governs point interpretation

		if frame is None:
			raise ValueError("A frame must be provided (via 'frame' or 'fits_header').")

		ast_frame = bridge._unwrap(frame) # frame sets and regions are legal: the current/encapsulated frame governs

		pts = bridge.to_frame_units(points, ast_frame) # -> (naxes, npoints), the parallel form Ast.Polygon wants

		# A polygon needs at least three vertices: two points do not a polygon
		# make, and the square two-point case is inherently ambiguous (D3 note —
		# polygons refuse it rather than interpret it).
		if pts.shape[1] < 3:
			raise ValueError(f"A polygon requires at least 3 vertices ({pts.shape[1]} provided); "
			                 f"note that two points are ambiguous (pairs vs. parallel axes) and are always refused.")

		super().__init__(ast_object=Ast.Polygon(ast_frame, pts))
		if wcs_frameset is not None:
			self.wcs = wcs_frameset # SPEC-04 §2

	@staticmethod
	def fromFITSFilepath(path:Union[str,os.PathLike]=None, hdu:int=1):
		'''
		Create a polygon bounding the region of a 2D image.

		:param path: the path to the FITS file
		:param hdu: the HDU to open, 1-based: first HDU = 1
		'''
		hdu = as_integer(hdu, "hdu")
		if hdu < 1:
			# guard the 0-based-indexing guess: hdu-1 == -1 would silently
			# wrap around to the LAST extension via Python list indexing
			raise ValueError(f"'hdu' is 1-based — the first HDU is 1 (got {hdu}).")
		with astropy.io.fits.open(path) as hdu_list:
			if hdu > len(hdu_list):
				raise ValueError(f"'hdu' is {hdu} but the file contains only {len(hdu_list)} HDU(s).")
			polygon = ASTPolygon.fromFITSHeader(hdu_list[hdu-1].header)
		return polygon

	@staticmethod
	def fromFITSHeader(header=None, maxerr:astropy.units.Quantity=1.0*u.arcsec, maxvert:int=200):
		'''
		Creates an ASTPolygon in a sky frame that bounds the field described by
		a FITS header. The header must describe a 2D image and contain WCS
		information.

		Implementation (SPEC-04 §3): the pixel bounding box is mapped directly
		into the sky frame with ``astMapRegion``; when AST returns an exact
		polygon that validates against probe points along the pixel boundary,
		that (typically 4-vertex) polygon is returned. For pathological WCS
		(all-sky projections, discontinuity-crossing fields, heavy distortion)
		the method falls back to the boundary-mesh + downsize recipe.

		The returned polygon carries the pixel<->world frame set in ``.wcs``.

		:param header: a FITS header (astropy, fitsio, an array of cards, a dict, or a string)
		:param maxerr: (fallback path) maximum deviation of the polygon from the true field outline
		:param maxvert: (fallback path) maximum number of vertices in the returned polygon
		'''
		if header is None:
			raise ValueError("A FITS header must be provided.")

		from ..channel import ASTFITSChannel # avoids circular import

		fitsChannel = ASTFITSChannel(header=header)

		# create an ASTFrameSet that contains two frames (pixel grid, WCS) and the mapping between them
		wcsFrameSet = fitsChannel.frameSet          # raises FrameNotFoundException when no WCS
		dims = fitsChannel.dimensions               # raises IncompleteHeader when NAXIS cards are missing

		return ASTPolygon._fromParsedWCS(wcsFrameSet, dims, maxerr=maxerr, maxvert=maxvert)

	@staticmethod
	def _fromParsedWCS(wcsFrameSet:ASTFrameSet, dims, maxerr:astropy.units.Quantity=1.0*u.arcsec, maxvert:int=200,
	                   _force_fallback:bool=False):
		'''
		Internal implementation behind :meth:`fromFITSHeader` (and
		:meth:`ASTFITSChannel.boundingPolygon`), taking an already-parsed
		pixel<->world frame set and the pixel dimensions.

		``_force_fallback`` skips the fast path; it exists so tests can exercise
		the mesh recipe on headers whose fast path would validate.
		'''
		if len(dims) != 2:
			raise ValueError(f"Only 2D images are supported (got {len(dims)} dimensions).")
		# validate up front rather than crash deep in the fallback path: maxerr
		# is a NEW parameter of this API, so an angular Quantity is simply
		# required (no bare-number legacy meaning exists to honor here)
		if not isinstance(maxerr, u.Quantity):
			raise TypeError(f"'maxerr' must be an angular astropy Quantity, e.g. 1*u.arcsec "
			                f"(got '{type(maxerr).__name__}').")
		if not maxerr.isscalar or maxerr.unit.physical_type != 'angle':
			raise ValueError(f"'maxerr' must be a scalar angular Quantity (got {maxerr!r}).")
		if not (np.isfinite(maxerr.value) and maxerr.value > 0):
			raise ValueError(f"'maxerr' must be positive and finite (got {maxerr!r}).")

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
		corner1 = [0.5, 0.5] # outer corner of lower left pixel
		corner2 = [dims[0] + 0.5, dims[1] + 0.5]
		pixelbox = ASTBox.fromCorners(frame=wcsFrameSet.baseFrame, corners=(corner1, corner2),
		                              uncertainty=None) # pixel frame: AST's internal default uncertainty

		# ------------------------------------------------------------------
		# Fast path (SPEC-04 §3.1): map the pixel box directly into the sky.
		# On AST 9.3 a plain mapregion returns an exact (typically 4-vertex)
		# sky polygon for well-behaved WCS in milliseconds.
		# ------------------------------------------------------------------
		fast_path_region = None
		if not _force_fallback:
			mapping = wcsFrameSet.astObject.getmapping()
			current_frame = wcsFrameSet.astObject.getframe(Ast.CURRENT)
			mapped = pixelbox.astObject.mapregion(mapping, current_frame)

			if isinstance(mapped, Ast.Polygon):
				fast_path_region = ASTPolygon(ast_object=mapped)
			elif isinstance(mapped, Ast.Box):
				fast_path_region = ASTBox(ast_object=mapped).toPolygon()
			# a Circle (or anything else) cannot represent the field as an exact
			# polygon; fall through to the mesh path with its controlled maxerr

		if fast_path_region is not None:
			# Validate (SPEC-04 §3.1): probe boundary-mesh points of the pixel
			# box, mapped to the sky, for membership. Any failure — including
			# unmappable (Ast.BAD) boundary pixels — falls back to the mesh
			# recipe, which is always correct, just slower.
			K = 32
			old_mesh_size = pixelbox.astObject.get("MeshSize")
			pixelbox.astObject.set(f"MeshSize={K}")
			try:
				probe_mesh_px = pixelbox.astObject.getregionmesh(1) # boundary mesh, pixel coords
			finally:
				pixelbox.astObject.set(f"MeshSize={old_mesh_size}")
			# pyast-internal radians end-to-end: these values never cross a
			# constructor/user boundary, so no bridge call is inserted (M20 rule)
			probe_sky = wcsFrameSet.astObject.tran(probe_mesh_px, True)
			values_ok = bool(np.all(np.isfinite(probe_sky)) and not np.any(probe_sky == Ast.BAD))
			if values_ok and all(mapped.pointinregion(probe_sky[:, i]) for i in range(probe_sky.shape[1])):
				fast_path_region.wcs = wcsFrameSet # SPEC-04 §2
				return fast_path_region

		# ------------------------------------------------------------------
		# Fallback (SPEC-04 §3.2 / SPEC-04A M19): boundary mesh + downsize.
		# The code below is adapted from code originally provided by David Berry.
		# ------------------------------------------------------------------

		#  Map the box into (RA,Dec) and get the (RA,Dec) at a large number of
		#  points evenly distributed around the boundary.
		skybox = pixelbox.regionWithMapping(map=wcsFrameSet, frame=wcsFrameSet)
		mesh_deg = skybox.boundaryPointMesh() # np.array of point pairs in degrees

		#  A field straddling RA=0 would produce a self-crossing flat-frame
		#  polygon; unwrap the longitudes onto a continuous branch first
		#  (SPEC-04 §3.2 discontinuity-awareness). The flat frame is unit-blind,
		#  so out-of-range longitudes are fine; the final sky polygon normalizes.
		lon = mesh_deg[:, 0]
		if lon.max() - lon.min() > 180.0:
			mesh_deg = mesh_deg.copy()
			mesh_deg[:, 0] = np.where(lon > 180.0, lon - 360.0, lon)

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

		# basic frame -> native pass-through: the flat frame HONESTLY holds degrees (M19)
		flatpoly = ASTPolygon(frame=degFlatFrame, points=mesh_deg)

		#  Remove mesh points where the polygon is close to a Cartesian straight
		#  line, and retain them where it deviates from a straight line.
		#  The flat frame is a degree frame, so degrees ARE its native maxerr
		#  units — this .to(u.deg) is a unit read-out, not a hidden conversion (M19).
		downsizedpoly = flatpoly.downsize(maxerr=maxerr.to(u.deg).value, maxvert=maxvert) # -> ASTPolygon

		# "downsizedpoly" is a polygon in a frame with axes in degrees, but is not a sky frame.
		# The bridge converts its degree points into the frame set's sky frame.
		sky_frame_polygon = ASTPolygon(frame=wcsFrameSet, points=downsizedpoly.points)

		#  The order in which the vertices are supplied to the polygon constructor above
		#  defines which side of the polygon boundary is the inside and which is the
		#  outside. Since we do not know the order of the points returned by the
		#  getregionpoints or boundarypointmesh methods, we check now that the
		#  central pixel in the FITS image is "inside" the polygon, and negate
		#  the polygon if it is not.
		(a, b) = wcsFrameSet.astObject.tran( [[dims[0]/2], [dims[1]/2]] )
		# pyast-internal radians end-to-end: this centre never crosses a
		# constructor/user boundary, so no bridge call is inserted (M20 rule)
		center = wcsFrameSet.astObject.norm(np.array([a, b]))

		if not sky_frame_polygon.astObject.pointinregion( center ):
			# The region as constructed covers all points *outside* the FITS area
			# (and is unbounded); rather than negate a flag, recreate the region
			# by reversing the order of the points.
			sky_frame_polygon = ASTPolygon(frame=wcsFrameSet, points=np.flip(downsizedpoly.points, axis=0))

		sky_frame_polygon.wcs = wcsFrameSet # SPEC-04 §2

		# Return a polygon with the same vertices but defined in a SkyFrame
		# rather than a flat Frame.
		return sky_frame_polygon

	@staticmethod
	def fromPointsOnSkyFrame(frame:ASTSkyFrame=None, points:np.ndarray=None, expand_by:astropy.units.quantity.Quantity=20*u.pix): # astropy.coordinates.BaseRADecFrame
		'''
		Create an ``ASTPolygon`` specifically in a sky frame from an array of points.

		Points can be provided in any of the bridge-accepted forms — SkyCoords,
		Quantities, bare values in degrees — either as coordinate pairs, e.g.

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
		if frame is None:
			raise ValueError("A sky frame must be provided.")

		ast_frame = bridge._unwrap(frame)
		if not bridge.is_sky(ast_frame):
			from ..exc import NotASkyRegion
			raise NotASkyRegion(f"The frame provided to fromPointsOnSkyFrame must be a sky frame (got '{ast_frame.Class}').")

		# the code below requires two parallel arrays of ra, dec in radians
		ra_list, dec_list = bridge.to_frame_units(points, ast_frame)

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

		#  Create a PointList holding the (RA,Dec) positions.
		plist = Ast.PointList( ast_frame, np.array([ra_list, dec_list]) )

		#  Get the centre and radius of the circle that bounds the points (in
		#  radians).
		(centre,radius) = plist.getregiondisc()

		#  Determine the number of pixels (M) along each size of the grid that
		#  will produce pixels equal is size of ACC. If a non-zero value for M
		#  has already been set, use it.
		if M == 0 :
		   M = int( 1 + 2.0*radius/ACC )

		#  Create a minimal set of FITS-WCS headers that describe a TAN
		#  projection that projects the above circle into a square of M.M
		#  pixels. The reference point is the centre of the circle and is put
		#  at the centre of the square grid. Put the headers into a FitsChan.
		#  (Ast.DR2D is AST's own radians->degrees constant: FITS cards are
		#  pyast-internal here, never crossing a user boundary.)
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
			new_ast_polygon = big_poly.mapregion( wcs, ast_frame )

		else:
			# Transform the Polygon into (RA,Dec)
			new_ast_polygon = pix_poly.mapregion( wcs, ast_frame )

		polygon = ASTPolygon(ast_object=new_ast_polygon)

		# check if we need to negate polygon; the centre stays pyast-internal
		# radians end-to-end (M20 rule), so test membership on the raw object
		if not polygon.astObject.pointinregion(centre):
			polygon.negate()

		return polygon


	def downsize(self, maxerr=None, maxvert=0):
		'''
		Returns a new ASTPolygon that contains a subset of the vertices of this polygon.

		The subset is chosen so that the returned polygon is a good approximation of this polygon,
		within the limits specified. The density of points in the new polygon is greater
		where the curvature of the boundary is greater.

		The 'maxerr' parameter sets the maximum allowed discrepancy between the original and
		new polygons as a geodesic distance within the polygon's coordinate frame:

		* an :class:`astropy.units.Quantity` (angular) — sky-frame polygons only;
		* a bare number on a NON-sky (e.g. flat degree) frame — the frame's native units;
		* a bare number on a SKY frame — currently read as RADIANS for backward
		  compatibility, which is DEPRECATED: pass a Quantity instead (this form
		  will become a ValueError in a future release).

		Setting maxerr to zero returns a polygon with exactly "maxvert" vertices.

		The 'maxvert' parameter sets the maximum number of vertices the new polygon can have. If this is
		less than 3, the number of vertices in the returned polygon will be the minimum needed
		to achieve the maximum discrepancy specified by "maxerr".

		:param maxerr: maximum allowed discrepancy between the original and new polygons (see above)
		:param maxvert: maximum allowed number of vertices in the returned polygon
		:returns: a new ASTPolygon.
		'''
		if maxerr is None or maxvert is None:
			raise ValueError("Both 'maxerr' and 'maxvert' must be specified.")
		maxvert = as_integer(maxvert, "maxvert")

		# pyast's own downsize silently accepts NaN/BAD/negative/bool maxerr
		# values, so the validation must live here (M33): bool -> TypeError;
		# NaN/inf/Ast.BAD/negative -> ValueError; zero stays legal (AST-documented
		# meaning: "return a polygon with exactly maxvert vertices").
		if isinstance(maxerr, bool):
			raise TypeError("A bool cannot be interpreted as a 'maxerr' distance.")
		if isinstance(maxerr, u.Quantity):
			# sky frames: converted to radians; basic frames: ValueError (native units unknowable, D12)
			maxerr_frame_units = bridge.to_frame_distance(maxerr, self.astObject)
		elif isinstance(maxerr, numbers.Real):
			value = float(maxerr)
			if not math.isfinite(value) or value == Ast.BAD:
				raise ValueError(f"'maxerr' must be finite (got {value!r}).")
			if bridge.is_sky(self.astObject):
				# NEVER silently reinterpret: a bare float here has always meant
				# radians (a documented, working meaning — unlike `unc`, D15's
				# no-compat argument does not apply). Deprecate toward Quantity.
				warnings.warn("Passing a bare number as 'maxerr' to downsize() on a sky-frame "
				              "polygon is read as RADIANS; this will raise a ValueError in a "
				              "future release — pass an astropy Quantity (e.g. maxerr * u.rad) instead.",
				              DeprecationWarning, stacklevel=2)
			# sky frames: legacy radians; basic frames: native frame units, permanently
			maxerr_frame_units = value
		else:
			raise TypeError(f"Cannot interpret an object of type '{type(maxerr).__name__}' as 'maxerr'.")
		if maxerr_frame_units < 0:
			raise ValueError(f"'maxerr' must be non-negative (got {maxerr!r}).")

		ast_polygon = self.astObject.downsize(maxerr_frame_units, maxvert)
		new_polygon = ASTPolygon(ast_object=ast_polygon)
		new_polygon.wcs = self.wcs # propagate the originating WCS, when known (SPEC-04 §2)
		return new_polygon

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
				# reading out the interior ANGLE in radians — a property of the
				# result, not a coordinate conversion (gate-allowlisted)
				angles.append(angle.to(u.rad).value)

			sum_of_polygon_angles = sum(angles) # radians
			area_sr = (sum_of_polygon_angles - (n-2) * math.pi) * u.sr
			return area_sr

		else:

			# Ref: https://mathworld.wolfram.com/PolygonArea.html
			raise NotImplementedError("The area calculation for a polygon in a non-sky frame has not been implemented.")

	def toPolygon(self, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec) -> ASTPolygon:
		'''
		Common interface to return a polygon from a region; here 'self' is returned.

		The parameters 'npoints' and 'maxerr' are ignored.
		'''
		return self
