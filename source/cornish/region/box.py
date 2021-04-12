
'''
Notes.

From the AST documentation: "The Box class does not define any new attributes beyond those
which are applicable to all Regions.", i.e. for AST, there is nothing special about a box
beyond a specific means to define it, i.e. a corner and center or two corner points.
'''

from __future__ import annotations # remove in Python 3.10

import math
import numpy as np
import logging
from typing import Union, Iterable
from collections.abc import Collection

import starlink
import starlink.Ast as Ast
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame, ASTFrameSet
#from ..channel.fits_channel import ASTFITSChannel

__all__ = ["ASTBox"]

CENTRE_CORNER = 0
CORNER_CORNER = 1

logger = logging.getLogger("cornish")

def point2array_rad(point:Union[Iterable,SkyCoord]): #, is_sky_coordinate:bool=True):
	'''
	Function that converts a coordinate point to a ``numpy.ndarray`` in radians (if a sky point).
	
	If no units are specific via Quantity, values are assumed to be degrees.
	
	:param point: a two element point in the form of an iterable (e.g. list, tuple) or a SkyCoord
	:param is_sky_coordinate: boolean to indicate whether the point is a celestial coordinate
	'''
	#if not is_sky_coordinate:
	#	return point #np.array([point[0], point[1]])
	
	if isinstance(point, SkyCoord):
		return np.array([point.ra.to(u.rad).value, point.dec.to(u.rad).value])
	elif isinstance(point, u.Quantity):
		# expecting something along the line of np.array([1,2]) * u.deg or [1,2] * u.deg
		return point.to(u.rad).value # get array out
	else:
		# expecting something iterable
		return np.deg2rad(point)

def point2array_deg(point:Union[Iterable,SkyCoord]):
	'''
	Function that converts a coordinate point to a ``numpy.ndarray`` in degrees (if a sky point).
	
	If no units are specific via Quantity, values are assumed to be degrees and returned.

	:param point: a two element point in the form of an iterable (e.g. list, tuple) or a SkyCoord
	'''
	if isinstance(point, SkyCoord):
		return np.array([point.ra.to(u.deg).value, point.dec.to(u.deg).value])
	elif isinstance(point, u.Quantity):
		# expecting something along the line of np.array([1,2]) * u.deg or [1,2] * u.deg
		return point.to(u.deg).value # get array out
	else:
		return point

class ASTBox(ASTRegion):
	'''
	ASTBox is an ASTRegion that represents a box with sides parallel to the axes of an ASTFrame.
	
	Accepted signatures for creating an ASTBox:
	
	.. code-block:: python
	
		b = ASTBox(frame, cornerPoint, cornerPoint2)
		b = ASTBox(frame, cornerPoint, centerPoint)
		b = ASTBox(frame, dimensions)
		b = ASTBox(ast_box) (where ast_box is an Ast.Box object)
	
	Points and dimensions can be any two element container, e.g.
	
	.. code-block:: python
	
		(1,2)
		[1,2]
		np.array([1,2])
		
	If ``dimensions`` is specified, a box enclosing the entire area will be defined.
	
	The 'frame' parameter may either be an ASTFrame object or a :class:`starlink.Ast.frame` object.
	
	A Box is similar to an Interval, the only real difference being that the Interval
	class allows some axis limits to be unspecified. Note, a Box will only look like a box
	if the Frame geometry is approximately flat. For instance, a Box centered close to a pole
	in a SkyFrame will look more like a fan than a box (the Polygon class can be used to
	create a box-like region close to a pole).
	
	:param ast_box: an existing object of type :class:`starlink.Ast.Box`
	:param frame: a frame the box is to be defined in, uses :class:`~cornish.ASTICRSFrame` if `None`
	:param cornerPoint:
	:param cornerPoint2:
	:param centerPoint:
	:param dimensions: dimensions of the box in pixels for use on a Cartesian frame (AST frame='Cartesian' and system='GRID')
	'''
	def __init__(self, ast_object:starlink.Ast.Box=None):
		
		if ast_object.isabox() == False:
			raise ValueError(f"The 'ast_frame' provided must be a starlink.Ast.Box object, got '{type(ast_object)}'.")
		
		super().__init__(ast_object)
		
		# if ast_object.isaframe():
		# 	if ast_object.isaskyframe():
		# 		self.frame = ASTSkyFrame(ast_object=ast_object)
		# 	elif ast_object.isaframeset():
		# 		self.frame = ASTFrameSet(ast_object=ast_object)
		# 	elif ast_object.isacmpframe():
		# 		self.frame = ASTCompoundFrame(ast_object=ast_object)
		# 	#elif ast_object.isafluxframe():
		# 	#	self.frame = ASTFluxFrama(ast_object=ast_object)
		# 	else:
		# 		self.frame = ASTFrame(ast_object=ast_object)
		# else:
		# 	raise ValueError(f"A frame could not be determined from the provided 'ast_object'.")
	
	def __init2__(self, ast_object:starlink.Ast.Box=None, \
		         frame:ASTFrame=None, \
                 cornerPoint:Union[Iterable,SkyCoord]=None, \
				 cornerPoint2:Union[Iterable,SkyCoord]=None, \
				 centerPoint:Union[Iterable,SkyCoord]=None, \
				 dimensions=None):
		#self.astFrame = frame
		self._uncertainty = 4.848e-6 # defaults to 1 arcsec
		#self._ast_box = None
		
		# I am assuming the box is immutable...
		# dimensions = pixel dimensions
		self.dimensions = None
		
		if ast_object is not None:
			if isinstance(ast_object, Ast.Box):
				if not any([cornerPoint, cornerPoint2, centerPoint, dimensions]):
					self.astObject = ast_object
					return
				else:
					raise Exception("ASTBox: Cannot specify both an 'ast_object' and any other parameter.")
			else:
				raise Exception("ASTBox: The 'ast_object' provided was not of type starlink.Ast.Box.")
		
		# input forms:
		#    0: box specified by center point and any corner point
		#    1: box specified by a corner and its opposite corner
		input_form = None
		
		# Get the frame from the FITS header
 		#if fits_header:
 		#	from ..channel import ASTFITSChannel
 		#	# get frame from header
 		#	fits_channel = ASTFITSChannel(header=fits_header)
 		#	
 		#	# does this channel contain a frame set?
 		#	frame_set = fits_channel.frameSet
 		#	#if frame_set is None:
 		#	#	raise ValueError("The provided FITS header does not describe a region (e.g. not an image, does not contain a WCS that AST can read).")
 		#	#else:
 		#
 		#	frame = frame_set.baseFrame
 		#				
 		#	# support n-dimensional boxes
 		#	
 		#	# define needed parameters for box creation below
 		#	dimensions = fits_channel.dimensions
 		#	n_dim = len(dimensions)
 		#	cornerPoint = [0.5 for x in range(n_dim)]
 		#	cornerPoint2 = [dimensions[x] + 0.5 for x in range(n_dim)]
 		#	#cornerPoint=[0.5,0.5], # center of lower left pixel
 		#	#cornerPoint2=[dimensions[0]+0.5, dimensions[1]+0.5])
 		#
 		#	if n_dim > 2:
 		#		raise NotImplementedError("the code below must be written to handle n-dim")
				
		# check valid combination of parameters
		# -------------------------------------
		if frame is None:
			raise Exception("ASTBox: A frame must be specified when creating an ASTBox.")
		else:
			if isinstance(frame, ASTFrame):
				self.frame = frame
			elif isinstance(frame, starlink.Ast.Frame):
				self.frame = ASTFrame(frame=frame)
			else:
				raise Exception("ASTBox: unexpected frame type specified ('{0}').".format(type(frame)))

		if all([x is not None for x in [cornerPoint,centerPoint]]) or \
		   all([x is not None for x in [cornerPoint,cornerPoint2]]) or \
		   dimensions is not None:
			if dimensions is not None:
				input_form = CORNER_CORNER
				p1 = [0.5,0.5] # use 0.5 to specify the center of each pixel
				p2 = [dimensions[0]+0.5, dimensions[1]+0.5]
			elif centerPoint is None:
				input_form = CORNER_CORNER
				p1 = [cornerPoint[0], cornerPoint[1]]
				p2 = [cornerPoint2[0], cornerPoint2[1]]
				dimensions = [math.fabs(cornerPoint[0] - cornerPoint2[0]),
							  math.fabs(cornerPoint[1] - cornerPoint2[1])]
			else:
				input_form = CENTRE_CORNER
				p1 = [centerPoint[0], centerPoint[1]]
				p2 = [cornerPoint[0], cornerPoint[1]]
				dimensions = [2.0 * math.fabs(centerPoint[0] - cornerPoint[0]),
							  2.0 * math.fabs(centerPoint[1] - cornerPoint[1])]

			self.dimensions = [dimensions[0], dimensions[1]]
			#logger.debug("Setting dims: {0}".format(self.dimensions))

		else:
			raise Exception("ASTBox: Either 'cornerPoint' and 'centerPoint' OR 'cornerPoint' " + \
							"and 'cornerPoint2' OR 'dimensions' must be specified when creating an ASTBox.")
	
# 		if input_form == CENTRE_CORNER:
# 			p1 = [centerPoint[0], centerPoint[1]]
# 			p2 = [cornerPoint[0], cornerPoint[1]]
# 		else:
# 			p1 = [cornerPoint[0], cornerPoint[1]]
# 			p2 = [cornerPoint2[0], cornerPoint2[1]]
		
		# Box args: :frame,form,point1,point2,unc=None,options=None  <-- note which are keyword args & which not
		# AstBox( starlink.Ast.Frame(2), [0,1], 
		
		raise Exception("convert points to radians before passing them here")
		
		self.astObject = Ast.Box( self.frame.astObject, input_form, p1, p2, unc=self.uncertainty )
		
 			#if fits_header is not None:
 			#	# create the mapping from pixel to sky (or whatever is there) if available
 			#	mapping = frame_set.astObject.getmapping() # defaults are good
 			#	current_frame = frame_set.astObject.getframe(starlink.Ast.CURRENT)
 			#	
 			#	# create a new region with new mapping
 			#	self.astObject = self.astObject.mapregion(mapping, current_frame)

	@classmethod
	def fromCentreAndCorner(cls, frame:Union[Ast.Frame, ASTFrame], centre:Iterable=None, corner:Iterable=None, center:Iterable=None, uncertainty:Union[ASTRegion, Ast.Region]=None) -> ASTBox:
		'''
		Create a new ASTBox object defined by the provided corner and centre points.
		
		:param frame: the frame the provided points lie in, accepts either :class:`ASTFrame` or :class:`starlink.Ast.Frame` objects
		:param centre: the coordinate of the point at the centre of the box in the frame provided
		:param corner: the coordinate of the point at any corner of the box in the frame provided
		:param center: synonym for 'centre', ignored if 'centre' is defined
		:param uncertainty: 
		'''
		if center and not centre:
			centre = center

		if frame is None:
			raise ValueError("A frame the region is defined in must be provided.")
		elif not all([x is not None for x in [centre, corner]]):
			raise ValueError("Both a corner and centre coordinates must be provided to define this box region.")
				
		if isinstance(frame, Ast.Frame):
			if frame.isaframeset():
				ast_frame = frame.frameAtIndex(Ast.CURRENT) # get current frame
			else:
				ast_frame = frame
		elif isinstance(frame, ASTFrameSet):
			ast_frame.currentFrame().astObject
		elif isinstance(frame, ASTFrame):
			ast_frame = frame.astObject
		else:
			raise Exception(f"Unexpected/unsupported frame type encountered: '{type(frame)}'.")
			
		# convert coordinate points to radians if they are in a sky frame
		if ast_frame.isaskyframe():
			centre_rad = point2array_rad(centre)
			corner_rad = point2array_rad(corner)
		
			if uncertainty:
				ast_box = Ast.Box(ast_frame, CENTRE_CORNER, centre_rad, corner_rad, uncertainty)
			else:
				ast_box = Ast.Box(ast_frame, CENTRE_CORNER, centre_rad, corner_rad)
		else:
			if uncertainty:
				ast_box = Ast.Box(ast_frame, CENTRE_CORNER, centre, corner, uncertainty)
			else:
				ast_box = Ast.Box(ast_frame, CENTRE_CORNER, centre, corner)
				
		return ASTBox(ast_object=ast_box)

	@classmethod
	def fromCorners(cls, frame:Union[Ast.Frame, ASTFrame], corners:Iterable[Iterable]=None, uncertainty:Union[ASTRegion, Ast.Region]=None) -> ASTBox:
		'''
		Create a new ASTBox object defined by two corner points.
		
		:param frame: the frame the provided points lie in, accepts either :class:`ASTFrame` or :class:`starlink.Ast.Frame` objects
		:param corners: a collection (list, tuple, array, etc.) of coordinates of two corners of the box in the frame provided
		:param uncertainty: 
		'''

		c1 = corners[0]
		c2 = corners[1]
		
		#print(corners)
		
		if frame is None:
			raise ValueError("A frame the region is defined in must be provided.")
		elif not all([x is not None for x in [c1, c2]]):
			raise ValueError("Tow corner coordinates must be provided to define this box region.")
				
		if isinstance(frame, Ast.Frame):
			if frame.isaframeset():
				ast_frame = frame.frameAtIndex(Ast.CURRENT) # get current frame
			else:
				ast_frame = frame
		elif isinstance(frame, ASTFrameSet):
			ast_frame.currentFrame().astObject
		elif isinstance(frame, ASTFrame):
			ast_frame = frame.astObject
		else:
			raise Exception(f"Unexpected/unsupported frame type encountered: '{type(frame)}'.")
			
		# convert coordinate points to radians if they are in a sky frame
		if ast_frame.isaskyframe():
			c1_rad = point2array_rad(c1)
			c2_rad = point2array_rad(c2)
		
			if uncertainty:
				ast_box = Ast.Box(ast_frame, CORNER_CORNER, c1_rad, c2_rad, uncertainty)
			else:
				ast_box = Ast.Box(ast_frame, CORNER_CORNER, c1_rad, c2_rad)
		else:
			if uncertainty:
				ast_box = Ast.Box(ast_frame, CORNER_CORNER, c1, c2, uncertainty)
			else:
				ast_box = Ast.Box(ast_frame, CORNER_CORNER, c1, c2)
			
		return ASTBox(ast_object=ast_box)

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
	def centre(self) -> np.ndarray:
		'''
		Returns the location of the Box's center as a coordinate pair, in degrees if a sky frame.
		
		:returns: a Numpy array of points (axis1, axis2)
		'''
		# 'getregionpoints' returns two points as a Numpy array: (center, a corner)
		#center_rad, a_corner_rad = rad2deg(box.astObject.norm(box.astObject.getregionpoints())).T
		#return np.rad2deg(center_rad)
		
		center_point, a_corner_point = self.astObject.getregionpoints().T
		if self.frame().isSkyFrame:
			return np.rad2deg(self.astObject.norm(center_point))
		else:
			return center_point
		
	@property
	def center(self) -> np.ndarray:
		'''
		Returns "self.centre". This is a British library, after all.
		'''
		return self.centre
		
	
	@property
	def corner(self) -> np.ndarray:
		'''
		Returns the location of one of the box's corners as a coordinate pair, in degrees if a sky frame.
		
		:returns: a Numpy array of points (axis1, axis2)
		'''
		
		# !! Is this off by 1 or 0.5  (or neither) due to lower left being [0.5,0.5] ?
		
		center_point, a_corner_point = self.astObject.getregionpoints().T
		
		if self.frame().isSkyFrame:
			return np.rad2deg(self.astObject.norm(a_corner_point))
		else:
			return a_corner_point
	
	def corners(self, mapping=None) -> list:
		'''
		Returns a list of all four corners of box.
		
		Parameters
		----------
		mapping : `ASTMapping`
			A mapping object.
		
		Returns
		-------
		list
			A list of points in degrees, e.g. ``[(p1,p2), (p3, p4), (p5, p6), (p7, p8)]``
		'''
		
		if mapping is None:
			raise Exception("ASTBox: A mapping must be provided to return a list of corner points.")
		
		# create a 2D array of shape points in pixel frame to transform
		
		d1, d2 = self.dimensions[0], self.dimensions[1]
		points = [[0.5    , 0.5],   	# center of lower left pixel
				  [0.5    , d2+1.0], 	# center of upper left pixel
				  [d1+1.0 , d2+1.0], 	# center of upper right pixel
				  [d1+1.0 , 0.5]]    	# center of lower right pixel
		
		# Need to refactor (transpose) points for format that AST expects:
		#     [[x1, x2, x3, x4], [y1, y2, y3, y4]]
		points = np.array(points).T
		
		#x_pos = [p[0] for p in points]
		#y_pos = [p[1] for p in points]
		
		forward = True # True = forward transformation, False = inverse
		corner_points = mapping.astObject.tran(points, forward) # transform points from one frame (pix) to another (WCS)
				
		#logger.debug("Back from tran: {0}".format(corner_points))
		
		# Transpose back to get a list of points: [[x1, y1], [x2, y2], ...]
		corner_points = corner_points.T
		
		#NOTE: The box may not necessary be on a sky frame!
		if self.frame.isSkyFrame():
			# AST returns values in radians; normalize and convert to degrees
			corner_points = np.rad2deg(frame.norm(corner_points))
		
		#logger.debug("box.corners: corner points = {}".format(corner_points))
		
		# normalize RA positions on [0,360)
		for point in corner_points:
			while point[0] < 0.0:
				point[0] += 360
			while point[0] >= 360.0:
				point[0] -= 360.0
		
		#logger.debug("corner_points={0}".format(corner_points))
		
		return corner_points
	
	def toPolygon(self) -> ASTPolygon: #, maxerr:astropy.units.Quantity=1.0*u.arcsec):
		'''
		Returns a four-vertex ASTPolygon that describes this box in the same frame.
		'''
		
		# Note: other region methods use the boundary mesh points technique to get a polygon.
		# A box is a simple enough shape to get a precise polygon from the provided points.
		
		from .polygon import ASTPolygon # avoid circular import
		if isinstance(self, ASTPolygon):
			raise Exception("why is this a polygon?")
			return self
		
		if self.frame().astObject.isaskyframe():
			centre, corner = np.rad2deg(self.astObject.getregionpoints()).T # in radians
		else:
			centre, corner = self.astObject.getregionpoints().T

		d_x, d_y = centre - corner # or d_ra, d_dec if a sky frame
		
		polygon_points = np.array([
			[centre[0] - d_x, centre[1] - d_y],
			[centre[0] - d_x, centre[1] + d_y],
			[centre[0] + d_x, centre[1] + d_y],
			[centre[0] + d_x, centre[1] - d_y]
		])
		
		polygon = ASTPolygon(frame=self.frame(), points=polygon_points)
		if polygon.containsPoint(centre) is False:
			polygon.negate()
		return polygon
				
		#boundary_mesh_points = self.boundaryPointMesh(npoints=npoints)
		#new_polygon = ASTPolygon(frame=self.frame(),
		#						 points=boundary_mesh_points).downsize(maxerr=maxerr.to(u.rad).value)
		#return new_polygon

	@property
	def area(self) -> u.Quantity:
		'''
		The area of the box within its frame (e.g. on a Cartesian plane or sphere). [Not yet implemented.]
		'''
		frame = self.frame() # create variable here as frame() creates a copy
		
		if frame.isSkyFrame:
			
			centre, corner = np.rad2deg(self.astObject.getregionpoints()).T
			
			angles = list()
			polygon_points = self.toPolygon().points
			n = len(polygon_points) # number of sides of polygon
			
			for idx, c in enumerate(polygon_points):
				if idx == 0:
					p1 = polygon_points[-1]
				else:
					p1 = polygon_points[idx-1]
				if idx == n-1:
					p2 = polygon_points[0]
				else:
					p2 = polygon_points[idx+1]
					
				angle = frame.angle(vertex=c, points=(p1,p2)) # -> Quantity
				angles.append(angle.to(u.deg).value)
			
			#print(angles)
			sum_of_polygon_angles = sum(angles)
			area = math.pi / 180 * (sum_of_polygon_angles - (n - 2) * 180)
			return area * u.deg * u.deg
			
		elif frame.system == 'Cartesian' and frame.domain == "GRID":
			
			centre, corner = self.astObject.getregionpoints().T
			d_x, d_y = centre - corner
			return d_x * d_y * u.pixel * u.pixel
		
		else:
			NotImplementedError("The area computation and units for a frame with .system='{frame.system}' and/or .domain='{frame.domain}' has not been written.")
			
		