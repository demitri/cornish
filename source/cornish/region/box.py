
import math
import numpy as np
import logging

import starlink
import starlink.Ast as Ast

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame
#from ..channel.fits_channel import ASTFITSChannel

__all__ = ["ASTBox"]

CENTER_CORNER = 0
CORNER_CORNER = 1

logger = logging.getLogger("cornish")

class ASTBox(ASTRegion):
	'''
	ASTBox is an ASTRegion that represents a box with sides parallel to the axes of an ASTFrame.
	
	Accepted signatures for creating an ASTBox:
	
	.. code-block:: python
	
		b = ASTBox(frame, cornerPoint, cornerPoint2)
		b = ASTBox(frame, cornerPoint, centerPoint)
		b = ASTBox(frame, dimensions)
		b = ASTBox(ast_box) (where ast_box is an Ast.Box object)
	
	Points and dimensions can be any two element container, e.g. (1,2), [1,2], np.array([1,2])
	If "dimensions" is specified, a box enclosing the entire area will be defined.
	
	The 'frame' parameter may either be an ASTFrame object or a starlink.Ast.frame object.
	
	A Box is similar to an Interval, the only real difference being that the Interval
	class allows some axis limits to be unspecified. Note, a Box will only look like a box
	if the Frame geometry is approximately flat. For instance, a Box centered close to a pole
	in a SkyFrame will look more like a fan than a box (the Polygon class can be used to
	create a box-like region close to a pole).
	
	:param ast_box: An existing object of type starlink.Ast.Box.
	:param frame:
	:param cornerPoint:
	:param cornerPoint2:
	:param centerPoint:
	:param dimensions:
	'''
	def __init__(self, ast_object:starlink.Ast.Box=None, frame=None, cornerPoint=None, cornerPoint2=None, centerPoint=None, dimensions=None):
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
# 		if fits_header:
# 			from ..channel import ASTFITSChannel
# 			# get frame from header
# 			fits_channel = ASTFITSChannel(header=fits_header)
# 			
# 			# does this channel contain a frame set?
# 			frame_set = fits_channel.frameSet
# 			#if frame_set is None:
# 			#	raise ValueError("The provided FITS header does not describe a region (e.g. not an image, does not contain a WCS that AST can read).")
# 			#else:
# 
# 			frame = frame_set.baseFrame
# 						
# 			# support n-dimensional boxes
# 			
# 			# define needed parameters for box creation below
# 			dimensions = fits_channel.dimensions
# 			n_dim = len(dimensions)
# 			cornerPoint = [0.5 for x in range(n_dim)]
# 			cornerPoint2 = [dimensions[x] + 0.5 for x in range(n_dim)]
# 			#cornerPoint=[0.5,0.5], # center of lower left pixel
# 			#cornerPoint2=[dimensions[0]+0.5, dimensions[1]+0.5])
# 	
# 			if n_dim > 2:
# 				raise NotImplementedError("the code below must be written to handle n-dim")
				
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

			if all([cornerPoint,centerPoint]) or all([cornerPoint,cornerPoint2]) or dimensions is not None:
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
					input_form = CENTER_CORNER
					p1 = [centerPoint[0], centerPoint[1]]
					p2 = [cornerPoint[0], cornerPoint[1]]
					dimensions = [2.0 * math.fabs(centerPoint[0] - cornerPoint[0]),
								  2.0 * math.fabs(centerPoint[1] - cornerPoint[1])]
	
				self.dimensions = [dimensions[0], dimensions[1]]
				#logger.debug("Setting dims: {0}".format(self.dimensions))
	
			else:
				raise Exception("ASTBox: Either 'cornerPoint' and 'centerPoint' OR 'cornerPoint' " + \
								"and 'cornerPoint2' OR 'dimensions' must be specified when creating an ASTBox.")
		
	# 		if input_form == CENTER_CORNER:
	# 			p1 = [centerPoint[0], centerPoint[1]]
	# 			p2 = [cornerPoint[0], cornerPoint[1]]
	# 		else:
	# 			p1 = [cornerPoint[0], cornerPoint[1]]
	# 			p2 = [cornerPoint2[0], cornerPoint2[1]]
			
			# Box args: :frame,form,point1,point2,unc=None,options=None  <-- note which are keyword args & which not
			# AstBox( starlink.Ast.Frame(2), [0,1], 
			self.astObject = Ast.Box( self.frame.astObject, input_form, p1, p2, unc=self.uncertainty )
			
# 			if fits_header is not None:
# 				# create the mapping from pixel to sky (or whatever is there) if available
# 				mapping = frame_set.astObject.getmapping() # defaults are good
# 				current_frame = frame_set.astObject.getframe(starlink.Ast.CURRENT)
# 				
# 				# create a new region with new mapping
# 				self.astObject = self.astObject.mapregion(mapping, current_frame)


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
	def center(self) -> np.ndarray:
		'''
		Returns the location of the Box's center as a coordinate pair (tuple).
		
		:returns: a Numpy array of points (axis1, axis2)
		'''
		return self.astObject.getregionpoints()[0] # returns two points as a Numpy array: (center, a corner)
	
	@property
	def corner(self) -> np.ndarray:
		'''
		Returns the location of one of the Box's corners as a coordinate pair.
		
		:returns: a Numpy array of points (axis1, axis2).
		'''
		
		# !! Is this off by 1 or 0.5  (or neither) due to lower left being [0.5,0.5] ?
		
		return self.astObject.getregionpoints()[1] # returns two points as a Numpy array: (center, a corner)
	
	def corners(self, mapping=None):
		'''
		Returns a list of all four corners of box.
		
		Parameters
		----------
		mapping : `ASTMapping`
			A mapping object.
		
		Returns
		-------
		list
			A list of points: [(p1,p2), (p3, p4), (p5, p6), (p7m p8)] in degrees
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
		
		# AST returns values in radians, convert to degrees
		corner_points = np.rad2deg(corner_points)
		
		#logger.debug("box.corners: corner points = {}".format(corner_points))
		
		# normalize RA positions on [0,360)
		for point in corner_points:
			while point[0] < 0.0:
				point[0] += 360
			while point[0] >= 360.0:
				point[0] -= 360.0
		
		#logger.debug("corner_points={0}".format(corner_points))
		
		return corner_points
	
	@property
	def area(self):
		raise NotImplementedError()
			
			
		