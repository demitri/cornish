#/usr/bin/env python

from typing import Union

import numpy as np
import starlink
import starlink.Ast as Ast

from ... import ASTObject
from .frame import ASTFrame
from ..mapping import ASTMapping
from ...exc import FrameNotFoundException
#from ...channel import ASTFITSChannel

__all__ = ['ASTFrameSet']

class ASTFrameSet(ASTFrame):
	'''
	Create a new AST frame set.
	Object can be created from an :class:`starlink.Ast.FrameSet` "primitive"
	(e.g. returned by another object).

	A set of inter-related coordinate systems made up of existing mapping's and frames.
	A FrameSet may be extended by adding a new Frame and associated Mapping.
	
	A FrameSet must have a "base" frame which represents the "native" coordinate system
	(for example, the pixel coordinates of an image). Similarly, one Frame is termed the
	current Frame and represents the "currently-selected" coordinates. It might typically
	be a celestial or spectral coordinate system and would be used during interactions
	with a user, as when plotting axes on a graph or producing a table of results.
	Other Frames within the FrameSet represent a library of alternative coordinate systems
	which a software user can select by making them current.
	
	Accepted signatures for creating an ``ASTFrameSet``:
	
	.. code-block::python
	
	    fs = ASTFrameSet(ast_object)
	    fs = ASTFrameSet(base_frame)
	
	:param ast_object: an :class:`Ast.astFrame` object from the `starlink-pyast` library
	:param base_frame: base frame to create the FrameSet from
	'''
	def __init__(self, ast_object:starlink.Ast.FrameSet=None, base_frame:Union[starlink.Ast.FrameSet,ASTFrame]=None): #, fits_header=None):
		# validate parameters
		if all([ast_object, base_frame]):
			raise ValueError("Both an 'ast_object' and a 'base_frame' cannot be specified.")
		elif all([x is None for x in [ast_object, base_frame]]):
			raise ValueError("One of 'ast_object' or 'base_frame' must be provided.")
		
		if ast_object:
			if isinstance(ast_object, starlink.Ast.FrameSet):
				super().__init__(ast_object=ast_object)
			else:
				Exception("ASTFrameSet: Unhandled ast_object type ('{0}')".format(ast_object))
		else:
			# construct from provided base_frame
			if isinstance(base_frame, starlink.Ast.Frame):
				fs = Ast.FrameSet(frame=base_frame, options=None)
			elif isinstance(base_frame, ASTFrame):
				fs = Ast.FrameSet(frame=base_frame.astObject, options=None)
			else:
				Exception("ASTFrameSet: Unhandled base_frame type ('{0}')".format(base_frame))

			super().__init__(ast_object=fs)
					
	@staticmethod
	def fromFrames(frame1:Union[ASTFrame, starlink.Ast.Frame], frame2:Union[ASTFrame, starlink.Ast.Frame]):
		'''
		Static method to create a frame set from two existing frames.
		'''
		# get Ast objects for each frame
		ast_frames = list()
		
		for frame in [frame1,frame2]:
			if isinstance(frame, ASTFrame):
				ast_frames.append(frame.astObject)
			elif isinstance(frame, starlink.Ast.AstFrame):
				ast_frames.append(frame)
			else:
				raise ValueError(f"The provided frames must either be of type ASTFrame or starlink.Ast.Frame (got '{type(frame)}'.")
		
		frame_set = Ast.Frame.convert(ast_frames[0], ast_frames[1])
		if frame_set is None:
			# I don't know if this is actually the failure mode
			raise Exception("An ASTFrameSet could not be created since the mapping between the two provided frames could not be determined (conversion may not be possible).")
		else:
			return ASTFrameSet(ast_object=frame_set)
	
	@staticmethod
	def fromFITSHeader(fits_header=None):
		'''
		Static method that returns a FrameSet object read from the provided FITS header.
		'''
		
		if fits_header is None:
			raise ValueError("A FITS header must be provided.")
		
		from ...channel import ASTFITSChannel # or "from cornish import ..." ?

		fits_channel = ASTFITSChannel(header=fits_header)
		
		# does this channel contain a frame set?
		frame_set = fits_channel.frameSet
		if frame_set:
			return frame_set
		else:
			raise FrameNotFoundException("A valid frame set could not be read from the provided FITS header.")

	@property
	def baseFrame(self):
		'''
		Return the base frame.
		'''
		return ASTFrame(ast_object=self.astObject.getframe(Ast.BASE))
	
	@property
	def currentFrame(self):
		''' Returns the current frame. '''
		return ASTFrame(ast_object=self.astObject.getframe(Ast.CURRENT))
	
	def removeCurrentFrame(self):
		''' Remove the current frame from the frame set. '''
		# todo: check if the frame is actually part of the frame set
		# Note: The "removeframe" documenation is a little unclear on how to remove
		#       a frame that is not current.
		# Better to make this "removeFrame(frame=...) or by frame name, or add a new method.
		self.astObject.removeframe(iframe=Ast.CURRENT)
	
	def centerCoordinates(self):
		''' Returns the coordinates at the center of the frame. '''
		raise Exception("Useful, not sure how to do this.")
	
	def addToBaseFrame(self, frame=None):
		'''
		Add a new frame to this frame set's base frame.
		
		'''
		
		# Ref: http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211ss5.html#xref_astAddFrame
		# addframe( iframe, map, frame )
		
		if frame is None:
			raise Exception(f"frame must be provided to 'addFrame'")
			
		iframe = Ast.BASE # add to base frame
		map = None
		self.astObject.addframe(iframe, map, frame)
	
	@property
	def mapping(self):
		'''
		Return an object that maps points from the base frame to the current frame of this frame set.
		'''
		return ASTMapping(ast_object=self.astObject.getmapping()) # default is from base to current
			
	def convert(self, x_coordinates=None, y_coordinates=None):
		'''
		
		'''
		#TODO validate inputs for error reporting
		assert len(x_coordinates) == len(y_coordinates), "Coordinate arrays of unequal lengths."
		forward = True # convert from frame1 to frame2
		#in_coords = np.array([x_coordinates, y_coordinates]).T
		in_coords = np.array(list(zip(x_coordinates, y_coordinates))).T
		print(in_coords)
		forward = True
		return self.astObject.tran(in_coords, forward)
	
	
	#def frame(self, name=None):
	#	'''
	#	Extract a frame from the frame set with the given name.
	#	'''
	#	self.

