#/usr/bin/env python
#from __future__ import (absolute_import, division, print_function, unicode_literals)

import numpy as np

import starlink
import starlink.Ast as Ast

from ... import ASTObject
from .frame import ASTFrame
from ..mapping import ASTMapping
#from ...channel import ASTFITSChannel

class ASTFrameSet(ASTFrame):
	'''
	A sets of inter-related coordinate systems made up of existing `Mapping`s and `Frame`s.
	A FrameSet may be extended by adding a new Frame and associated Mapping.
	
	A `FrameSet` must have a "base" frame which represents the "native" coordinate system
	(for example, the pixel coordinates of an image). Similarly, one Frame is termed the
	current Frame and represents the "currently-selected" coordinates. It might typically
	be a celestial or spectral coordinate system and would be used during interactions
	with a user, as when plotting axes on a graph or producing a table of results.
	Other Frames within the FrameSet represent a library of alternative coordinate systems
	which a software user can select by making them current.
	
	
	'''
	def __init__(self, ast_object=None, base_frame=None, ast_frame_set=None, fits_header=None):
		'''
		Create a new AST frame set.
		Object can be created from an starlink.Ast.FrameSet "primitive"
		(e.g. returned by another object).
		
		self.astObject is of type starlink.Ast.FrameSet
		'''
		# TODO: properly handle when ast_object is set
		#super().__init__(ast_object=ast_object)
		if ast_object:
			assert all([x is None for x in [base_frame, ast_frame_set, fits_header]]), "too many parameters set" 
			self.astObject = ast_object
			return
		
		if all([base_frame, ast_frame_set, fits_header]):
			raise Exception("Only 'base_frame' or 'ast_frame_set' may be used to initialize this object.")
		elif not any([base_frame, ast_frame_set, fits_header]):
			raise Exception("This object must be initialized with either a base frame or AST frame set object")
		
		# todo: fully test for invalid combinations of parameters
		
		if fits_header:
			from cornish import ASTFITSChannel # TODO: reorg later
			fits_channel = ASTFITSChannel(header=fits_header)
			self.astObject = fits_channel.frameSet.astObject
			return
		
		if ast_frame_set is not None:
			if isinstance(ast_frame_set, starlink.Ast.FrameSet):
				self.astObject = ast_frame_set
			else:
				Exception("ASTFrameSet: Unhandled ast_frame_set type ('{0}')".format(ast_frame_set))
		else:
			if isinstance(base_frame, starlink.Ast.Frame):
				self.astObject = Ast.FrameSet(frame=base_frame, options=None)
			else:
				Exception("ASTFrameSet: Unhandled base_frame type ('{0}')".format(base_frame))
	
		
	
	def setFromFrames(frame1:ASTFrame, frame2:ASTFrame):
		'''
		Static method to create a frame set from two existing frames.
		'''
		frame_set = Ast.Frame.convert(frame1.astObject, frame2.astObject)
		if frame_set is None:
			# I don't know if this is actually the failure mode
			raise Exception("An ASTFrameSet could not be created since the mapping between the two provided frames could not be determined (conversion may not be possible).")
		else:
			return ASTFrameSet(ast_frame_set=frame_set)
	setFromFrames = staticmethod(setFromFrames)
	
	@property
	def baseFrame(self):
		'''
		Return the base frame.
		'''
		return ASTFrame(ast_frame=self.astObject.getframe(Ast.BASE))
	
	@property
	def currentFrame(self):
		''' Returns the current frame. '''
		return ASTFrame(ast_frame=self.astObject.getframe(Ast.CURRENT))
	
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
#		in_coords = np.array([x_coordinates, y_coordinates]).T
		in_coords = np.array(list(zip(x_coordinates, y_coordinates))).T
		print(in_coords)
		forward = True
		return self.astObject.tran(in_coords, forward)
	
	
	#def frame(self, name=None):
	#	'''
	#	Extract a frame from the frame set with the given name.
	#	'''
	#	self.

