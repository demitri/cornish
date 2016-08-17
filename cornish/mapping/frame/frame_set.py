#/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink.Ast as Ast

from ... import ASTObject
from .frame import ASTFrame

class ASTFrameSet(ASTObject):
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
	def __init__(self, base_frame=None, ast_frame_set=None):
		'''
		Create a new AST frame set.
		Object can be created from an starlink.Ast.FrameSet "primitive"
		(e.g. returned by another object).
		'''
		if all([base_frame, ast_frame_set]):
			raise Exception("Only 'base_frame' or 'ast_frame_set' may be used to initialize this object.")
		elif not any([base_frame, ast_frame_set]):
			raise Exception("This object must be initialized with either a base frame or AST frame set object")
			
		if ast_frame_set:
			self.astFrameSet = ast_frame_set
		else:
			self.astFrameSet = Ast.FrameSet(frame=base_frame, options=None)
		
		#self.base_frame = None
		#self.current_frame = None
	
	def baseFrame(self):
		'''
		Return the base frame.
		'''
		return ASTFrame(ast_frame=self.astFrameSet.getframe(Ast.BASE))
	
	def currentFrame(self):
		''' Returns the current frame. '''
		return ASTFrame(ast_frame=self.astFrameSet.getframe(Ast.CURRENT))
	
	def removeCurrentFrame(self):
		''' Remove the current frame from the frame set. '''
		# todo: check if the frame is actually part of the frame set
		# Note: The "removeframe" documenation is a little unclear on how to remove
		#       a frame that is not current.
		# Better to make this "removeFrame(frame=...) or by frame name, or add a new method.
		self.astFrameSet.removeframe(iframe=Ast.CURRENT)
	
	def centerCoordinates(self):
		''' Returns the coordinates at the center of the frame. '''
		raise Exception("Useful, not sure how to do this.")
	
	#def frame(self, name=None):
	#	'''
	#	Extract a frame from the frame set with the given name.
	#	'''
	#	self.

