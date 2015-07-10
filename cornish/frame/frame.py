#/usr/bin/env python

from astropy.units import u
from ..ast_object import ASTObject

class Frame(ASTObject):
	'''
	A Frame is a representation of a coordinate system, e.g. Cartesian, RA/dec.
	It contains information about the labels which appear on the axes, the axis units,
	a title, knowledge of how to format the coordinate values on each axis, etc.
	'''
	def __init__(self):
		self.title = None
		self._number_of_axes = 0
		self.axis_labels = list() # a list of labels indexed by axis number, first=0
		self.axis_units = list() # a list of axis units indexed by axis number, first=0
		
	@property
	def number_axes(self):
		pass
		
class CompoundFrame(ASTObject):
	'''
	A compound frame is the merging of two existing `Frame`s.
	For example, a `CompoundFrame` could have celestial coordinates on two axes
	and an unrelated coordinate (wavelength, perhaps) on a third.
	Knowledge of the relationships between the axes is preserved internally
	by the process of constructing the `CompoundFrame` which represents them.
	'''
	def __init__(self):
		pass
	
class FrameSet(ASTObject):
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
	def __init__(self):
		self.base_frame = None
		self.current_frame = None