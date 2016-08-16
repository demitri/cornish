#/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

from astropy.units import u
import starlink.Ast as Ast

from ... import ASTObject
from ..mapping import ASTMapping

class ASTFrame(ASTMapping):
	'''
	A Frame is a representation of a coordinate system, e.g. Cartesian, RA/dec.
	It contains information about the labels which appear on the axes, the axis units,
	a title, knowledge of how to format the coordinate values on each axis, etc.
	
	List and description of starlink.Ast.Frame attributes in documentation: Section 7.5.
	Ref:
	http://www.starlink.rl.ac.uk/docs/sun95.htx/sun95se27.html
	http://www.strw.leidenuniv.nl/docs/starlink/sun210.htx/node71.html
	'''
	def __init__(self, naxes=None, ast_frame=None):
		#self.title = None
		#self._number_of_axes = 0
		#self.axis_labels = list() # a list of labels indexed by axis number, first=0
		#self.axis_units = list() # a list of axis units indexed by axis number, first=0
		
		if all([naxes, ast_frame]):
			raise Exception("The number of axes (naxes) argument cannot be specified with a provided ast_frame.")
		
		if ast_frame is None:
			self.astFrame = Ast.Frame(naxes=naxes)
		else:
			self.astFrame = ast_frame
		
	@property
	def naxes(self):
		''' Returns the number of axes for the frame. '''
		return int(self.astFrame.get("Naxes"))
		
	@property
	def title(self):
		''' Returns the frame title, a string describing the coordinate system which the frame represents. '''
		return self.astFrame.get("Title")
	
	@title.setter
	def title(self, newTitle):
		self.astFrame.set("Title={0}".format(newTitle))
	
	def label(self, axis=None):
		''' Return the label for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.astFrame.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		return self.astFrame.get("Label({0})".format(axis))

	def setLabelForAxis(axis=None, label=None):
		''' Set the label for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.astFrame.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		elif label is None:
			raise Exception("A new label must be specified.")
		self.astFrame.set("Label({0}={1}".format(axis, newLabel))
		
	def unit(self, axis=None):
		''' Return the unit for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.astFrame.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		self.astFrame.get("Unit({0})".format(axis))

	def setUnitForAxis(self, axis=None, newUnit=None):
		''' Set the unit as a string value for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.astFrame.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		elif newUnit is None:
			raise Exception("A new unit must be specified.")
		self.astFrame.set("Unit({0})={1}".format(axis, newUnit))
	
	@property
	def system(self):
		'''
		String which identifies the coordinate system represented by the Frame.
		The system is "Cartesian" by default, but can have other values for subclasses of Frame,
		e.g. "FK4", "Galactic".
		'''
		return self.astFrame.get("System")
	
	@system.setter
	def system(self, system=None):
		''' Set the a label for the frame's system. '''
		if system is None:
			raise Exception("A 'system' parameter must be specified.")
		elif not isinstance(system, str):
			raise Exception("A string value was expected for 'system' (got '{0}').".format(type(system)))
		self.astFrame.set("System={0}".format(system))
		
	@property
	def domain(self):
		'''
		The physical domain of the coordinate system (string value).
		The Domain attribute also controls how Frames align with each other.
		If the Domain value in a Frame is set, then only Frames with the same Domain value can be aligned with it.
		
		Example values: GRID, FRACTION, PIXEL, AXIS, SKY, SPECTRUM, CURPIC, NDC, BASEPIC, GRAPHICS
		
		Frames created by the user (for instance, using WCSADD) can have any Domain value, but the standard
		domain names should be avoided unless the standard meanings are appropriate for the Frame being created.
		'''
		return self.astFrame.get("Domain")
		
	@domain.setter
	def domain(self, newDomain=None):
		if newDomain is None or isinstance(newDomain, str) == False:
			raise Exception("The domain value must be set to a string.")
		self.astFrame.set("Domain={0}".format(newDomain))
	
class ASTCompoundFrame(ASTObject):
	'''
	A compound frame is the merging of two existing `Frame`s.
	For example, a `CompoundFrame` could have celestial coordinates on two axes
	and an unrelated coordinate (wavelength, perhaps) on a third.
	Knowledge of the relationships between the axes is preserved internally
	by the process of constructing the `CompoundFrame` which represents them.
	'''
	def __init__(self):
		pass
	
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
		return ASTFrame(ast_frame=self.astFrameSet.getFrame(Ast.CURRENT))
	
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

