#/usr/bin/env python

from astropy.units import u
import starlink.Ast as Ast

from ... import ASTObject
from ..mapping import ASTMapping

__all__ = ['ASTFrame', 'ASTCompoundFrame']

class ASTFrame(ASTMapping):
	'''
	A Frame is a representation of a coordinate system, e.g. Cartesian, RA/dec.
	It contains information about the labels which appear on the axes, the axis units,
	a title, knowledge of how to format the coordinate values on each axis, etc.
	
	List and description of starlink.Ast.Frame attributes in documentation: Section 7.5.
	Ref:
	http://www.starlink.rl.ac.uk/docs/sun95.htx/sun95se27.html
	http://www.strw.leidenuniv.nl/docs/starlink/sun210.htx/node71.html
	
	self.astObject is of type starlink.Ast.Frame.
	'''
	def __init__(self, ast_object:Ast.Frame=None, naxes:int=None):

		if all([ast_object, naxes]):
			raise ValueError("Cannot initialize ASTFrame with both 'ast_object' and 'naxes' set.")
		
		if all([x is None for x in [ast_object, naxes]]):
			raise ValueError("Either 'ast_object' or 'naxes' must be specified to create an ASTFrame.")
		
		if ast_object:
			if isinstance(ast_object, Ast.Frame):
				super().__init__(ast_object=ast_object)
			else:
				raise ValueError(f"The provided 'ast_object' is not an Ast.Frame object (got class '{type(ast_object)}').")
		else:
			self.astObject = Ast.Frame(naxes)

		#self.axis_labels = list() # a list of labels indexed by axis number, first=0
		#self.axis_units = list() # a list of axis units indexed by axis number, first=0
# 		
# 		if ast_frame is None:
# 			self.astObject = Ast.Frame(naxes)
# 		else:
# 			self.astObject = ast_frame
	
	@staticmethod
	def frameFromAstObject(ast_object:Ast.Frame=None):
		'''
		Factory method that returns the appropriate Cornish frame object (e.g. ASTSkyFrame) for a given frame.
		
		:param ast_object: an :py:class:`Ast.Frame` object
		'''
		if ast_object is None:
			raise ValueError("An ast_object must be specified.")
		elif not isinstance(ast_object, Ast.Frame):
			raise Exception(f"Expected 'ast_object' to be some kind of Ast frame, but got {type(ast_object)}")
		
		# the order might need to be rearranged depending on the object hierarchy
		
		if isinstance(ast_object, Ast.SkyFrame):
			return ASTSkyFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.SpecFrame):
			return ASTSpecFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.TimeFrame):
			return ASTTimeFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.CmpFrame):
			return ASTCompositeFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.FluxFrame):
			return ASTFluxFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.SpecFluxFrame):
			return ASTSpecFluxFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.TimeFrame):
			return ASTTimeFrame(ast_object=ast_object)
		else:
			return ASTFrame(ast_object=ast_object)
	
	@property
	def naxes(self):
		''' Returns the number of axes for the frame. '''
		return int(self.astObject.get("Naxes"))
		
	@property
	def title(self):
		''' Returns the frame title, a string describing the coordinate system which the frame represents. '''
		return self.astObject.get("Title")
	
	@title.setter
	def title(self, newTitle):
		self.astObject.set("Title={0}".format(newTitle))
	
	def label(self, axis=None):
		''' Return the label for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		return self.astObject.get("Label({0})".format(axis))

	def setLabelForAxis(axis=None, label=None):
		''' Set the label for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		elif label is None:
			raise Exception("A new label must be specified.")
		self.astObject.set("Label({0}={1}".format(axis, newLabel))
		
	def unit(self, axis=None):
		''' Return the unit for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		self.astObject.get("Unit({0})".format(axis))

	def setUnitForAxis(self, axis=None, unit=None):
		''' Set the unit as a string value for the specified axis. '''
		if axis is None:
			raise Exception("An axis number must be specified.")
		elif not isinstance(axis, int):
			raise Exception("The parameter 'axis' must be an integer (a '{0}' was provided).".format(type(axis)))
		elif axis > self.naxes:
			raise Exception("The axis provided ({0}) is larger than the number of axes ({1}).".format(axis, self.naxes))
		elif unit is None:
			raise Exception("A new unit must be specified.")
		self.astObject.set("Unit({0})={1}".format(axis, unit))
	
	@property
	def system(self):
		'''
		String which identifies the coordinate system represented by the Frame.
		The system is "Cartesian" by default, but can have other values for subclasses of Frame,
		e.g. "FK4", "Galactic".
		'''
		return self.astObject.get("System")
	
	@system.setter
	def system(self, system=None):
		''' Set the a label for the frame's system. '''
		if system is None:
			raise Exception("A 'system' parameter must be specified.")
		elif not isinstance(system, str):
			raise Exception("A string value was expected for 'system' (got '{0}').".format(type(system)))
		self.astObject.set("System={0}".format(system))
		
	@property
	def isSkyFrame(self):
		'''
		Returns 'True' if this is a SkyFrame, 'False' otherwise.
		'''
		return self.astObject.isaskyframe()
	
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
		return self.astObject.get("Domain")
		
	@domain.setter
	def domain(self, newDomain=None):
		if newDomain is None or isinstance(newDomain, str) == False:
			raise Exception("The domain value must be set to a string.")
		self.astObject.set("Domain={0}".format(newDomain))
		
	def distance(self, point1, point2):
		'''
		Distance between two points in this frame.
		
		Parameters
		----------
		point1 : list
			A two element list/tuple/Numpy array of the first point coordinates.
		point2 : list
			A two element list/tuple/Numpy array of the second point coordinates.
		
		Returns
		-------
		float, distance between two provided points
		'''
		return self.astObject.distance(point1, point2)
		
	def framesetWithMappingTo(self, template_frame=None): # maybe think of better name?!
		'''
		Searches 
		:param template_frame:
		:returns: ... , None if no mapping can be found
		'''
		from .frame_set import ASTFrameSet
		if template_frame is None:
			raise ValueError("A template frame must be provided.")
		return ASTFrameSet(ast_object=self.astObject.findframe(template_frame.astObject))
		
class ASTCompoundFrame(ASTObject):
	'''
	A compound frame is the merging of two existing `Frame`s.
	For example, a `CompoundFrame` could have celestial coordinates on two axes
	and an unrelated coordinate (wavelength, perhaps) on a third.
	Knowledge of the relationships between the axes is preserved internally
	by the process of constructing the `CompoundFrame` which represents them.
	'''
	def __init__(self, ast_object=None):
		raise NotImplementedError()
	
