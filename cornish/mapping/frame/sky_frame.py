#/usr/bin/env python

from astropy.units import u
from .frame import ASTFrame
import starlink.Ast as Ast

# this is a list of sky coordinate systems supported by AST,
# see: http://starlink.eao.hawaii.edu/docs/sun211.htx/sun211ss424.html
sky_systems = ["ICRS", "J2000", "AZEL", "ECLIPTIC", "FK4", "FK4-NO-E",
			   "FK4_NO_E", "FK5", "EQUATORIAL",
			   "GALACTIC", "GAPPT", "GEOCENTRIC", "APPARENT",
			   "HELIOECLIPTIC", "SUPERGALACTIC"]

equinoxes

class ASTSkyFrame(ASTFrame):
	'''
	
	self.astObject is of type starlink.Ast.SkyFrame.
	
	Systems available in AST:
	    "ICRS", "J2000", "AZEL", "ECLIPTIC", "FK4", "FK4-NO-E",
	    "FK4_NO_E", "FK5", "EQUATORIAL",
	    "GALACTIC", "GAPPT", "GEOCENTRIC", "APPARENT",
	    "HELIOECLIPTIC", "SUPERGALACTIC"
	
	:param ast_object:
	:param equinox: 
	:param system: coordinate system used to describe positions within the domain, see `AST System <http://starlink.eao.hawaii.edu/docs/sun211.htx/sun211ss424.html>`
	:param epoch: epoch of the mean equinox as a string value, e.g. "J2000.0", "B1950.0"
	'''
	def __init__(self, ast_object:starlink.Ast.SkyFrame=None, equinox:str=None, system:str=None, epoch:str=None):

		#if all([naxes, ast_frame]):
		#	raise Exception("The number of axes (naxes) argument cannot be specified with a provided ast_frame.")
		
		# TODO: if ast_frame provided, check it is a sky frame (see below)
		
		if ast_object is None:
			self.astObject = Ast.SkyFrame()
		else:
			if ast_object.isaskyframe():
				self.astObject = ast_frame
			else:
				raise ValueError(f"The provided 'ast_object' value is not an Ast.SkyFrame (got '{type(ast_object)}').")
		if system:
			if system in sky_systems:
				self.system = system
			else:
				raise ValueError(f"The provided system must be one of: [{sky_systems}].")
		if equinox:
			self.equinox = equinox
		if epoch:
			self.epoch = epoch
		
	@property
	def equinox(self):
		'''
		
		..todo :: how to evaluate a valid equinox string?
		'''
		self.astObject.get("Equinox") 
		
	@equinox.setter
	def equinox(self, equinox=None):
		'''	Set the equinox for the frame. '''
		if equinox is None:
			raise Exception("An 'equinox' parameter must be specified.")
		elif not isinstance(equinox, str):
			raise Exception("A string value was expected for 'equinox' (got '{0}').".format(type(equinox)))
		self.astObject.set(f"Equinox={equinox}")

	@property
	def epoch(self):
		'''
		
		'''
		self.astObject.get("Epoch") 
		
	@equinox.setter
	def epoch(self, epoch=None):
		'''	Set the epoch for the frame. '''
		if epoch is None:
			raise Exception("An 'epoch' parameter must be specified.")
		else:
			if isinstance(epoch, str):
				try:
					float(epoch)
				except ValueError as e:
					raise ValueError("'epoch' much be a numeric value (or a string that can be converted to a numeric value")
		self.astObject.set(f"Epoch={epoch}")

class ASTICRSFrame(ASTSkyFrame):
	
	def __init__(self, ast_frame=None): # this should be --> ast_object=None):
		#super().__init__(ast_object=ast_object)
		super().__init__(ast_frame=ast_frame)
		self.system = "ICRS"




