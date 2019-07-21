#/usr/bin/env python

from astropy.units import u
from .frame import ASTFrame
import starlink.Ast as Ast

class ASTSkyFrame(ASTFrame):
	'''
	
	self.astObject is of type starlink.Ast.SkyFrame.
	'''
	def __init__(self, ast_frame=None, equinox=None, system=None):

		#if all([naxes, ast_frame]):
		#	raise Exception("The number of axes (naxes) argument cannot be specified with a provided ast_frame.")
		
		if ast_frame is None:
			self.astObject = Ast.SkyFrame()
		else:
			self.astObject = ast_frame

		if system:
			self.system = system
		if equinox:
			self.equinox = equinox
		
	@property
	def equinox(self):
		'''
		
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






