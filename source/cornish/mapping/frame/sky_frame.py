#/usr/bin/env python

import starlink.Ast as Ast

from .frame import ASTFrame

__all__ = ['ASTSkyFrame', 'ASTICRSFrame']

# this is a list of sky coordinate systems supported by AST,
# see: http://starlink.eao.hawaii.edu/docs/sun211.htx/sun211ss424.html
sky_systems = ["ICRS", "J2000", "AZEL", "ECLIPTIC", "FK4", "FK4-NO-E",
			   "FK4_NO_E", "FK5", "EQUATORIAL",
			   "GALACTIC", "GAPPT", "GEOCENTRIC", "APPARENT",
			   "HELIOECLIPTIC", "SUPERGALACTIC"]

class ASTSkyFrame(ASTFrame):
	'''

	A SkyFrame is a specialised form of Frame which describes celestial longitude/latitude coordinate systems.

	Systems available in AST:

	``ICRS``, ``J2000``, ``AZEL``, ``ECLIPTIC``, ``FK4``, ``FK4-NO-E``,
	``FK4_NO_E``, ``FK5``, ``EQUATORIAL``,
	``GALACTIC``, ``GAPPT``, ``GEOCENTRIC``, ``APPARENT``,
	``HELIOECLIPTIC``, ``SUPERGALACTIC``

	:param ast_object: an existing :class:`starlink.Ast.SkyFrame` object
	:param equinox: frame equinox, default value ``2000.0``
	:param system: coordinate system used to describe positions within the domain, see `AST System documentation <http://starlink.eao.hawaii.edu/docs/sun211.htx/sun211ss424.html>`_, default value = ``ICRS``
	:param epoch: epoch of the mean equinox as a string value, e.g. ``J2000.0``, ``B1950.0``, default = ``2000.0``
	'''
	def __init__(self, ast_object:Ast.SkyFrame=None, equinox:str=None, system:str=None, epoch:str=None):

		#if all([naxes, ast_frame]):
		#	raise Exception("The number of axes (naxes) argument cannot be specified with a provided ast_frame.")

		# TODO: if ast_frame provided, check it is a sky frame (see below)

		if ast_object and any([equinox, system, epoch]):
			raise ValueError("If 'ast_object' is provided, none of the other parameters ('equinox', 'system', 'epoch') can be specified.")

		if ast_object:
			if ast_object.isaskyframe():
				super().__init__(ast_object=ast_object)
			else:
				raise ValueError(f"The provided 'ast_object' value is not an Ast.SkyFrame (got '{type(ast_object)}').")
		else:
			self.astObject = Ast.SkyFrame()

		if system:
			if system.upper() in sky_systems:
				self.system = system
			else:
				raise ValueError(f"The provided system must be one of: [{sky_systems}].")
		if equinox:
			self.equinox = equinox
		if epoch:
			self.epoch = epoch

	@classmethod
	def fromFITSHeader(cls, header) -> "ASTSkyFrame":
		'''
		Returns a SkyFrame built from the WCS in the provided header, if found.
		Creates a sky frame from the provided FITS header.
		:raises: exc.NoWCSFound: if no sky frame WCS found
		'''
		frame = ASTFrameSet.fromFITSHeader(fits_header=header).currentFrame # -> ASTFrame

		# double check it's the right frame
		if frame.isSkyFrame:
			return frame
		else:
			# .. todo:: could perform a more rigorous to find a sky frame
			raise exc.NoWCSFoumd("A WCS corresponding to a sky frame was not be found.")

	@property
	def equinox(self):
		'''

		.. todo:: how to evaluate a valid equinox string?
		'''
		return self.astObject.get("Equinox")

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
		return float(self.astObject.get("Epoch"))

	@epoch.setter
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

# .. todo:: make this a factory class
class ASTICRSFrame(ASTSkyFrame):
	'''
	Factory class that returns an :class:`ASTSkyFrame` automatically set to ``System=ICRS``, ``equinox=2000.0``, ``epoch=2000.0``.
	'''
	def __init__(self, equinox:str="2000.0", epoch:str="2000.0") -> ASTSkyFrame:

		ast_object = Ast.SkyFrame()
		super().__init__(ast_object=ast_object)

		self.system = "ICRS"
		self.equinox = equinox
		self.epoch = epoch





