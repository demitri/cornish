
'''
Custom exceptions for the cornish package.
'''

class FrameNotFoundException(Exception):
	pass

class NotA2DRegion(Exception):
	pass

class CoordinateSystemsCouldNotBeMapped(Exception):
	pass

class NoWCSFound(Exception):
	pass

class NotASkyRegion(Exception):
	''' Raised when a sky region is required, e.g. converting a region to a MOC. '''
	pass

NoWCSFoumd = NoWCSFound # deprecated alias (original misspelling), kept for compatibility
