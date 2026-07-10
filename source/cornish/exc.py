
'''
Custom exceptions for the cornish package.

All cornish-specific exceptions inherit from :class:`CornishError`, so
``except CornishError`` catches everything this package raises on its own
behalf. Plain :class:`TypeError`/:class:`ValueError` are still used per the
package-wide policy: an argument of the wrong *type* raises ``TypeError``; an
argument of the correct type with a bad *value* raises ``ValueError`` or a
more specific :class:`CornishError` subclass.
'''

__all__ = ['CornishError', 'FrameNotFoundException', 'NotA2DRegion', 'NotASkyRegion',
           'CoordinateSystemsCouldNotBeMapped', 'NoWCSFound', 'IncompleteHeader',
           'UnsupportedASTClass', 'SerializationNotPossible', 'RegionNotConnected']

class CornishError(Exception):
	''' Root of the cornish exception hierarchy. '''
	pass

class FrameNotFoundException(CornishError):
	''' Raised when a required frame (e.g. a WCS frame set from a FITS header) cannot be found. '''
	pass

class NotA2DRegion(CornishError):
	''' Raised when an operation requires a 2D region (e.g. a bounding circle). '''
	pass

class NotASkyRegion(CornishError):
	''' Raised when a sky region is required, e.g. converting a region to a MOC. '''
	pass

class CoordinateSystemsCouldNotBeMapped(CornishError):
	''' Raised when AST cannot determine a mapping between two coordinate systems. '''
	pass

class NoWCSFound(CornishError):
	''' Raised when no usable WCS could be read, e.g. from a FITS header. '''
	pass

class IncompleteHeader(CornishError):
	''' Raised when a FITS header lacks cards required for the requested operation (e.g. NAXIS1/NAXIS2). '''
	pass

class UnsupportedASTClass(CornishError):
	''' Raised when an AST object of a class not (yet) wrapped by cornish is encountered. '''
	pass

class SerializationNotPossible(CornishError):
	''' Raised when an object cannot be represented in the requested serialization format. '''
	pass

class RegionNotConnected(CornishError):
	''' Raised when an operation requires a simply-connected region but the region has multiple patches. '''
	pass

NoWCSFoumd = NoWCSFound # deprecated alias (original misspelling), kept for compatibility
