#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

from .channel import ASTChannel

_astropy_available = True
_fitsio_available = True
try:
	import astropy
	import astropy.io
except ImportError:
	_astropy_available = False
try:
	import fitsio
except ImportError:
	_fitsio_available = False

class ASTFITSHeader(object):
	
	def __init__(self, hdu):
		
		header = None
		if _astropy_available:
			if isinstance(fits_file, astropy.io.fits.HDU):
				header = hdu.header # -> astropy.io.fits.Header
		
		if header is None and isinstance(fits_file, fitsio.):
			header = hdu.read_header() # -> fitsio.fitslib.FITSHDR
		
		
	

