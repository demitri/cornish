#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

from ..ast_object import ASTObject

_astropy_available = True
_fitsio_available = True
try:
	import astropy
except ImportError:
	_astropy_available = False
try:
	import fitsio
except ImportError:
	_fitsio_available = False

class Channel(ASTObject):
	'''
	A Channel is an class that handles data input/output into and from different data formats.
	
	Channel objects encapsulate data serialization, i.e. they are used to read data from some
	format into an AST object, or can write an AST object to a file (e.g. a FITS file).
	In AST parlance, the data comes from a "source" and is written to a "sink".
	'''
	pass

