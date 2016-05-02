#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

from .ast_object import ASTObject

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
	pass

