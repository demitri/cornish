#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

from ..ast_object import ASTObject

class ASTChannel(ASTObject):
	'''
	An ASTChannel is an class that handles data input/output into and from different data formats.
	
	ASTChannel objects encapsulate data serialization, i.e. they are used to read data from some
	format into an AST object, or can write an AST object to a file (e.g. a FITS file).
	In AST parlance, the data comes from a "source" and is written to a "sink".
	'''
	def __init__(self):
		# defines internal AST object
		super(ASTChannel, self).__init__()
		# super().__init__() # Python 3
	
	#def __repr__(self):
	#	'''
	#	
	#	'''
		