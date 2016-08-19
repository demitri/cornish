#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

class ASTObject(object):
	'''
	This is the root class for all AST objects.
	
	Every object subclassed from this abstract superclass is a wrapper for a starlink.Ast object.
	This object is stored in the property "astObject".
	'''
	def __init__(self):
		self.astObject = None
		pass
	
	def ast_description(self):
		return "This is an AST object (override the 'ast_description' method for a more descriptive output)."
