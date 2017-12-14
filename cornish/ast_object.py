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

	@property
	def id(self):
		'''
		String which may be used to identify this object.
		
		NOTE: Unlike most other attributes, the value of the ID attribute is not transferred when
		an Object is copied. Instead, its value is undefined (and therefore defaults to an empty string)
		in any copy. However, it is retained in any external representation of an Object produced by
		the astWrite function.
		
		Not sure how to handle the above in this class.
		
		@returns string identifier that can be used to uniquely identify this object
		'''
		return self.astObject.get("ID")
		
