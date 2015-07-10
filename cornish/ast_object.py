#!/usr/bin/env python

class ASTObject(object):
	'''
	This is the root class for all AST objects.
	'''
	def __init__(self):
		pass
	
	def ast_description(self):
		return "This is an AST object (override the 'ast_description' method for a more descriptive output)."
