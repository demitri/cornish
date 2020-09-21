#!/usr/bin/env python

#from typing import Union
import starlink.Ast as Ast
from ..ast_object import ASTObject

__all__ = ['ASTChannel']

class ASTChannel(ASTObject):
	'''
	An ASTChannel is a class that handles data input/output into and from different data formats.
	
	ASTChannel objects encapsulate data serialization, i.e. they are used to read data from some
	format into an AST object, or can write an AST object to a file (e.g. a FITS file).
	In AST parlance, the data comes from a "source" and is written to a "sink".
	
	:param ast_object: an existing `starlink.Ast.Channel` object
	'''
	def __init__(self, ast_object=None):
		
		if ast_object:
			if isinstance(ast_object, Ast.Channel):
				super().__init__(ast_object=ast_object)
			else:
				raise ValueError(f"The parameter 'ast_object' must be an existing starlink.Ast.Channel object; got '{type(ast_object)}'.")
		
		raise NotImplementedError()
		# under construction...
	
