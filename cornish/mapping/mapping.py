#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

from ..ast_object import ASTObject
import starlink.Ast as Ast

class ASTMapping(ASTObject):
	'''
	
	self.astObject is of type starlink.Ast.Mapping.
	'''
	def __init__(self):
		pass
	
	@property
	def number_of_input_coordinates(self):
		'''
		Number of dimensions of the space in which the Mapping’s input points reside.
		This property gives the number of coordinate values required to specify an input point for a Mapping.
		
		@returns number of dimensions described by this mapper
		'''
		return self.astObject.get("Nin")
