#!/usr/bin/env python

from typing import Iterable

from ..ast_object import ASTObject
import starlink.Ast as Ast

class ASTMapping(ASTObject):
	'''
	
	:param ast_object: is of type :py:class:``starlink.Ast.Mapping``.
	'''
	def __init__(self, ast_object=None):
		super().__init__(ast_object=ast_object)
	
	@property
	def number_of_input_coordinates(self):
		'''
		Number of dimensions of the space in which the Mapping’s input points reside.
		This property gives the number of coordinate values required to specify an input point for a Mapping.
		
		:returns: number of input dimensions described by this mapper
		'''
		return self.astObject.get("Nin")

	@property
	def number_of_output_coordinates(self):
		'''
		Number of dimensions of the space in which the Mapping’s output points reside.
		This property gives the number of coordinate values required to specify an output point for a Mapping.
		
		:returns: number of output dimensions described by this mapper
		'''
		return self.astObject.get("Nout")

	def inverseMapping(self):
		'''
		Returns a new mapping object that is the inverse of this mapping.
		
		For example, if the forward transformation of this mapping is pixel to sky,
		then the forward transformation of the returned mapping
		will be sky to pixel.
		'''
		inv_mapping = self.astObject.copy() # creates a new mapping
		return ASTMapping(ast_object=inv_mapping.invert())
	
	def transform(self, points:Iterable=None): # TODO, test "out" param, out:Iterable=None): # ? use forward=T/F
		'''
		Transform the coordinates of a set of points provided according the mapping defined by this object.
		
		:param in: input list of coordinates as numpy.ndarray,
				   any iterable list accepted
				   2-dimensional array with shape (nin,npoint)
		:param out:
		:returns: numpy.ndarray
		
		Format of points: [ [ values on axis 1 ], [ values on axis 2], ... ]
		e.g. sky to pixel: [ [ra1, ra2, ...], [dec1, dec2, ...] ]
		'''
		return self.astObject.tran(points)
		
		