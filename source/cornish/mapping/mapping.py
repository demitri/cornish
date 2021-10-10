#!/usr/bin/env python

from typing import Iterable

import numpy as np

from ..ast_object import ASTObject

class ASTMapping(ASTObject):
	'''


	:param ast_object: an existing :class:`starlink.Ast.Mapping` object
	'''
	def __init__(self, ast_object=None):
		super().__init__(ast_object=ast_object)

	@property
	def numberOfInputCoordinates(self):
		'''
		Number of dimensions of the space in which the Mapping’s input points reside.
		This property gives the number of coordinate values required to specify an input point for a Mapping.

		:returns: number of input dimensions described by this mapper
		'''
		return self.astObject.get("Nin")

	@property
	def numberOfOutputCoordinates(self):
		'''
		Number of dimensions of the space in which the Mapping’s output points reside.
		This property gives the number of coordinate values required to specify an output point for a Mapping.

		:returns: number of output dimensions described by this mapper
		'''
		return self.astObject.get("Nout")

	@property
	def isLinear(self) -> bool:
		'''
		Returns True if the mapping is linear
		'''
		return self.astObject.get("IsLinear") == 1

	@property
	def isSimple(self) -> bool:
		'''
		Returns True if the mapping has been simplified.
		'''
		return self.astObject.get("IsSimple") == 1

	def inverseMapping(self):
		'''
		Returns a new mapping object that is the inverse of this mapping.

		For example, if the forward transformation of this mapping is pixel to sky,
		then the forward transformation of the returned mapping
		will be sky to pixel.
		'''
		inv_mapping = self.astObject.copy() # creates a new mapping
		return ASTMapping(ast_object=inv_mapping.invert())

	# TODO, test "out" param, out:Iterable=None): # ? use forward=T/F
	def transform(self, points:Iterable=None) -> np.ndarray:
		'''
		Transform the coordinates of a set of points provided according the mapping defined by this object.

		:param in: input list of coordinates as numpy.ndarray,
				   any iterable list accepted
				   2-dimensional array with shape (nin,npoint)
		:param out: output coordinates

		Format of points:

		.. code-block::

		    [ [ values on axis 1 ], [ values on axis 2 ], ... ]

		e.g. sky to pixel:

		.. code-block::

		    [ [ra1, ra2, ...], [dec1, dec2, ...] ]

		'''
		#raise NotImplementedError("This needs to be tested.") # use this instead of "convert" in ASTFrameSet?
		return self.astObject.tran(points)

