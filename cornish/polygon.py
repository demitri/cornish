#!/usr/bin/env python

from ..ast_object import ASTObject
import starlink.Ast as Ast

class Polygon(ASTObject):
	'''
	
	self.astObject is of type starlink.Ast.Polygon.
	'''
	def __init__(self):
		self.vertices = None
		self.title = None
		self.domain = None
		self.epoch = None
		self.system = None
		self.axes_labels = ["Axis 1", "Axis 2"]
		self.axes_units = ["pixel", "pixel"]
		self.frame = None