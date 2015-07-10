#!/usr/bin/env python

from ..ast_object import ASTObject

class Polygon(ASTObject):
	
	def __init__(self):
		self.vertices = None
		self.title = None
		self.domain = None
		self.epoch = None
		self.system = None
		self.axes_labels = ["Axis 1", "Axis 2"]
		self.axes_units = ["pixel", "pixel"]
		self.frame = None