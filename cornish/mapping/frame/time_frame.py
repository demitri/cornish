#/usr/bin/env python

from astropy.units import u
from ...ast_object import ASTObject
import starlink.Ast as Ast

__all__ = ['ASTTimeFrame']

class ASTTimeFrame(ASTObject):
	'''
	
	self.astObject is of type TimeFrame.
	'''
	def __init__(self, ast_object=None):
		raise NotImplementedError()