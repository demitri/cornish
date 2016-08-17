from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink.Ast as Ast

from .region import ASTRegion
from ..mapping import ASTMapping
from ..mapping import ASTFrame

__all__ = ["ASTCircle"]

class ASTCircle(ASTRegion):
	
	def __init__(self):
		raise Exception("ASTCircle not yet implemented.")
	