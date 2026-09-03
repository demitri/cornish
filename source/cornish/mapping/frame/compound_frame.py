
from ... import ASTObject

from .frame import ASTFrame
from ..._registry import register_frame

__all__ = ['ASTCompoundFrame']

@register_frame("CmpFrame")
class ASTCompoundFrame(ASTFrame):
	'''
	A compound frame is the merging of two existing :class:`ASTFrame` objects.

	For example, a compound frame could have celestial coordinates on two axes
	and an unrelated coordinate (wavelength, perhaps) on a third.
	Knowledge of the relationships between the axes is preserved internally
	by the process of constructing the frames which represents them.
	'''
	pass # no specific behaviour is yet defined for this subclass
	#def __init__(self, ast_object=None):
		
		
	
