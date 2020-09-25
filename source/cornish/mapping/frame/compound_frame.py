
from ... import ASTObject
		
class ASTCompoundFrame(ASTObject):
	'''
	A compound frame is the merging of two existing :class:`ASTFrame` objects.

	For example, a compound frame could have celestial coordinates on two axes
	and an unrelated coordinate (wavelength, perhaps) on a third.
	Knowledge of the relationships between the axes is preserved internally
	by the process of constructing the frames which represents them.
	'''
	def __init__(self, ast_object=None):
		raise NotImplementedError()
	
