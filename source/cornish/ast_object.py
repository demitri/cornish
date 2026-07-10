
import copy as _copy
from abc import ABCMeta

class ASTObject(metaclass=ABCMeta):
	'''
	This is the root class for all AST objects.

	Every object subclassed from this abstract superclass is a wrapper for a starlink.Ast object.
	This object is stored in the property "astObject".
	'''
	def __init__(self, ast_object=None):
		self.astObject = ast_object

	def ast_description(self):
		''' A string description of this object, customized for each subclass of :class:`ASTObject`. '''
		return "This is an AST object (override the 'ast_description' method for a more descriptive output)."

	@property
	def id(self) -> str:
		'''
		String which may be used to identify this object.

		NOTE: Unlike most other attributes, the value of the ID attribute is not transferred when
		an Object is copied. Instead, its value is undefined (and therefore defaults to an empty string)
		in any copy. However, it is retained in any external representation of an Object produced by
		the astWrite function.

		:returns: string identifier that can be used to uniquely identify this object
		'''
		return self.astObject.get("ID")

	def __repr__(self):
		# a one-liner; the full AST serialization (60+ lines) remains available as .astString
		try:
			ast_class = self.astObject.Class
		except AttributeError:
			# a repr must never raise: astObject may be None/unset on a
			# partially-constructed wrapper being inspected in a debugger
			ast_class = "no AST object attached"
		return f"<{type(self).__module__}.{type(self).__name__} ({ast_class}) at {hex(id(self))}>"

	def __reduce__(self):
		# Pickle support (DECISIONS D6/D16): persist the AST-native text dump plus
		# the Python-side attributes; `reconstruct` rehydrates via Ast.Channel and
		# cls.__new__ (deliberately bypassing __init__ — constructor validation is
		# for users, not for rehydration).
		from . import _pyast_bridge
		state = {k: v for k, v in self.__dict__.items() if k != "astObject"}
		return (_pyast_bridge.reconstruct, (type(self), self.astString), state)

	def copy(self):
		'''
		Return a new instance of this class wrapping a deep copy of the underlying AST object.

		Python-side attributes are carried over by reference (use :func:`copy.deepcopy`
		for a fully independent copy).
		'''
		cloned = type(self).__new__(type(self))
		cloned.__dict__.update({k: v for k, v in self.__dict__.items() if k != "astObject"})
		cloned.astObject = self.astObject.copy()
		return cloned

	def __deepcopy__(self, memo):
		cloned = type(self).__new__(type(self))
		memo[id(self)] = cloned
		for key, value in self.__dict__.items():
			if key == "astObject":
				continue
			setattr(cloned, key, _copy.deepcopy(value, memo))
		cloned.astObject = self.astObject.copy()
		return cloned

	@property
	def astString(self):
		'''
		Return the AST serialization of this object.
		'''
		return str(self.astObject)
