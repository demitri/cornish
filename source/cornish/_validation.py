
'''
Tiny shared parameter-validation helpers (private module).

The integer policy (decided once, rule-of-two, 2026-07-10): integer-valued
parameters accept anything that implements ``__index__`` — Python ints and
NumPy integer scalars alike — and loudly reject bools (which pass
``operator.index``) and floats (no silent truncation). Every integer seam in
the package routes through :func:`as_integer` so the policy cannot drift
site by site.
'''

import operator

__all__ = ['as_integer']

def as_integer(value, name:str) -> int:
	'''
	Coerce an integer-valued parameter per the package policy.

	:param value: the user-provided value
	:param name: the parameter name, for error messages
	:returns: the normalized Python int
	:raises TypeError: for bools and anything not implementing ``__index__`` (e.g. floats)
	'''
	if isinstance(value, bool):
		raise TypeError(f"'{name}' must be an integer, not a bool.")
	try:
		return operator.index(value) # accepts NumPy ints; rejects floats — no silent truncation
	except TypeError:
		raise TypeError(f"'{name}' must be an integer (got '{type(value).__name__}').") from None
