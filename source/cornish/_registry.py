'''
Class-name -> wrapper-class dispatch for AST objects (SPEC-03A §1).

Registration is an import side effect: each wrapper module decorates its class
with :func:`register_region` / :func:`register_frame`. The cornish package
``__init__`` imports every wrapper module, so any ``import cornish`` guarantees
a fully-populated registry. A test (the T6 gate) asserts the exact expected
table contents, so a module that silently stops registering fails the suite.

This module imports NOTHING from the wrapper modules (they import it), so no
circular imports are possible.

Dispatch is on ``ast_object.Class`` — the authoritative AST class name — and is
therefore immune to the isinstance-ordering pitfalls of the hand-written
ladders it replaces (e.g. ``Ast.DSBSpecFrame`` is an ``Ast.SpecFrame``
subclass, so an isinstance ladder testing SpecFrame first misroutes it).
'''

from __future__ import annotations

import logging
from typing import Callable, TYPE_CHECKING

import starlink.Ast as Ast

from .exc import UnsupportedASTClass

if TYPE_CHECKING:  # pragma: no cover - typing only
	from .region.region import ASTRegion
	from .mapping.frame.frame import ASTFrame

__all__ = ['register_region', 'register_frame', 'region_from_ast', 'frame_from_ast']

logger = logging.getLogger("cornish")

#: ``.Class`` name -> cornish region wrapper, e.g. ``'Polygon' -> ASTPolygon``.
_AST_REGION_CLASSES: dict[str, type] = {}

#: ``.Class`` name -> cornish frame wrapper, e.g. ``'SkyFrame' -> ASTSkyFrame``.
_AST_FRAME_CLASSES: dict[str, type] = {}


def _register(table: dict, table_name: str, ast_class_name: str) -> Callable[[type], type]:
	'''
	Shared body of the two registration decorators.

	:raises ValueError: a DIFFERENT class is already registered under this name
	'''
	if not isinstance(ast_class_name, str) or not ast_class_name:
		raise TypeError(f"An AST class name (a non-empty string) is required (got {ast_class_name!r}).")

	def decorator(wrapper_class: type) -> type:
		if not isinstance(wrapper_class, type):
			raise TypeError(f"Only classes can be registered (got '{type(wrapper_class).__name__}').")
		existing = table.get(ast_class_name)
		if existing is not None and existing is not wrapper_class:
			# double registration of the SAME class is a no-op (module re-imports
			# are harmless); a different class is a genuine conflict
			raise ValueError(
				f"The AST {table_name} class '{ast_class_name}' is already registered to "
				f"'{existing.__name__}'; cannot also register '{wrapper_class.__name__}'.")
		table[ast_class_name] = wrapper_class
		return wrapper_class

	return decorator


def register_region(ast_class_name: str) -> Callable[[type], type]:
	'''
	Class decorator: register a wrapper for the AST region class whose
	``.Class`` attribute equals ``ast_class_name``.

	The decorated class MUST accept ``cls(ast_object=<pyast object of that
	class>)`` — that is the registry's construction contract.

	:raises ValueError: a DIFFERENT class is already registered under this name
		(double registration of the same class is a no-op, so module re-imports
		are harmless)
	'''
	return _register(_AST_REGION_CLASSES, "region", ast_class_name)


def register_frame(ast_class_name: str) -> Callable[[type], type]:
	''' Same contract as :func:`register_region`, for frame wrappers. '''
	return _register(_AST_FRAME_CLASSES, "frame", ast_class_name)


def _reject_wrapper(ast_object) -> None:
	'''
	Raise if ``ast_object`` is a cornish wrapper rather than a raw pyast object.

	Imported lazily: ``ast_object.py`` is imported by everything, and this
	module must stay import-free of the wrapper layer at module level.
	'''
	from .ast_object import ASTObject
	if isinstance(ast_object, ASTObject):
		raise TypeError(
			f"'{type(ast_object).__name__}' is already a cornish wrapper; "
			f"pass its .astObject if you need to re-wrap it.")


def region_from_ast(ast_object: Ast.Region) -> "ASTRegion":
	'''
	The one place a pyast Region becomes a cornish region.

	Dispatch is on ``ast_object.Class`` (authoritative; immune to the
	isinstance-order pitfalls of the hand-written ladders). The returned
	wrapper holds ``ast_object`` BY REFERENCE — callers that need isolation
	from live pointers (e.g. :attr:`ASTCompoundRegion.components`) copy first.

	:raises TypeError: ``ast_object`` is not a ``starlink.Ast.Region`` (a
		cornish wrapper is also a TypeError — "already wrapped"), with the
		received type named
	:raises cornish.exc.UnsupportedASTClass: the class has no registered
		wrapper; never a silent generic wrapper
	'''
	_reject_wrapper(ast_object)
	if not isinstance(ast_object, Ast.Region):
		raise TypeError(
			f"Expected a starlink.Ast.Region to wrap (got '{type(ast_object).__name__}').")

	ast_class = ast_object.Class
	wrapper_class = _AST_REGION_CLASSES.get(ast_class)
	if wrapper_class is None:
		raise UnsupportedASTClass(
			f"No cornish region wrapper is registered for the AST class '{ast_class}'. "
			f"Wrappers for the remaining region classes (Ellipse, Interval, NullRegion, "
			f"PointList, Prism, ...) are specified in SPEC-05.")
	return wrapper_class(ast_object=ast_object)


def frame_from_ast(ast_object: Ast.Frame) -> "ASTFrame":
	'''
	The one place a pyast Frame becomes a cornish frame.

	Unlike regions, unregistered PLAIN frame classes fall back to the generic
	:class:`~cornish.mapping.frame.frame.ASTFrame` wrapper — a generic ASTFrame
	is fully functional for any frame (attribute access, unit handling),
	whereas no generic concrete ASTRegion exists. The fallback is logged at
	DEBUG level with the class name (visible, not silent).

	:raises TypeError: not an ``Ast.Frame``; or an ``Ast.Region`` (a Region
		passes ``isinstance(x, Ast.Frame)`` by MRO but wrapping one as a frame
		is a category error — call ``.getregionframe()`` first); or a cornish
		wrapper
	'''
	_reject_wrapper(ast_object)
	if isinstance(ast_object, Ast.Region):
		raise TypeError(
			f"A region ('{ast_object.Class}') is not a frame, even though it satisfies "
			f"isinstance(x, Ast.Frame); call .getregionframe() to obtain the frame it is "
			f"defined in and wrap that instead.")
	if not isinstance(ast_object, Ast.Frame):
		raise TypeError(
			f"Expected a starlink.Ast.Frame to wrap (got '{type(ast_object).__name__}').")

	ast_class = ast_object.Class
	wrapper_class = _AST_FRAME_CLASSES.get(ast_class)
	if wrapper_class is None:
		# the generic ASTFrame answers every question a caller can ask of a
		# frame; the fallback is recorded, never silent
		logger.debug("No cornish frame wrapper is registered for the AST class '%s'; "
		             "returning a generic ASTFrame.", ast_class)
		wrapper_class = _AST_FRAME_CLASSES.get("Frame")
		if wrapper_class is None:
			# only reachable if frame.py stopped registering ASTFrame; loud, never
			# a bare KeyError from the fallback lookup
			raise UnsupportedASTClass(
				"The generic frame wrapper ('Frame') is not registered — the cornish frame "
				"registry is not populated (was cornish imported?).")
	return wrapper_class(ast_object=ast_object)
