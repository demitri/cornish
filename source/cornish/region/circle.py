
from __future__ import annotations # remove in Python 3.10
# Needed for forward references, see:
# https://stackoverflow.com/a/33533514/2712652

import logging
import numbers
from typing import Union, Iterator, Tuple

import astropy
import astropy.units as u

import numpy as np
import starlink
import starlink.Ast as Ast

from .region import ASTRegion, DEFAULT_UNCERTAINTY
from .polygon import ASTPolygon
from ..mapping.frame import ASTFrame
from ..mapping.frame.sky_frame import ASTICRSFrame
from .. import _pyast_bridge as bridge

__all__ = ["ASTCircle"]

CENTER_EDGE = 0
CENTER_RADIUS = 1

logger = logging.getLogger("cornish")

class ASTCircle(ASTRegion):
	'''
	ASTCircle is an :class:`ASTRegion` that represents a circle.

	Accepted signatures for creating an ASTCircle

	.. code-block::python

		c = ASTCircle(ast_object)     # where ast_object is a starlink.Ast.Circle object
		c = ASTCircle(frame, center, edge_point)
		c = ASTCircle(frame, center, radius)

	Points and the radius may be given in any of the bridge-accepted forms:
	SkyCoord (any system — converted), Quantities, or bare values read as
	degrees on sky frames and native units otherwise.

	:param ast_object: a circle object from the ``starlink-pyast`` module
	:param frame: a frame the circle is to be defined in — a cornish or ``starlink.Ast`` frame, frame set (current frame governs), or region; uses :class:`ASTICRSFrame` if `None`
	:param center: two elements that describe the center point of the circle in the provided frame
	:param edge_point: two elements that describe a point on the circumference of the circle in the provided frame
	:param radius: radius of the circle: a Quantity (angular, sky frames only) or a number (degrees on sky frames, native units otherwise)
	:param uncertainty: uncertainty of the boundary: a Quantity (angular, sky frames only), a number (degrees on sky frames), or a Region; omit for the 1 arcsec default (sky frames), or pass ``None`` explicitly for AST's internal default
	'''
	def __init__(self, \
				 ast_object:starlink.Ast.Circle=None, \
				 frame:Union[ASTFrame,Ast.Frame,astropy.coordinates.ICRS]=None, \
				 center:Union[astropy.coordinates.SkyCoord,Iterator]=None, \
				 edge_point:Union[np.ndarray,Iterator]=None, \
				 radius:Union[float,astropy.units.quantity.Quantity]=None, \
				 uncertainty:Union[None, float, u.Quantity, ASTRegion, Ast.Region]=DEFAULT_UNCERTAINTY):

		if ast_object is not None:
			if any([x is not None for x in [frame, center, edge_point, radius]]):
				raise ValueError("Cannot specify both 'ast_object' and any other parameter.")

			if isinstance(ast_object, Ast.Circle):
				# an externally-created circle: its uncertainty (if any) is unknown here
				super().__init__(ast_object=ast_object)
				return
			else:
				raise TypeError("The 'ast_object' provided was not of type starlink.Ast.Circle.")

		# check valid combination of parameters
		# -------------------------------------

		# make sure we have a frame we can work with
		if frame is None:
			ast_frame = ASTICRSFrame().astObject
		elif isinstance(frame, astropy.coordinates.ICRS):
			ast_frame = ASTICRSFrame().astObject
		else:
			ast_frame = bridge._unwrap(frame) # frame sets and regions are legal: the current/encapsulated frame governs

		if all([x is not None for x in [edge_point, radius]]):
			raise ValueError("Both 'edge_point' and 'radius' cannot be simultaneously specified.")
		if center is None:
			raise ValueError("The 'center' parameter must be set.")
		if all([x is None for x in [edge_point, radius]]):
			raise ValueError("Along with 'center', a 'radius' or 'edge_point' must be specified.")

		# input forms:
		#	CENTER_EDGE   (0) : circle specified by center point and any point on the circumference (p1 = [float,float], p2 = [float,float])
		#	CENTER_RADIUS (1) : circle specified by center point and radius                         (p1 = [float,float], p2 = float)
		p1 = bridge.to_frame_units(center, ast_frame, squeeze=True)
		if edge_point is not None:
			input_form = CENTER_EDGE
			p2 = bridge.to_frame_units(edge_point, ast_frame, squeeze=True)
		else:
			input_form = CENTER_RADIUS
			p2 = [bridge.to_frame_distance(radius, ast_frame)]

		if uncertainty is DEFAULT_UNCERTAINTY and not bridge.is_sky(ast_frame):
			uncertainty = None # the angular 1 arcsec default is meaningless on a basic frame; AST's internal default governs
		unc_region = bridge.as_uncertainty_region(uncertainty, ast_frame, p1)

		# pyast SILENTLY IGNORES the `unc=` keyword on region constructors and
		# rejects a positional None, so the uncertainty region must be appended
		# POSITIONALLY and only when non-None (SPEC-04A §8 / D15). Do not "clean
		# this up" into the kwarg form — it parses fine and silently drops the
		# uncertainty.
		args = (ast_frame, input_form, p1, p2) + ((unc_region,) if unc_region is not None else ())
		super().__init__(ast_object=Ast.Circle(*args), uncertainty=uncertainty)

	def __repr__(self):
		return "<{0}.{1} {2}: center={3}, r={4:0.6}>".format(self.__class__.__module__, self.__class__.__name__, hex(id(self)),
														self.center, self.radius)

	def _validateDilationOperand(self, value, operation:str):
		'''
		Shared validation for the dilation operators below; returns the numeric operand.
		'''
		if isinstance(value, (ASTRegion, Ast.Region)):
			raise TypeError(f"A region cannot be used to {operation} a circle; "
			                f"pass a number (or an angular Quantity on a sky frame).")
		if isinstance(value, bool) or not isinstance(value, (numbers.Real, u.Quantity)):
			raise TypeError(f"Cannot {operation} a circle's radius by an object of type '{type(value).__name__}'.")
		return value

	def __add__(self, value):
		'''
		Return a new ASTCircle with its radius increased by the provided distance.

		On a sky-frame circle the value may be an angular :class:`astropy.units.Quantity`
		or a bare number read as degrees. On a basic-frame circle the value must be
		a bare number in the frame's native units (a Quantity raises ValueError —
		its unit cannot be reconciled with unknowable native units).
		'''
		value = self._validateDilationOperand(value, "dilate")
		if bridge.is_sky(self.astObject):
			if not isinstance(value, u.Quantity):
				value = float(value) * u.deg
			return ASTCircle(frame=self.frame(), center=self.center, radius=self.radius + value)
		else:
			if isinstance(value, u.Quantity):
				raise ValueError("A Quantity cannot dilate a circle on a non-sky frame: its unit "
				                 "cannot be reconciled with unknowable native frame units — pass a bare number.")
			return ASTCircle(frame=self.frame(), center=self.center, radius=self.radius + float(value))

	def __mul__(self, value):
		'''
		Return a new ASTCircle with a radius scaled by the value provided.
		:param value: numeric value to scale the radius by, e.g. 1.2 increases the radius by 20%
		'''
		value = self._validateDilationOperand(value, "scale")
		if isinstance(value, u.Quantity):
			raise TypeError("A circle's radius can only be scaled by a plain number, not a Quantity.")
		if not (np.isfinite(value) and value > 0):
			raise ValueError(f"A circle's radius can only be scaled by a positive, finite value (got {value!r}).")
		return ASTCircle(frame=self.frame(), center=self.center, radius=self.radius * float(value))

	def __truediv__(self, value):
		'''
		Return a new ASTCircle with a radius divided by the value provided.
		:param value: numeric value to divide the radius by, e.g. 2 decreases the radius by 50%
		'''
		value = self._validateDilationOperand(value, "scale")
		if isinstance(value, u.Quantity):
			raise TypeError("A circle's radius can only be scaled by a plain number, not a Quantity.")
		if not (np.isfinite(value) and value > 0):
			raise ValueError(f"A circle's radius can only be divided by a positive, finite value (got {value!r}).")
		return ASTCircle(frame=self.frame(), center=self.center, radius=self.radius / float(value))

	@property
	def radius(self) -> Union[astropy.units.quantity.Quantity, float]:
		'''
		The radius of this circle region.

		:returns: for a circle on a sky frame, an :class:`astropy.units.Quantity`
			in degrees (a geodesic distance); for a circle on a basic (non-sky)
			frame, a bare float in the frame's native units — never labeled with
			a unit that AST does not know
		'''
		( center, radius, some_point_on_circumference ) = self.astObject.circlepars()

		r = bridge.from_frame_distance(radius, self.astObject)
		if bridge.is_sky(self.astObject):
			return r * u.deg
		return r

	@property
	def center(self) -> Tuple[float]:
		'''
		The center of this circle region (a synonym for :func:`self.centre`" for the Americans).

		Returns
		-------
		:returns: a tuple of points (x,y) that describe the centre of the circle, in degrees if a sky frame
		'''
		return self.centre

	@property
	def centre(self) -> Tuple[float]:
		'''
		The centre of this circle region, in degrees if a sky frame (native frame units otherwise).

		Returns
		-------
		:returns: a tuple of points (x,y) that describe the centre of the circle
		'''
		( center, radius, some_point_on_circumference ) = self.astObject.circlepars()

		return tuple(bridge.from_frame_units(center, self.astObject))

	def toPolygon(self, npoints=200, maxerr:astropy.units.Quantity=1.0*u.arcsec):
		'''
		Returns a new polygon region that approximates this circle in the same frame.

		The algorithm used in this method leads to the new polygon being fully inscribed by the
		originating circle; all points generated are on the circle's circumference. Although the
		default number of points is 200, typically a much smaller number (e.g. 20) is
		sufficient.

		:param npoints: number of points to sample from the circle for the resulting polygon
		:param maxerr:
		'''
		points = self.boundaryPointMesh(npoints=npoints)

		polygon = ASTPolygon(frame=self, points=points) # can get the frame from the astObject
		polygon.wcs = self.wcs # propagate the originating WCS, when known (SPEC-04 §2)
		return polygon

	def boundingCircle(self) -> ASTCircle:
		'''
		This method returns itself; a circle region is its own bounding circle.
		'''
		return self

	@property
	def area(self):
		'''
		The area of the circle within its frame (e.g. on a Cartesian plane or sphere). [Not yet implemented.]
		'''
		# see: https://math.stackexchange.com/questions/1832110/area-of-a-circle-on-sphere
		raise NotImplementedError()
