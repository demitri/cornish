#/usr/bin/env python

from typing import Union, Iterable, Container

import starlink.Ast as Ast

from ... import ASTObject
from ..mapping import ASTMapping
from ... import _pyast_bridge as bridge
from ..._validation import as_integer

import numpy as np
import astropy
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

__all__ = ['ASTFrame']

class ASTFrame(ASTMapping):
	'''
	A Frame is a representation of a coordinate system, e.g. Cartesian, RA/dec.
	It contains information about the labels which appear on the axes, the axis units,
	a title, knowledge of how to format the coordinate values on each axis, etc.

	List and description of ``starlink.Ast.Frame`` attributes in documentation: Section 7.5.

	Ref:
	http://www.starlink.rl.ac.uk/docs/sun95.htx/sun95se27.html
	http://www.strw.leidenuniv.nl/docs/starlink/sun210.htx/node71.html

	:param ast_object: an existing :class:`starlink.Ast.Frame` object
	'''
	def __init__(self, ast_object:Ast.Frame=None, naxes:int=None):

		if ast_object is not None and naxes is not None:
			raise ValueError("Cannot initialize ASTFrame with both 'ast_object' and 'naxes' set.")

		if all([x is None for x in [ast_object, naxes]]):
			raise ValueError("Either 'ast_object' or 'naxes' must be specified to create an ASTFrame.")

		if ast_object is not None:
			if isinstance(ast_object, Ast.Frame):
				super().__init__(ast_object=ast_object)
			else:
				raise TypeError(f"The provided 'ast_object' is not an Ast.Frame object (got class '{type(ast_object)}').")
		else:
			super().__init__(ast_object=Ast.Frame(naxes))

	@staticmethod
	def frameFromAstObject(ast_object:Ast.Frame=None):
		'''
		Factory method that returns the appropriate Cornish frame object (e.g. :class:`ASTSkyFrame`) for a given frame.

		:param ast_object: an :py:class:`Ast.Frame` object
		'''
		if ast_object is None:
			raise ValueError("An ast_object must be specified.")
		elif not isinstance(ast_object, Ast.Frame):
			raise TypeError(f"Expected 'ast_object' to be some kind of Ast frame, but got {type(ast_object)}")

		# the order might need to be rearranged depending on the object hierarchy

		if isinstance(ast_object, Ast.SkyFrame):
			from .sky_frame import ASTSkyFrame
			return ASTSkyFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.SpecFrame):
			from .spec_frame import ASTSpecFrame
			return ASTSpecFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.TimeFrame):
			from .time_frame import ASTTimeFrame
			return ASTTimeFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.CmpFrame):
			from .compound_frame import ASTCompoundFrame
			return ASTCompoundFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.FluxFrame):
			raise NotImplementedError("ASTFluxFrame not yet implemented.")
			#return ASTFluxFrame(ast_object=ast_object)
		elif isinstance(ast_object, Ast.SpecFluxFrame):
			raise NotImplementedError("ASTSpecFluxFrame not yet implemented.")
			#return ASTSpecFluxFrame(ast_object=ast_object)
		else:
			return ASTFrame(ast_object=ast_object)

	@staticmethod
	def frameFromFITSHeader(header):
		'''
		Factory method that returns a new ASTFrame from the provided FITS header.
		'''
		from ...channel import ASTFITSChannel
		fitschannel = ASTFITSChannel(header=header)
		return fitschannel.frameSet

	@property
	def naxes(self) -> int:
		''' Returns the number of axes for the frame. '''
		return int(self.astObject.get("Naxes"))

	@property
	def title(self) -> str:
		''' Returns the frame title, a string describing the coordinate system which the frame represents. '''
		return self.astObject.get("Title")

	@title.setter
	def title(self, newTitle):
		self._setAttribute("Title", newTitle)

	def _validate_axis(self, axis) -> int:
		''' Shared validation for the per-axis accessors below; returns the normalized axis. '''
		if axis is None:
			raise ValueError("An axis number must be specified.")
		axis = as_integer(axis, "axis")
		if not (1 <= axis <= self.naxes):
			# axes are 1-based; out-of-range values would leak a raw Ast.AXIIN
			raise ValueError("The axis provided ({0}) must be between 1 and the number of axes ({1}).".format(axis, self.naxes))
		return axis

	def label(self, axis=None) -> str:
		''' Return the label for the specified axis. '''
		axis = self._validate_axis(axis)
		return self.astObject.get("Label({0})".format(axis))

	def setLabelForAxis(self, axis=None, label=None):
		''' Set the label for the specified axis. '''
		axis = self._validate_axis(axis)
		if label is None:
			raise ValueError("A new label must be specified.")
		self._setAttribute(f"Label({axis})", label)

	def unit(self, axis=None):
		''' Return the unit for the specified axis. '''
		axis = self._validate_axis(axis)
		return self.astObject.get("Unit({0})".format(axis))

	def setUnitForAxis(self, axis=None, unit=None):
		''' Set the unit as a string value for the specified axis. '''
		axis = self._validate_axis(axis)
		if unit is None:
			raise ValueError("A new unit must be specified.")
		self._setAttribute(f"Unit({axis})", unit)

	@property
	def system(self):
		'''
		String which identifies the coordinate system represented by the Frame.
		The system is ``Cartesian`` by default, but can have other values for subclasses of Frame,
		e.g. ``FK4``, ``Galactic``.
		'''
		return self.astObject.get("System")

	@system.setter
	def system(self, system=None):
		''' Set the a label for the frame's system. '''
		if system is None:
			raise ValueError("A 'system' parameter must be specified.")
		elif not isinstance(system, str):
			raise TypeError("A string value was expected for 'system' (got '{0}').".format(type(system)))
		self._setAttribute("System", system)

	@property
	def isSkyFrame(self) -> bool:
		'''
		Returns ``True`` if this is a SkyFrame, ``False`` otherwise.
		'''
		return self.astObject.isaskyframe()

	@property
	def domain(self) -> str:
		'''
		The physical domain of the coordinate system (string value).
		The Domain attribute also controls how Frames align with each other.
		If the Domain value in a Frame is set, then only Frames with the same Domain value can be aligned with it.

		Example values: ``GRID``, ``FRACTION``, ``PIXEL``, ``AXIS``, ``SKY``, ``SPECTRUM``, ``CURPIC``, ``NDC``, ``BASEPIC``, ``GRAPHICS``

		Frames created by the user (for instance, using WCSADD) can have any Domain value, but the standard
		domain names should be avoided unless the standard meanings are appropriate for the Frame being created.
		'''
		return self.astObject.get("Domain")

	@domain.setter
	def domain(self, newDomain=None):
		if newDomain is None or isinstance(newDomain, str) == False:
			raise TypeError("The domain value must be set to a string.")
		self._setAttribute("Domain", newDomain)

	def distance(self, point1:Union[Iterable, SkyCoord], point2:Union[Iterable,SkyCoord]) -> Quantity:
		'''
		Distance between two points in this frame.

		On a sky frame, points may be given as SkyCoords (any system — converted),
		Quantities, or bare values read as degrees; the returned distance is a
		Quantity in degrees. On a non-sky frame, points must be bare values in
		the frame's native units (a Quantity raises ValueError — its unit cannot
		be reconciled with unknowable native units), and the returned distance
		is native (a pixel Quantity where the frame is identifiably a pixel grid).

		:param point1: a two element list/tuple/Numpy array, SkyCoord, or Quantity pair of the first point coordinates
		:param point2: a two element list/tuple/Numpy array, SkyCoord, or Quantity pair of the second point coordinates
		'''
		p1 = bridge.to_frame_units(point1, self.astObject, squeeze=True)
		p2 = bridge.to_frame_units(point2, self.astObject, squeeze=True)

		distance = self.astObject.distance(p1, p2)

		if bridge.is_sky(self.astObject):
			return bridge.from_frame_distance(distance, self.astObject) * u.deg

		# try to add a unit if it can be deduced
		if self.system == 'Cartesian' and self.domain == "GRID":
			distance = distance * u.pixel
		return distance

	def angle(self, vertex:Iterable=None, points:Container[Union[SkyCoord,Iterable]]=None) -> Quantity:
		'''
		Calculate the angle in this frame between two line segments connected by a point.

		Let A = point1, C = point2, and B = the vertex point. This method calculates the
		angle between the line segments :math:`\\overline{AB}` and :math:`\\overline{CB}`.

		If the frame is a sky frame, lines are drawn on great circles.
		On sky frames, units are assumed to be degrees if not provided with units,
		e.g. as an astropy.coordinates.SkyCoord or astropy.units.Quantity values;
		on non-sky frames, bare values are read in the frame's native units.

		:param vertex: a two element list/tuple/Numpy array or a SkyCoord of the vertex
		:param points: a two element list/tuple/etc. containing two points in this frame
		'''
		point1, point2 = points

		p1 = bridge.to_frame_units(point1, self.astObject, squeeze=True)
		p2 = bridge.to_frame_units(point2, self.astObject, squeeze=True)
		v = bridge.to_frame_units(vertex, self.astObject, squeeze=True)

		angle_rad = self.astObject.angle(p1, v, p2)

		# The RETURNED ANGLE is genuinely radians in both frame kinds — this is
		# geometry of the result, not a coordinate conversion, so it stays at
		# this call site (allowlisted in the conversions-outside-bridge gate).
		return np.rad2deg(angle_rad) * u.deg

	def offsetAlongGeodesicCurve(self, point1:Iterable, point2:Iterable, offset:Union[float, u.Quantity]):
		'''
		Return the point offset along the geodesic curve from ``point1`` toward ``point2``.

		In a sky frame, the line will be curved (a great circle) and the offset
		may be a Quantity (angular) or a bare number read as degrees. In a basic
		frame, the line will be straight and the offset must be a bare number in
		the frame's native units (a Quantity raises ValueError).

		:param point1: a two element list/tuple/NumPy array of the first point coordinates
		:param point2: a two element list/tuple/NumPy array of the second point coordinates
		:param offset: a distance along the geodesic connecting the two points
		'''
		p1 = bridge.to_frame_units(point1, self.astObject, squeeze=True)
		p2 = bridge.to_frame_units(point2, self.astObject, squeeze=True)
		off = bridge.to_frame_distance(offset, self.astObject)

		out = self.astObject.offset(p1, p2, off)
		result = bridge.from_frame_units(out, self.astObject)

		if bridge.is_sky(self.astObject):
			return result * u.deg
		return result

	def framesetWithMappingTo(self, template_frame:"ASTFrame"=None) -> Union["ASTFrame", None]: # maybe think of better name?!
		'''
		Search this frame (or set) to identify a frame that shares the same coordinate system as the provided template frame.

		For example, this method can be used to see if this frame
		(or frame set) contains a sky frame.

		Returns ``None`` if no mapping can be found.

		:param template_frame: an instance of the type of frame
		  being searched for
		:returns: a frame that matches the template
		'''
		from .frame_set import ASTFrameSet
		if template_frame is None:
			raise ValueError("A template frame must be provided.")
		return ASTFrameSet(ast_object=self.astObject.findframe(template_frame.astObject))
