#/usr/bin/env python

import logging
import warnings
from typing import Union, Iterable, Optional

import numpy as np
import starlink
import starlink.Ast as Ast
import astropy
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from ... import ASTObject
from .frame import ASTFrame
from ..mapping import ASTMapping
from ...exc import FrameNotFoundException
from ... import _pyast_bridge as bridge

__all__ = ['ASTFrameSet']

logger = logging.getLogger("cornish") # cornish logger

class ASTFrameSet(ASTFrame):
	'''
	Create a new AST frame set.
	Object can be created from an :class:`starlink.Ast.FrameSet` "primitive"
	(e.g. returned by another object).

	A set of inter-related coordinate systems made up of existing mapping's and frames.
	A FrameSet may be extended by adding a new Frame and associated Mapping.

	A FrameSet must have a "base" frame which represents the "native" coordinate system
	(for example, the pixel coordinates of an image). Similarly, one Frame is termed the
	current Frame and represents the "currently-selected" coordinates. It might typically
	be a celestial or spectral coordinate system and would be used during interactions
	with a user, as when plotting axes on a graph or producing a table of results.
	Other Frames within the FrameSet represent a library of alternative coordinate systems
	which a software user can select by making them current.

	Accepted signatures for creating an ``ASTFrameSet``:

	.. code-block::python

	    fs = ASTFrameSet(ast_object)
	    fs = ASTFrameSet(base_frame)

	:param ast_object: an :class:`Ast.astFrame` object from the `starlink-pyast` library
	:param base_frame: base frame to create the FrameSet from
	'''
	def __init__(self, ast_object:starlink.Ast.FrameSet=None, base_frame:Union[starlink.Ast.FrameSet,ASTFrame]=None): #, fits_header=None):
		# validate parameters
		if all([ast_object, base_frame]):
			raise ValueError("Both an 'ast_object' and a 'base_frame' cannot be specified.")
		elif all([x is None for x in [ast_object, base_frame]]):
			raise ValueError("One of 'ast_object' or 'base_frame' must be provided.")

		if ast_object:
			if isinstance(ast_object, starlink.Ast.FrameSet):
				super().__init__(ast_object=ast_object)
			else:
				raise TypeError("Unhandled ast_object type ('{0}')".format(ast_object))
		else:
			# construct from provided base_frame
			# NOTE: the frame MUST be passed positionally — pyast silently swallows
			# the `frame=` keyword (same bug class as `unc=`/`oper=`), leaving the
			# constructor with zero arguments (verified 2026-07-10: the kwarg form
			# raises "takes at least 1 argument (0 given)").
			if isinstance(base_frame, starlink.Ast.Frame):
				fs = Ast.FrameSet(base_frame)
			elif isinstance(base_frame, ASTFrame):
				fs = Ast.FrameSet(base_frame.astObject)
			else:
				raise TypeError("Unhandled base_frame type ('{0}')".format(base_frame))

			super().__init__(ast_object=fs)

	@staticmethod
	def fromFrames(frame1:Union[ASTFrame, starlink.Ast.Frame], frame2:Union[ASTFrame, starlink.Ast.Frame]):
		'''
		Static method to create a frame set from two existing frames.

		A frame set is a mapping between two frames that can convert coordinates from
		one frame (the "base" frame) to the other frame (the "current" frame).

		:param frame1: the "base" frame (frame to convert coordinates from)
		:param frame2: the "current" frame (frame to convert coordinates to)
		'''
		# get Ast objects for each frame
		ast_frames = list()

		for frame in [frame1,frame2]:
			if isinstance(frame, ASTFrame):
				ast_frames.append(frame.astObject)
			elif isinstance(frame, starlink.Ast.Frame):
				ast_frames.append(frame)
			else:
				raise TypeError(f"The provided frames must either be of type ASTFrame or starlink.Ast.Frame (got '{type(frame)}'.")

		frame_set = Ast.Frame.convert(ast_frames[0], ast_frames[1])
		if frame_set is None:
			from ...exc import CoordinateSystemsCouldNotBeMapped
			raise CoordinateSystemsCouldNotBeMapped("An ASTFrameSet could not be created since the mapping between the two provided frames could not be determined (conversion may not be possible).")
		else:
			return ASTFrameSet(ast_object=frame_set)

	@staticmethod
	def fromFITSHeader(fits_header=None):
		'''
		Static method that returns a FrameSet object read from the provided FITS header.
		'''

		if fits_header is None:
			raise ValueError("A FITS header must be provided.")

		from ...channel import ASTFITSChannel # or "from cornish import ..." ?

		fits_channel = ASTFITSChannel(header=fits_header)

		# does this channel contain a frame set?
		frame_set = fits_channel.frameSet
		if frame_set:
			return frame_set
		else:
			raise FrameNotFoundException("A valid frame set could not be read from the provided FITS header (no WCS?).")

	def _get_frame(self, frame_index:int) -> starlink.Ast.Frame:
		'''
		Internal method to retrieve a ``starlink.Ast.Frame`` frame by its index value within this frame set.
		'''
		return self.astObject.getframe(frame_index)

	def frameAtIndex(self, frame_index:int) -> ASTFrame:
		'''
		Return the frame at the specified index within this frame set.
		'''
		frame = self.astObject.getframe(frame_index)

		if frame.isaskyframe():
			from .sky_frame import ASTSkyFrame
			return ASTSkyFrame(ast_object=frame)
		elif frame.isacmpframe():
			#from .compound_frame import ASTCompoundFrame
			#return ASTCompoundFrame(ast_object=frame)
			logger.info("A compound frame is being returned as an ASTFrame instead of the ASTCompoundFrame subclass (not yet implemented).")
			pass
		elif frame.isadsbspecframe(): # dual sideband instrument
			logger.info("A compound frame is being returned as an ASTFrame instead of the ASTDSBSpectrumFrame subclass (not yet implemented).")
			pass
		elif frame.isafluxframe():
			logger.info("A compound frame is being returned as an ASTFrame instead of the ASTFluxFrame subclass (not yet implemented).")
			pass
		elif frame.isaspecfluxframe():
			logger.info("A compound frame is being returned as an ASTFrame instead of the ASTSpectrumFluxFrame subclass (not yet implemented).")
			pass
		elif frame.isaspecframe():
			logger.info("A compound frame is being returned as an ASTFrame instead of the ASTSpectrumFrame subclass (not yet implemented).")
			pass
		elif frame.isatimeframe():
			#from .time_frame import ASTTimeFrame
			logger.info("A compound frame is being returned as an ASTFrame instead of the ASTTimeFrame subclass (not yet implemented).")
			pass

		return ASTFrame(ast_object=frame)

	@property
	def baseFrame(self):
		'''
		Return the base frame.
		'''
		return self.frameAtIndex(Ast.BASE)

	@property
	def currentFrame(self):
		''' Returns the current frame. '''
		return self.frameAtIndex(Ast.CURRENT)

	def removeCurrentFrame(self):
		''' Remove the current frame from the frame set. '''
		# todo: check if the frame is actually part of the frame set
		# Note: The "removeframe" documenation is a little unclear on how to remove
		#       a frame that is not current.
		# Better to make this "removeFrame(frame=...) or by frame name, or add a new method.
		self.astObject.removeframe(iframe=Ast.CURRENT)

	def centerCoordinates(self):
		''' Returns the coordinates at the center of the frame. '''
		raise NotImplementedError("Useful, not sure how to do this.")

	def addToBaseFrame(self, frame=None):
		'''
		Add a new frame to this frame set's base frame.

		'''

		# Ref: http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211ss5.html#xref_astAddFrame
		# addframe( iframe, map, frame )

		if frame is None:
			raise ValueError("frame must be provided to 'addToBaseFrame'")

		iframe = Ast.BASE # add to base frame
		map = None
		self.astObject.addframe(iframe, map, frame)

	@property
	def mapping(self):
		'''
		Return an object that maps points from the base frame to the current frame of this frame set.
		'''
		return ASTMapping(ast_object=self.astObject.getmapping()) # default is from base to current

	def _tran(self, points, forward:bool=True, *,
	          parallel_axes:Optional[bool]=None, bad:str="raise") -> np.ndarray:
		'''
		The single private implementation behind pix2world/world2pix/convertPoints:
		bridge in -> Ast.tran -> bridge out, preserving single-point-in ->
		single-point-out return shapes by FORM (SPEC-04A §5).
		'''
		src = self.astObject.getframe(Ast.BASE if forward else Ast.CURRENT)
		dst = self.astObject.getframe(Ast.CURRENT if forward else Ast.BASE)
		single = bridge.is_single_point(points, src)          # decided by FORM, before conversion
		pts = bridge.to_frame_units(points, src, parallel_axes=parallel_axes)   # always (naxes, npoints) — tran requires 2-D (V5)
		out = bridge.from_frame_units(self.astObject.tran(pts, forward), dst, bad=bad)
		return out[0] if single else out               # single-point form in -> (naxes,) out

	def convertPoints(self, points:Iterable, forward:bool=True, *, parallel_axes:Optional[bool]=None):
		'''
		Convert the provided coordinate points from this frame set's base frame to
		its current frame (``forward=True``) or the reverse (``forward=False``).

		Accepted input forms are those of the bridge decision table: SkyCoord
		(scalar, array, or a sequence of SkyCoords — sky source frames only),
		Quantity arrays (sky source frames only), bare arrays/lists read in the
		source frame's units (degrees when the source frame is a sky frame),
		in either pairs ``(n, naxes)`` or parallel ``(naxes, n)`` orientation.
		A square ``(naxes, naxes)`` array reads as PAIRS (DECISIONS D3); pass
		``parallel_axes=True`` for the parallel reading.

		Returns a bare array of converted points — shape ``(npoints, naxes)``, or
		``(naxes,)`` for single-point-form input — in the destination frame's
		units (degrees for sky). As a convenience, when the input was SkyCoord-form
		and the destination frame is a sky frame, SkyCoords are returned instead,
		built in the destination frame's actual system (never ICRS-by-fiat):
		scalar in -> scalar out, array in -> array out, sequence in -> list of
		scalar SkyCoords out.

		:param points: coordinate points in the source frame
		:param forward: True: base -> current; False: current -> base
		:param parallel_axes: orientation override for square array input (see above)
		'''
		skycoord_form = isinstance(points, SkyCoord) or \
			(isinstance(points, (list, tuple)) and len(points) > 0 and all(isinstance(p, SkyCoord) for p in points))

		out = self._tran(points, forward, parallel_axes=parallel_axes)

		if not skycoord_form:
			return out

		dst = self.astObject.getframe(Ast.CURRENT if forward else Ast.BASE)
		if not bridge.is_sky(dst):
			return out  # nothing to wrap: SkyCoords cannot represent non-sky coordinates

		# build output SkyCoords in the DESTINATION's actual system (§2.3, both
		# directions); an unmappable sky destination raises the §2.3 ValueError
		target = bridge.astropy_frame_for(dst)

		if isinstance(points, SkyCoord):
			if points.isscalar:
				return SkyCoord(out[0] * u.deg, out[1] * u.deg, frame=target)
			return SkyCoord(out[:, 0] * u.deg, out[:, 1] * u.deg, frame=target)
		return [SkyCoord(row[0] * u.deg, row[1] * u.deg, frame=target) for row in out]

	def convertRaDec(self, ra, dec, forward:bool=True):
		'''
		Convert the provided ra,dec values using the mapping defined in this object from one frame to another.

		Bare numbers are read as degrees; Quantities are converted. When either
		input is a Quantity, the results carry units — degrees when the
		destination frame is a sky frame, or bare native values when it is not
		(e.g. pixel results are never labeled with any unit).

		:param ra: right ascension value(s) in degrees or astropy.units.Quantity
		:param dec: declination value(s) in degrees or astropy.units.Quantity
		:param forward: True: base -> current; False: current -> base
		'''
		quantity_input = isinstance(ra, Quantity) or isinstance(dec, Quantity)

		if quantity_input:
			# coerce each independently: bare numbers keep their documented
			# degrees meaning, so mixed Quantity/bare pairs stay legal
			ra = u.Quantity(ra, u.deg)
			dec = u.Quantity(dec, u.deg)
			pts = u.Quantity([ra, dec])
		else:
			pts = np.stack([np.asarray(ra, dtype=float), np.asarray(dec, dtype=float)])

		# scalars stack to rank-1 (2,), a single-point form where parallel_axes
		# must stay None (§2.2); arrays stack to (2, n), parallel by declaration
		out = self._tran(pts, forward, parallel_axes=(True if np.ndim(ra) > 0 else None))

		if np.ndim(ra) == 0:
			ra_out, dec_out = out          # single-point form -> (naxes,)
		else:
			ra_out, dec_out = out[:, 0], out[:, 1]

		if quantity_input:
			dst = self.astObject.getframe(Ast.CURRENT if forward else Ast.BASE)
			if bridge.is_sky(dst):
				return ra_out * u.deg, dec_out * u.deg
			# non-sky destination: bare native values — there is no honest unit
			# to attach (NEVER u.deg unconditionally; that was bug N2)
			return ra_out, dec_out
		return ra_out, dec_out

	def convert(self, points=None, forward:bool=True, *, parallel_axes:Optional[bool]=None):
		'''
		Deprecated alias for :meth:`convertPoints` (kept for one release).
		'''
		warnings.warn("ASTFrameSet.convert() is deprecated; use convertPoints() "
		              "(or convertRaDec() for separate ra/dec arrays).",
		              DeprecationWarning, stacklevel=2)
		return self.convertPoints(points, forward, parallel_axes=parallel_axes)

	# move to ASTMapping?
	def pix2world(self, points:Iterable, *, parallel_axes:Optional[bool]=None) -> np.ndarray:
		'''
		Convert provided coordinates from the pixel (base) frame to the world (current) frame.

		Points may be provided as coordinate pairs, e.g.

		.. code-block::

			[ [x1, y1], [x2, y2], ... ]           # shape (n, 2)

		or as parallel axis arrays, e.g.

		.. code-block::

			[ [x1, x2, ...], [y1, y2, ...] ]      # shape (2, n)

		A square ``(2, 2)`` array is read as PAIRS (DECISIONS D3); pass
		``parallel_axes=True`` for the parallel reading. A single point may also
		be specified alone, e.g. ``[a, b]``, and returns a single point.

		Returns shape ``(npoints, 2)`` — or ``(2,)`` for single-point input — in
		the current frame's units (degrees for sky frames, normalized).

		:param points: input coordinates in the base (pixel) frame
		:param parallel_axes: orientation override for square input (see above)
		'''
		return self._tran(points, True, parallel_axes=parallel_axes)

	# move to ASTMapping?
	def world2pix(self, points:Union[Iterable, SkyCoord], *, parallel_axes:Optional[bool]=None) -> np.ndarray:
		'''
		Convert provided coordinates from the world (current) frame to the pixel (base) frame.

		Points may be provided as coordinate pairs ``(n, 2)``, parallel axis
		arrays ``(2, n)``, a single point ``[a, b]``, SkyCoords (any system —
		converted to the world frame's system), or Quantities; bare values are
		read as degrees when the world frame is a sky frame. A square ``(2, 2)``
		array is read as PAIRS (DECISIONS D3); pass ``parallel_axes=True`` for
		the parallel reading.

		Returns shape ``(npoints, 2)`` — or ``(2,)`` for single-point-form
		input — in the base frame's native units.

		:param points: input coordinates in the current (world) frame
		:param parallel_axes: orientation override for square input (see above)
		'''
		return self._tran(points, False, parallel_axes=parallel_axes)
