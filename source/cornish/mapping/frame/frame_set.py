#/usr/bin/env python

import logging
from typing import Union, Iterable

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
#from ...channel import ASTFITSChannel

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
				Exception("ASTFrameSet: Unhandled ast_object type ('{0}')".format(ast_object))
		else:
			# construct from provided base_frame
			if isinstance(base_frame, starlink.Ast.Frame):
				fs = Ast.FrameSet(frame=base_frame, options=None)
			elif isinstance(base_frame, ASTFrame):
				fs = Ast.FrameSet(frame=base_frame.astObject, options=None)
			else:
				Exception("ASTFrameSet: Unhandled base_frame type ('{0}')".format(base_frame))

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
			elif isinstance(frame, starlink.Ast.AstFrame):
				ast_frames.append(frame)
			else:
				raise ValueError(f"The provided frames must either be of type ASTFrame or starlink.Ast.Frame (got '{type(frame)}'.")

		frame_set = Ast.Frame.convert(ast_frames[0], ast_frames[1])
		if frame_set is None:
			# I don't know if this is actually the failure mode
			raise Exception("An ASTFrameSet could not be created since the mapping between the two provided frames could not be determined (conversion may not be possible).")
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
			raise Exception(f"frame must be provided to 'addFrame'")

		iframe = Ast.BASE # add to base frame
		map = None
		self.astObject.addframe(iframe, map, frame)

	@property
	def mapping(self):
		'''
		Return an object that maps points from the base frame to the current frame of this frame set.
		'''
		return ASTMapping(ast_object=self.astObject.getmapping()) # default is from base to current

	def convertPoints(self, points:Iterable, forward:bool=True):
		'''
		Convert the provided coordinate points using the mapping defined in this object from one frame to another.

		:param points:
		:param forward:
		'''
		if isinstance(points[0], SkyCoord):
			# handle SkyCoord objects
			ra_rad = [c.ra.to(u.rad).value for c in points]
			dec_rad = [c.dec.to(u.rad).value for c in points]
			out = self.astObject.tran(np.array([ra_rad, dec_rad]), forward)
			return_points = list()
			for p in out.T:
				return_points.append(  SkyCoord(ra=p[0]*u.rad, dec=p[1]*u.rad)  )
			return return_points

		else:
			# handle arrays

			if points.shape[0] == 2:
				#
				# check the form [ [ra1, ra2, ...], [dec1, dec2, ...] ]

				#if len(points[0]) == 1:
				#	raise ValueError("")
				ra_rad = np.deg2rad(points[0])
				dec_rad = np.deg2rad(points[1])
				out = self.astObject.tran(np.array([ra_rad, dec_rad]), forward)
				return np.rad2deg(out)
			elif points.shape[1] == 2:
				#
				# check the form [ [ra1, dec1], [ra2, dec2], ... ]

				points = points.T
				ra_rad = np.deg2rad(points[0])
				dec_rad = np.deg2rad(points[1])
				out = self.astObject.tran(np.array([ra_rad, dec_rad]), forward)
				return np.rad2deg(out.T)



		# elif isinstance(points, np.ndarray):
		# 	raise NotImplementedError("handle numpy arrays")
		# else:
		# 	# handle other iterables
		# 	#TODO validate inputs for error reporting
		# 	assert len(x_coordinates) == len(y_coordinates), "Coordinate arrays of unequal lengths."
		# 	#forward = True # convert from frame1 to frame2
		# 	#in_coords = np.array([x_coordinates, y_coordinates]).T
		# 	in_coords = np.array(list(zip(x_coordinates, y_coordinates))).T
		# 	print(in_coords)
		# 	forward = True
		# 	return self.astObject.tran(in_coords, forward)

	def convertRaDec(self, ra, dec, forward:bool=True):
		'''
		Convert the provided ra,dec values using the mapping defined in this object from one frame to another.

		:param ra: right ascension value in degrees or astropy.units.Quantity
		:param dec: declination value in degrees or astropy.units.Quantity
		:param forward:
		'''

		provided_quantity = False # return the same type as provided
		if isinstance(ra, Quantity):
			ra = ra.to(u.deg).value
			provided_quantity = True
		if isinstance(dec, Quantity):
			dec = dec.to(u.deg).value
			provided_quantity = True

		try:
			# when ra,dec single values
			out = self.astObject.tran(np.array(np.deg2rad([[ra], [dec]])), forward)
		except ValueError as e:
			# when ra, dec are already arrays
			out = self.astObject.tran(np.array(np.deg2rad([ra, dec])), forward)

		if provided_quantity:
			return out[0]*u.deg, out[1]*u.deg
		else:
			return np.rad2deg(out[0]), np.rad2deg(out[1])


	def convert(self, points=None, ra:Iterable=None, dec:Iterable=None, forward:bool=True): #x_coordinates=None, y_coordinates=None):
		'''
		Convert the provided points using the mapping defined in this object from one frame to another.

		If no units are provided (e.g. a bare NumPy array), degrees are assumed.

		:param points:
		:param ra:
		:param dec:
		:param forward:
		'''

		# # parameter validation
		# if all([x is None in [points, ra, dec]):
		# 	raise ValueError(f"Either 'points' must be provided, or both 'ra' and 'dec' parameters.")
#
		# if points is None:
		# 	# specify ra & dec?
		# 	if any([x is None in [ra,dec]]):
		# 		raise ValueError("If one of 'ra' or 'dec' is specified, both must be provided.")
		# 	elif all(
		# else:
		# 	# specifying points?
		# 	if any([x is not None for x in [ra,dec]]):
		# 		raise ValueError("'ra' and 'dec' should not be provided")


		# handle SkyCoords
		if isinstance(points, SkyCoord) or \
		  (isinstance(points, list) and len(points) > 0 and isinstance(points[0], SkyCoord)):
			if isinstance(points, SkyCoord):
				p = points
				out = self.astObject.tran(np.array([[p.ra.to(u.rad).value], [p.dec.to(u.rad).value]]), forward)
				return SkyCoord(ra=out[0]*u.rad, dec=out[1]*u.rad)
			else:
				raise NotImplementedError("Need to implement arrays of SkyCoord objects.")
		elif isinstance(points, np.ndarray):
			raise NotImplementedError("handle numpy arrays")
		else:
			# handle other iterables
			#TODO validate inputs for error reporting
			assert len(x_coordinates) == len(y_coordinates), "Coordinate arrays of unequal lengths."
			#forward = True # convert from frame1 to frame2
			#in_coords = np.array([x_coordinates, y_coordinates]).T
			in_coords = np.array(list(zip(x_coordinates, y_coordinates))).T
			print(in_coords)
			forward = True
			return self.astObject.tran(in_coords, forward)

	# move to ASTMapping?
	def pix2world(self, points:Iterable) -> np.ndarray:
		'''
		Convert provided coordinates from a world frame to a pixel frame.

		This method will throw a ``cornish.exc.FrameNotAvailable`` exception if
		the frame set does not contain both a pixel and world frame.

		Format of points:

		.. code-block::

			[ [ values on axis 1 ], [ values on axis 2 ], ... ]

		e.g. pixel to sky:

		.. code-block::

			[ [x1, x2, ...], [y1, y2, ...] ]

		A single point may also be specified alone, e.g. ``[a,b]`` or ``np.array([a,b])``.

		:param points: input list of coordinates as numpy.ndarray, 2-dimensional array with shape (2,npoint)
		'''
		pix2world = True # AST flag ("forward" parameter)

		if isinstance(points, (list, tuple)):
			points = np.array(points)

		# check for single point, e.g. pix2world([12, 33])
		single_point = False
		if points.shape == np.array([1,2]).shape:
			single_point = True
			points = np.array([points])

		# Note that tran() takes points in the form [ [ra2, ra2, ...], [dec1, dec2, ...] ]
		out = np.rad2deg(self.astObject.norm((self.astObject.tran(points.T, pix2world)))).T

		if single_point:
			return out[0]
		else:
			return out

	# move to ASTMapping?
	def world2pix(self, points:Union[Iterable, SkyCoord]) -> np.ndarray:
		'''
		Convert provided coordinates from a world frame to a pixel frame.

		This method will throw a ``cornish.exc.FrameNotAvailable`` exception if
		the frame set does not contain both a pixel and world frame.

		Points must have the shape (2,n), e.g.:

		.. code-block::

			[ [ra1, ra2, ...], [dec1, dec2, ...] ]

		A single point may also be specified alone, e.g. ``[a,b]`` or ``np.array([a,b])``.

		:param points: input list of coordinates as numpy.ndarray, 2-dimensional array with shape (2,npoints); units are assumed to be degrees if not specified via ``astropy.units.Quantity``
		'''
		world2pix = False # AST flag ("forward" parameter)
		single_point = False

		if isinstance(points, (list, tuple)):
			points = np.array(points)
		elif isinstance(points, astropy.coordinates.SkyCoord):
			points = np.array([[points.ra.to(u.rad).value, points.dec.to(u.rad).value]]) * u.rad
			single_point = True

		# check for single point, e.g. world2pix([12.34, 56.78])
		if points.shape == np.array([1,2]).shape:
			single_point = True
			points = np.array([points])

		# AST expects radians
		if isinstance(points, astropy.units.quantity.Quantity):
			points_rad = points.to(u.rad).value
		else:
			points_rad = np.deg2rad(points)

		# Note that tran() takes points in the form [ [ra2, ra2, ...], [dec1, dec2, ...] ]
		out = self.astObject.tran(points_rad.T, world2pix).T
		if single_point:
			return out[0]
		else:
			return out

	#def frame(self, name=None):
	#	'''
	#	Extract a frame from the frame set with the given name.
	#	'''
	#	self.


'''
Transform the coordinates of a set of points provided according the mapping defined by this object.

:param in: input list of coordinates as numpy.ndarray,
		   any iterable list accepted
		   2-dimensional array with shape (nin,npoint)
:param out: output coordinates

Format of points:

.. code-block::

	[ [ values on axis 1 ], [ values on axis 2 ], ... ]

e.g. sky to pixel:

.. code-block::

	[ [ra1, ra2, ...], [dec1, dec2, ...] ]

'''
