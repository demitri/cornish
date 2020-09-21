#!/usr/bin/env python

import math
import logging
from typing import Union

#import ast
import numpy as np
import starlink
import starlink.Ast as Ast

from .ast_channel import ASTChannel
from ..mapping.frame import ASTFrame, ASTFrameSet
from ..region import ASTBox, ASTPolygon, ASTCircle
from ..exc import FrameNotFoundException

_astropy_available = True
_fitsio_available = True
try:
	import astropy
	import astropy.io.fits
except ImportError:
	_astropy_available = False
try:
	import fitsio
except ImportError:
	_fitsio_available = False

__all__ = ['ASTFITSChannel']

logger = logging.getLogger("cornish") # cornish logger

class ASTFITSChannel(ASTChannel):
	'''
	A representation of a FITS header. Use the property "astObject" for AST functions.
	self.astObject is of type starlink.Ast.FitsChan.
	
	:param ast_object: an existing `starlink.Ast.FitsChan` object
	:param hdu: a FITS HDU of type `astropy.io.fits.hdu.base.ExtensionHDU` or `fitsio.fitslib.HDUBase`
	:param header: a FITS header, type `astropy.io.fits.header.Header` or `fitsio.fitslib.FITSHDR`, a list of tuples/arrays (keyword,value), or a single string
	'''
	# note: quote 'fitsio.hdu.base.HDUBase' as is it an optional package
	def __init__(self, ast_object:Ast.FitsChan=None, hdu:Union[astropy.io.fits.hdu.base.ExtensionHDU,'fitsio.hdu.base.HDUBase']=None, header=None):
		'''
		Initialize object with either an HDU or header from fitsio or astropy.io.fits.
		
		Parameters
		----------
		header : astropy.io.fits.header.Header or fitsio.fitslib.FITSHDR or list of tuples/arrays (keyword,value)
		
		@param header FITS header as a dictionary of strings (keyword,value) or one of these types: astropy.io.fits.header.Header, fitsio.fitslib.FITSHDR, or a plain string divisible by 80 characters.
		'''
		self._frameSet = None

		if ast_object:
			if any([hdu, header]):
				raise ValueError("If 'ast_object' is provided, 'hdu' or 'header' should not be set.")
			else:
				if isinstance(ast_object, starlink.Ast.FitsChan):
					super().__init__(ast_object=ast_object)
				else:
					raise ValueError("Expected first parameter to be type 'ast_object'; did you mean to use 'header=' or 'hdu='?")
		
		# ----------
		# Validation
		# ----------
		if all([hdu, header]):
			raise Exception("Only specify an HDU or header to create a ASTFITSChannel object.")
			
		# get header from HDU
		if hdu:
			if hdu and _astropy_available and isinstance(hdu, (astropy.io.fits.hdu.image.PrimaryHDU, astropy.io.fits.hdu.base.ExtensionHDU)):
				header = hdu.header        # type: astropy.io.fits.header.Header
			elif hdu and _fitsio_available and isinstance(hdu, fitsio.hdu.base.HDUBase):
				header = hdu.read_header() # type: fitsio.fitslib.FITSHDR
			else:
				raise Exception(f"ASTFITSChannel: unknown HDU type specified ('{0}').".format(type(hdu)))
		# ----------
		
		self._dimensions = None
		self.astObject = Ast.FitsChan() # sink=self
		
		if header is not None:
			# Note that the starlink.Ast.Channel.read() operation is destructive.
			# Save the header so it can be restored/reused.
			self.header = header
			if isinstance(header, astropy.io.fits.header.Header) or \
			   (_fitsio_available and isinstance(header, fitsio.fitslib.FITSHDR)) or \
			   (isinstance(header, list) and isinstance(header[0], str)):
			   self._readHeader()
			elif isinstance(header, list):
				logger.debug("Found list header")
				self._readHeader()
			elif isinstance(header, str) and (len(header) % 80 == 0):
				# split into list, then process
				self.header = [header[i:i+80] for i in range(0, len(header), 80)]
				self._readHeader()
			elif isinstance(header, dict):
				for keyword in header:
					self.astObject[keyword] = header[keyword]
					#self.addHeader(keyword=key, value=header[key])
				#self._readHeader()
			else:
				raise Exception("Could not work with the type of header provided ('{0}').".format(type(header)))
		else:
			pass # work with an empty Ast.FitsChan()
			logger.warning("ASTFITSChannel: no data found: working with an empty FITSChannel.")
	
#	def astsink(self, header):
#		logger.debug("ASTFITSChan._sink: header = {0}".format(header))
	
	def _readHeader(self):
		'''
		Internal method to read a FITS header from 'self.header'.
		Accepts astropy.io.fits.header.Header, fitsio.fitslib.FITSHDR, a dictionary (keyword,value), a list of strings,
		or a plain string divisible by 80 characters (full header as a single string).
		'''
		assert self.header is not None, "Attempting to read a header before it has been set."
		
		# try to read the header from an Astropy header object
		if _astropy_available and isinstance(self.header, astropy.io.fits.header.Header):
			#self.fitsChan = Ast.FitsChan(source="".join([c.image for x in hdu.header.cards]))
			for card in self.header.cards:
				self.addHeader(card=card)
	
		# try to read the header from an fitsio header object
		elif _fitsio_available and isinstance(self.header, fitsio.fitslib.FITSHDR):
			all_cards = list()
			all_cards.extend([record["card_string"] for record in self.header.records()])
			#self.fitsChan = Ast.FitsChan(source=all_cards) # don't know why this doesn't work
			for card in all_cards:
				self.addHeader(card=card)
		
		#elif isinstance(self.header, dict):
		#	for keyword in self.header.keys():
		#		self.astObject[keyword] = self.header[keyword]
		
		elif isinstance(self.header, list) and len(self.header) > 0:
			#
			# read header from provided list of strings or keyword/value pairs
			#
			first_item = self.header[0]
			if isinstance(first_item, str):
				for card in self.header:
					self.addHeader(card=card)
			elif isinstance(first_item, (list,tuple)) and len(first_item) == 2:
				for keyword,value in self.header:
					self.addHeader(keyword=keyword, value=value)
		
		else:
			raise Exception("ASTFITSChannel: unable to parse header.")
	
	def cardForKeyword(self, keyword=None):
		'''
		Get FITS card for given keyword.
		'''
		return self.astObject.findfits(keyword, inc=False)
	
	def clearAllCards(self):
		'''
		Delete all cards.
		'''
		self.astObject.emptyfits()
	
	def addHeader(self, card=None, keyword=None, value=None, comment=None, overwrite:bool=False):
		'''
		Add a card to this header.
		
		:param card: fitsio.fitslib.FITSRecord or astropy.io.fits.card.Card or str. A single FITS header.
		:param overwrite: If 'True', the provided card overwrites any existing one with the same keyword.
		'''
				
		# NOTE: putfits(card, overwrite=<anything>) results in
		#       TypeError: putfits() takes no keyword arguments
		#       which is not what the documentation says.
		
		if card and any([keyword, value, comment]):
			raise Exception("If a card is specified, don't set 'keyword', 'value', or 'comment'.")
		if not any([card, keyword, value, comment]):
			raise Exception("Specify either a 'card' value or else 'keyword', 'value', or 'comment'.")
				
		fits_header_card = None
		
		# This next block parses the input and defines 'fits_header_card', the value to be
		# passed to 'putfits()'.
		
		if card:
			#
			# Parse provided FITS header
			#
			#fitsio.fitslib.FITSRecord doesn't exist anymore?
			#if _fitsio_available and isinstance(card, fitsio.fitslib.FITSRecord):
			#	#self.astObject.putfits(card["card_string"])
			#	fits_header_card = card["card_string"]
			#	#return
			
			if _astropy_available and isinstance(card, astropy.io.fits.card.Card):
				#self.astObject.putfits(card.image)
				if len(card.image) > 80 and "&'CONTINUE" in card.image:
					# Astropy will concatenate long strings into one card, e.g.
					#    "SEXCAT  = '/home/eisa2/fltops/pipe/01-vsn/03000-MISDR1_24278_0266/d/00-visits/&'CONTINUE  '0001-img/07-try/MISDR1_24278_0266_0001-nd-cat.fits'
					# Split them back up again.
					cards = card.image.split("CONTINUE")
					self.astObject.putfits(cards.pop(0), overwrite)
					for c in cards:
						self.astObject.putfits("CONTINUE"+c, overwrite)
					self.astObject.retainfits()
					return
				else:
					fits_header_card = card.image
					#return
			
			elif isinstance(card, str):
				if len(card) <= 80:
					fits_header_card = card
				else:
					raise NotImplementedError("handle long cards")
			else:
				raise Exception("Unknown card type: '{0}'.".format(type(card)))
		
		else:
			#
			# Parse keyword/value/comment
			#
			# reconstruct from keyword, value, comment
			self.astObject[keyword] = value
			return
						
			# THE BLOCK BELOW DOES NOT WORK.
			# At the very least, data types are not being handled correctly.
			# The best plan is to extract the cards close to the source.
			#
			# not sure why types were introduced
			
#			if False:
#				raise Exception("Don't use this section of code without rigorous testing - it's critical to get the data type correct.")
#				
#				if keyword in ["HISTORY", "COMMENT", ""]:
#					raise Exception("The keywords 'HISTORY', 'COMMENT', or blank ('') are not yet implemented.")
#				
#				if keyword.startswith('NAXIS'):
#					value = int(value)
#				elif keyword.startswith('CRVAL2'):
#					value = float(value)
#				elif keyword.startswith('CRPIX'):
#					value = float(value)
#				
#				# reconstruct from keyword, value, comment
#				if isinstance(value, str):
#					if comment:
#						fits_header_card = "{0:8}= '{1}' / {2}".format(keyword, value, comment)
#					else:
#						fits_header_card = "{0:8}= '{1}'".format(keyword, value)
#				
#				elif isinstance(value, (int, float)):
#					if comment:
#						fits_header_card = "{0:8}= {1:-20} / {2}".format(keyword, value, comment)
#					else:
#						fits_header_card = "{0:8}= {1:-20}".format(keyword, value)
#				
#				else:
#					raise Exception("Not implemented: handle value of type '{0}'.".format(type(value)))
		
		assert len(fits_header_card) <= 80, "The FITS card is incorrectly formatted (>80 char) or a bad value has been given: \"{0}\"".format(fits_header_card)
		
		self.astObject.putfits(fits_header_card, overwrite)
		self.astObject.retainfits() # this a single flag -> turn into a top level property
		
	def valueForKeyword(self, keyword=None):
		'''
		Value from header for given keyword.
		'''
		if keyword is None:
			raise Exception("A keyword must be specified.")
		
		return self.astObject[keyword]
	
	#@property
	#def allCards(self):
	#	'''
	#	Return all cards for this FITS header.
	#	'''
	#	return self.astObject.writefits() # nope, it doesn't work this way
	
	def purgeWCS(self):
		'''
		Remove all keywords that describe WCS information.
		'''
		self.astObject.purgewcs()
	
	@property
	def frameSet(self):
		'''
		Read the FITS header and convert it into an AST frame set.
		'''
		# NOTE! read() is a destructive operation! This logic will have to be changed
		#       if it needs to be used again anywhere else in this class (which is why the header is being saved).
		if self._frameSet is None:
			self.astObject.clear("card") # move pointer to start of fitschan to read from there
			ast_frame_set = self.astObject.read() # deletes all WCS cards from the starlink.Ast.FitsChan object
			#self._readHeader()					  # reread for future use
			#raise Exception()
			if ast_frame_set is None:
				raise FrameNotFoundException("Could not create frame set from FITS header.")
			self._frameSet =  ASTFrameSet(ast_object=ast_frame_set)
			
			# .. todo:: the result could be:
			# 	* NULL : wasn't able to create any AST object from the headers
			#	* some other AST object (other than frame set) -> very rare, but possible
			#		- there is a scheme where arbitrary AST objects can be stored in a FITS header
		return self._frameSet
	
	@property
	def dimensions(self):
		''' Returns a list of dimensions as a NumPy array. '''
		if self._dimensions is None:
			dims = list()
			try:
				naxis = self.valueForKeyword("NAXIS")
			except TypeError as e:
				# TypeError: 'int' object is not subscriptable
				if "'int' object is not subscriptable" in str(e):
					raise Exception()
			if isinstance(naxis, str):
				naxis = int(naxis)
			for i in range(naxis):
				dims.append(int(self.valueForKeyword("NAXIS{0}".format(i+1))))
			self._dimensions = np.array(dims)
		return self._dimensions
	
	# .. todo:: move to polygon class
	def boundingPolygon(self):
		'''
		Returns an ASTPolygon that bounds the field described by the FITS header. (Must be a 2D image with WCS present.)
		
		'''
		
		# The code below is adapted from code originally provided by David Berry.
		
		# create an ASTFrameSet that contains two frames (pixel grid, WCS) and the mapping between them
		wcsFrameSet = self.frameSet
				
		# Create a Box describing the extent of the image in pixel coordinates.
		#
		# From David Berry:
		#       "Because of the FITS-WCS standard, the base Frame in a FrameSet read
		#       from a FITS header will always represent FITS pixel coordinates, which
		#       are defined by the FITS-WCS standard so that the bottom left (i.e.
		#       first) pixel in a 2D array is centred at (1,1). That means that
		#       (0.5,0.5) is the bottom left corner of the bottom left pixel, and
		#       (dim1+0.5,dim2+0.5) is the top right corner of the top right pixel.
		#       This results in the Box covering the whole image area."
		#
		dims = self.dimensions
		pixelbox = ASTBox(frame=wcsFrameSet.baseFrame,
						  cornerPoint=[0.5,0.5], # center of lower left pixel
						  cornerPoint2=[dims[0]+0.5, dims[1]+0.5])
		
		#  Map this box into (RA,Dec)
		#
		skybox = pixelbox.regionWithMapping(map=wcsFrameSet, frame=wcsFrameSet) # -> ASTRegion

		#  Get the (RA,Dec) at a large number of points evenly distributed around
		#  the polygon. The number of points created is controlled by the
		#  MeshSize attribute of the polygon.
		#
		mesh = skybox.boundaryPointMesh() # np.array of points

		#  Create a polygon using the vertices in the mesh. This polygon is
		#  defined in a basic Frame (flat geometry) - not a SkyFrame (spherical
		#  geometry). If we used a SkyFrame, then all the mesh points along each
		#  edge of the box would fall exactly on a geodesic (i.e. a great circle),
		#  and so the subsequent call to the downsize function would remove them all
		#  (except the corners). Using a basic Frame means that the downsize function
		#  will use geodesics that are Cartesian straight lines. So points that
		#  deviate by more than the required error form a Cartesian straight line
		#  will be retained by the downsize function.
		#
		degFlatFrame = ASTFrame(naxes=2)
		degFlatFrame.setUnitForAxis(axis=1, unit="deg")
		degFlatFrame.setUnitForAxis(axis=2, unit="deg")

		flatpoly = ASTPolygon(frame=degFlatFrame, points=mesh)

		#  Remove mesh points where the polygon is close to a Cartesian straight
		#  line, and retain them where it deviates from a stright line, in order
		#  to achieve an max error of 1 arc-second (4.8E-6 rads).
		#
		return flatpoly.downsize(maxerr=4.848e-6) # -> ASTPolygon 

	def boundingCircle(self):
		'''
		Returns the smallest circle that bounds the HDU represented by this FITS header.
		
		It is up to the called to know that this is a 2D image (only minimal checks are made).
		'''
		
		# contains two frames (pixels grid, WCS) and mapping between them
		#    - baseFrame    - native coordinate system (pixels)
		#    - currentFrame - WCS (ra, dec)
		wcsFrameSet = self.frameSet
		baseFrame = wcsFrameSet.baseFrame   # pixel coordinates frame
		wcsFrame = wcsFrameSet.currentFrame # (AST SkyFrame)
		
		#print(baseFrame.astObject)
		#print(wcsFrame.astObject)
		
		if len(self.dimensions) != 2:
			raise Exception("ASTFITSChannel: Requesting bounding circle on an HDU that is not 2D.")
	
		#logger.debug("dimensions: {0}".format(dims))
	
		#dims = self.dimensions
		# pixelbox = ASTBox(frame=wcsFrameSet.baseFrame,
		# 				  cornerPoint=[0.5,0.5], # center of lower left pixel
		# 				  cornerPoint2=[dims[0]+0.5, dims[1]+0.5])
	
		# the base frame (pixels) is used since we're using the pixel dimensions to define the area
		pixelbox = ASTBox(frame=baseFrame, dimensions=self.dimensions)
	
		# Use frame set to map between pixels and WCS coords (SkyFrame).
		# Result coordinates are degrees on sky.
		#
		corner_points = np.array(pixelbox.corners(mapping=wcsFrameSet)) # -> unit: degrees
	
		logger.debug("corner points (deg): {0}".format(corner_points))
		logger.debug("corner points (rad): {0}".format(np.deg2rad(corner_points)))
	
		# convert back to rad to use in AST functions
		corner_points = np.deg2rad(corner_points) # -> unit: radians
	
		# calculate the distance between two corner points
		c1 = corner_points[0]
		c2 = corner_points[2]
		diagonal_distance = wcsFrame.distance(c1, c2)
	
		ed = math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2) # distance in Euclidean space
		logger.debug("diagonal distance: {0} (rad)".format(diagonal_distance))
		logger.debug("Euclidian value:   {0} ".format(ed)) # should not match the value above
	
		# find point halfway on diagonal line
		center = wcsFrame.astObject.offset(c1, c2, diagonal_distance/2.0) # -> unit: radians
	
		#print("center: ", center)
		#print("center: ", np.rad2deg(center))
	
		# Use this point as the circle center.
		# Calculate the distance from the center to each corner.
		distances = [wcsFrame.distance(center, p) for p in corner_points]
	
		logger.debug("distances from center to each corner (rad): {}".format(distances))
		logger.debug("distances from center to each corner (deg): {}".format(np.rad2deg(distances)))
	
		radius = max(distances) # -> unit: radians
		radius *= 1.000001 # increase by a tiny amount to account for possible rounding errors
		
	
		logger.debug("radius: {0} (radians), {1} (deg)".format(radius, np.rad2deg(radius)))
		
		return ASTCircle(frame=wcsFrame, center=np.rad2deg(center), radius=np.rad2deg(radius))

