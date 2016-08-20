#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

import ast
import numpy as np
import starlink.Ast as Ast

from .ast_channel import ASTChannel
from ..mapping.frame import ASTFrame, ASTFrameSet
from ..region import ASTBox, ASTPolygon

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

class ASTFITSChannel(ASTChannel):
	'''
	A representation of a FITS header. Use the property "astObject" for AST functions.
	self.astObject is of type starlink.Ast.FitsChan.
	'''
	def __init__(self, hdu=None, header=None):
		'''
		Initialize object with either an HDU or header from fitsio or astropy.io.fits.
		@param header FITS header as an array of strings or one of these types: astropy.io.fits.header.Header, fitsio.fitslib.FITSHDR
		TODO: accept some form of text string.
		'''
		# defines internal AST object
		super(ASTFITSChannel, self).__init__()
		# super().__init__() # Python 3
		
		self._frameSet = None
		
		# ----------
		# Validation
		# ----------
		if all([hdu, header]):
			raise Exception("Only specify an HDU or header to create a ASTFITSChannel object.")
		
		if hdu and _astropy_available and isinstance(hdu, astropy.io.fits.hdu.base.ExtensionHDU):
			header = hdu.header        # type: astropy.io.fits.header.Header
		elif hdu and _fitsio_available and isinstance(hdu, fitsio.fitslib.HDUBase):
			header = hdu.read_header() # type: fitsio.fitslib.FITSHDR
		# ----------
		
		self._dimensions = None
		self.astObject = Ast.FitsChan()
		
		if header is not None:
			# Note that the starlink.Ast.Chanel.read() operation is destructive.
			# Save the header so it can be restored/reused.
			self.header = header
			if isinstance(header, (astropy.io.fits.header.Header, fitsio.fitslib.FITSHDR)) or \
			   isinstance(header, dict) or \
			   (isinstance(header, list) and isinstance(header[0], str)):
			   self._readHeader()
			else:
				raise Exception("Could not work with the type of header provided ('{0}').".format(type(header)))
		else:
			pass # work with an empty Ast.FitsChan()
	
	def _readHeader(self):
		'''
		Internal method to read a FITS header.
		Accepts astropy.io.fits.header.Header, fitsio.fitslib.FITSHDR, a dictionary (keyword,value), or a list of strings.
		'''
		print("ASTFitsChannel._reading header.")

		assert self.header is not None, "Attempting to read a header before it has been set."
		
		# try to read the header from an Astropy header object
		if _astropy_available and isinstance(self.header, astropy.io.fits.header.Header):
			#self.fitsChan = Ast.FitsChan(source="".join([c.image for x in hdu.header.cards]))
			for card in self.header.cards:
				self.addHeader(card=card)
	
		# try to read the header from an fitsio header object
		elif _fitsio_available and isinstance(self.header, fitsio.fitslib.FITSHDR):
			all_cards = list()
			[all_cards.append(record["card_string"]) for record in self.header.records()]
			#self.fitsChan = Ast.FitsChan(source=all_cards) # don't know why this doesn't work
			for card in all_cards:
				self.addHeader(card=card)
		
		elif isinstance(self.header, dict):
			for keyword in self.header.keys():
				self.astObject[keyword] = self.header[keyword]
		
		elif isinstance(self.header, list) and isinstance(self.header[0], str):
			for card in self.header:
				self.addHeader(card=card)
	
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
	
	def addHeader(self, card=None, keyword=None, value=None, comment=None):
		'''
		Add a card to this header.
		
		Currently values are overwritten if they exist in the header.
		@param card Can be of type: fitsio.fitslib.FITSRecord, astropy.io.fits.card.Card, or a text string
		'''
		
		#raise Exception("putfits doesn't work")
		
		overwrite = True # can implement option later
		# NOTE: putfits(card, overwrite=<anything>) results in
		#       TypeError: putfits() takes no keyword arguments
		#       which is not what the documentation says.
		
		if card and any([keyword, value, comment]):
			raise Exception("If a card is specified, don't set 'keyword', 'value', or 'comment'.")
		if not any([card, keyword, value, comment]):
			raise Exception("Specify either a 'card' value or else 'keyword', 'value', or 'comment'.")
				
		if card:
			if _fitsio_available and isinstance(card, fitsio.fitslib.FITSRecord):
				self.astObject.putfits(card["card_string"])
				return
			
			elif _astropy_available and isinstance(card, astropy.io.fits.card.Card):
				self.astObject.putfits(card.image)
				return
			
			elif isinstance(card, str) and len(card) <= 80:
				# parse individual strings
				if "=" in card:
					keyword, string_value = card.split("=")
					keyword = keyword.strip()
					string_value = string_value.strip()
					try:
						# handles (float, int), converted to the correct type
						value = ast.literal_eval(string_value)
					except ValueError:
						# string value
						if s == "T":
							value = True
						elif s == "F":
							value = False
						else:
							# a regular string - remove leading and trailing quotes if present
							if s[0] == "'" and s[-1] == "'":
								value = s[1:-1]
							else:
								value = s
					self.astObject[keyword] = value
				return
			
			else:
				raise Exception("Unknown card type: '{0}'.".format(type(card)))
		
		else:
			#
			if keyword in ["HISTORY", "COMMENT", ""]:
				raise Exception("The keywords 'HISTORY', 'COMMENT', or blank ('') are not yet implemented.")
		
			# reconstruct from keyword, value, comment
			if type(value, str):
				if comment:
					card = "{0:8}= '{1}' / {2}".format(keyword, value, comment)
				else:
					card = "{0:8}= '{1}'".format(keyword, value)
			
			elif type(value) in [int, float]:
				if comment:
					card = "{0:8}= {1:-20} / {2}".format(keyword, value, comment)
				else:
					card = "{0:8}= {1:-20}".format(keyword, value)
			
			else:
				raise Exception("Not implemented: handle value of type '{0}'.".format(type(value)))
		
		assert len(card) <= 80, "The FITS card is incorrectly formatted (or a bad value has been given)."
		
		self.astObject.putfits(card=card, overwrite=True)
		
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
		
	def frameSet(self):
		'''
		Reads the FITS header and convert it into an AST frame set.
		'''
		# NOTE! read() is a destructive operation! This logic will have to be changed
		#       if is needs to be used again anywhere else in this class (whic is why the header is being saved).
		if self._frameSet is None:
			ast_frame_set = self.astObject.read() # deletes all WCS cards from the starlink.Ast.FitsChan object
			#self._readHeader()					  # reread for future use
			#raise Exception()
			self._frameSet =  ASTFrameSet(ast_frame_set=ast_frame_set)
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
			for i in range(naxis):
				dims.append(self.valueForKeyword("NAXIS{0}".format(i+1)))
			self._dimensions = np.array(dims)
		return self._dimensions
	
	def boundingPolygon(self):
		'''
		Returns an ASTPolygon that bounds the field described by the FITS header. (Must be a 2D image with WCS present.)
		
		'''
		
		# The code below is adapted from code originally provided by David Berry.
		
		print("A")
		
		# create an ASTFrameSet that contains two frames (pixel grid, WCS) and the mapping between them
		wcsFrameSet = self.frameSet()
				
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
		pixelbox = ASTBox(frame=wcsFrameSet.baseFrame(),
						  cornerPoint=[0.5,0.5], # center of lower left pixel
						  cornerPoint2=[dims[0]+0.5, dims[1]+0.5])
		
		print("B")

		#  Map this box into (RA,Dec)
		#
		skybox = pixelbox.regionWithMapping(map=wcsFrameSet, frame=wcsFrameSet) # -> ASTRegion
		raise Exception()
		print("C")

		#  Get the (RA,Dec) at a large number of points evenly distributed around
		#  the polygon. The number of points created is controlled by the
		#  MeshSize attribute of the polygon.
		#
		mesh = skybox.boundaryPointMesh() # np.array of points

		print("D")

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

		print("E")

		#  Remove mesh points where the polygon is close to a Cartesian straight
		#  line, and retain them where it deviates from a stright line, in order
		#  to achieve an max error of 1 arc-second (4.8E-6 rads).
		#
		return flatpoly.downsize(maxerr=4.848e-6) # -> ASTPolygon 

