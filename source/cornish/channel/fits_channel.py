#!/usr/bin/env python

import logging
from typing import Union

import numpy as np
import starlink
import starlink.Ast as Ast

from .ast_channel import ASTChannel
from ..mapping.frame import ASTFrame, ASTFrameSet
from ..region import ASTPolygon, ASTCircle
from ..exc import FrameNotFoundException, IncompleteHeader

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
	:param header: a FITS header, type `astropy.io.fits.header.Header` or `fitsio.fitslib.FITSHDR`, a list of tuples/arrays (keyword,value), a list of strings, or a single string
	'''
	# note: quote 'fitsio.hdu.base.HDUBase' as is it an optional package
	def __init__(self, ast_object:Ast.FitsChan=None, hdu:Union[astropy.io.fits.hdu.base.ExtensionHDU,'fitsio.hdu.base.HDUBase']=None, header=None):
		'''
		Initialize object with either an HDU or header from fitsio or astropy.io.fits.

		Parameters
		----------
		header : astropy.io.fits.header.Header or fitsio.fitslib.FITSHDR or list of tuples/arrays (keyword,value)

		:param header: FITS header as a dictionary of strings (keyword,value) or one of these types: astropy.io.fits.header.Header, fitsio.fitslib.FITSHDR, or a plain string divisible by 80 characters.
		'''
		self._frameSet = None

		# NOTE: `is not None`, not truthiness — an EMPTY Ast.FitsChan is falsy
		# (it exposes a card-count length), and truthiness here silently
		# discarded a provided empty channel (found by review round 3)
		if ast_object is not None:
			if hdu is not None or header is not None:
				raise ValueError("If 'ast_object' is provided, 'hdu' or 'header' should not be set.")
			if isinstance(ast_object, starlink.Ast.FitsChan):
				self._dimensions = None # pixel dimensions
				self.header = None
				super().__init__(ast_object=ast_object)
				return
			else:
				raise TypeError("Expected first parameter to be type 'ast_object'; did you mean to use 'header=' or 'hdu='?")

		# ----------
		# Validation
		# ----------
		# `is not None` throughout: an empty astropy Header (like an empty
		# FitsChan) is falsy, and truthiness would silently ignore it
		if hdu is not None and header is not None:
			raise ValueError("Only specify an HDU or header to create a ASTFITSChannel object.")

		# get header from HDU
		if hdu is not None:
			if _astropy_available and isinstance(hdu, (astropy.io.fits.hdu.image.PrimaryHDU, astropy.io.fits.hdu.base.ExtensionHDU)):
				header = hdu.header        # type: astropy.io.fits.header.Header
			elif _fitsio_available and isinstance(hdu, fitsio.hdu.base.HDUBase):
				header = hdu.read_header() # type: fitsio.fitslib.FITSHDR
			else:
				raise TypeError("Unknown HDU type specified ('{0}').".format(type(hdu)))
		# ----------

		self._dimensions = None # pixel dimensions
		self.astObject = Ast.FitsChan() # sink=self

		if header is not None:
			if isinstance(header, (list, str, dict)) and len(header) == 0:
				# an empty container can't describe a header — refuse loudly
				# rather than fall into the type probes below (an empty list
				# previously crashed with a raw IndexError)
				raise ValueError("An empty header was provided; pass header content or omit the parameter.")
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
				raise TypeError("Could not work with the type of header provided ('{0}').".format(type(header)))
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
			raise TypeError("Unable to parse the provided header.")

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
			raise ValueError("If a card is specified, don't set 'keyword', 'value', or 'comment'.")
		if not any([card, keyword, value, comment]):
			raise ValueError("Specify either a 'card' value or else 'keyword', 'value', or 'comment'.")

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
				raise TypeError("Unknown card type: '{0}'.".format(type(card)))

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
			raise ValueError("A keyword must be specified.")

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
			#logger.debug(self.astObject)
			self.astObject.clear("card") # move pointer to start of fitschan to read from there
			ast_frame_set = self.astObject.read() # deletes all WCS cards from the starlink.Ast.FitsChan object
			#self._readHeader()					  # reread for future use
			#raise Exception()
			if ast_frame_set is None:
				raise FrameNotFoundException("Could not create frame set from FITS header.")
			self._frameSet = ASTFrameSet(ast_object=ast_frame_set)

			# .. todo:: the result could be:
			# 	* NULL : wasn't able to create any AST object from the headers
			#	* some other AST object (other than frame set) -> very rare, but possible
			#		- there is a scheme where arbitrary AST objects can be stored in a FITS header
		return self._frameSet

	@property
	def dimensions(self):
		'''
		Returns pixel dimensions as a NumPy array.

		:raises cornish.exc.IncompleteHeader: when the NAXIS/NAXISn cards are missing
			(common in primary HDUs and some cutout services)
		'''
		if self._dimensions is None:
			dims = list()
			try:
				naxis = self.valueForKeyword("NAXIS")
			except KeyError as e:
				raise IncompleteHeader("The header does not contain an 'NAXIS' card, "
				                       "so the pixel dimensions cannot be determined.") from e
			if isinstance(naxis, str):
				naxis = int(naxis)
			for i in range(naxis):
				keyword = "NAXIS{0}".format(i+1)
				try:
					dims.append(int(self.valueForKeyword(keyword)))
				except KeyError as e:
					raise IncompleteHeader(f"The header declares NAXIS={naxis} but does not "
					                       f"contain an '{keyword}' card.") from e
			self._dimensions = np.array(dims)
		return self._dimensions

	# .. todo:: move to polygon class
	def boundingPolygon(self) -> ASTPolygon:
		'''
		Returns an ASTPolygon in the WCS sky frame that bounds the field described
		by this FITS header. (Must be a 2D image with WCS present.)

		Note (SPEC-04A M31): this method previously documented (but never
		delivered — it was dead code) a polygon in a flat degree frame; it now
		returns the field's bounding polygon in the sky frame, delegating to the
		same implementation as :meth:`ASTPolygon.fromFITSHeader`.
		'''
		return ASTPolygon._fromParsedWCS(self.frameSet, self.dimensions)

	def boundingCircle(self) -> ASTCircle:
		'''
		Returns the smallest circle that bounds the field described by this FITS header.

		The geometry is AST-exact (``astGetRegionDisc`` on the field's bounding
		polygon), replacing the earlier hand-rolled corner heuristics (SPEC-04A M32).
		'''
		if len(self.dimensions) != 2:
			raise ValueError("Requesting bounding circle on an HDU that is not 2D.")

		return self.boundingPolygon().boundingCircle()

