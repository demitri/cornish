#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)

import starlink.Ast as Ast

from .channel import Channel
from ..mapping.frame import ASTFrameSet

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

class FITSChannel(Channel):
	'''
	A representation of a FITS header. Use the property "astObject" for AST functions.
	self.astObject is of type starlink.Ast.FitsChan.
	'''
	def __init__(self, hdu=None, header=None):
		'''
		Initialize object with either an HDU or header from fitsio or astropy.io.fits.
		TODO: accept some form of text string.
		'''
		# internal AST object
		#self.fitsChan = None
		
		if all([hdu, header]):
			raise Exception("Only specify an HDU or header to create a FITSChannel object.")
		
		if hdu and _astropy_available and isinstance(hdu, astropy.io.fits.hdu.base.ExtensionHDU):
			header = hdu.header        # type: astropy.io.fits.header.Header
		elif hdu and _fitsio_available and isinstance(hdu, fitsio.fitslib.HDUBase):
			header = hdu.read_header() # type: fitsio.fitslib.FITSHDR
		
		if header is not None:
			# try to read the header from an Astropy header object
			if _astropy_available and isinstance(header, astropy.io.fits.header.Header):
				#self.fitsChan = Ast.FitsChan(source="".join([c.image for x in hdu.header.cards]))
				self.astObject = Ast.FitsChan()
				for card in hdu.header.cards:
					self.addHeader(card=card)
		
			# try to read the header from an fitsio header object
			elif _fitsio_available and isinstance(header, fitsio.fitslib.FITSHDR):
				all_cards = list()
				[all_cards.append(record["card_string"]) for record in header.records()]
				#print(all_cards)
				#self.fitsChan = Ast.FitsChan(source=all_cards) # don't know why this doesn't work
				self.astObject = Ast.FitsChan()
				for card in all_cards:
					#self.fitsChan[record["name"]] = record["value"]
					self.addHeader(card=card)
			
			else:
				raise Exception("Could not work with the type of HDU provided ('{0}').".format(type(hdu)))

		if self.astObject is None:
			self.astObject = Ast.FitsChan() # make an empty FITSChannel
		
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
				self.astObject.putfits(card)
				return
			
			else:
				raise Exception("Unknown card type: '{0}'.".format(type(card)))
		
		else:
			#
			if keyword in ["KEYWORD", "HISTORY", "COMMENT"]:
				raise Exception("The keywords 'KEYWORD', 'HISTORY', 'COMMENT' are not yet implemented.")
		
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
		return ASTFrameSet(ast_frame_set=self.astObject.read())
	
	@property
	def dimensions(self):
		''' Returns a list of dimensions. '''
		dims = list()
		for i in range(self.valueForKeyword("NAXIS")):
			dims.append(self.valueForKeyword("NAXIS{0}".format(i+1)))
		return dims
	
		
