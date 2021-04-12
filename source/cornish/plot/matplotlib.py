
from typing import Union, Tuple

import matplotlib.pyplot as plt
import astropy.units as u
import starlink.Grf as Grf
import starlink.Ast as Ast
#import starlink.Atl as Atl

import numpy as np
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from .cornish import CornishPlot
from ..region.region import ASTRegion
from ..region.circle import ASTCircle

from cornish import ASTFITSChannel, ASTFrameSet

markerstr2value = {
	"circle" : 1,
	"cross" : 2,
	"star" : 3,
	"circle" : 4,
	"x" : 5,
	"dot" : 6,
	"triangle" : 7,
	"triangle down" : 8,
	"triangle left" : 9,
	"triangle right" : 10
}

class SkyPlot(CornishPlot):
	'''
	A convenience class providing a high level interface for creating sky plots in Matplotlib.
	
	:param extent: an ASTRegion that encompasses the full area to plot
	:param figsize: width,height of the plot figure in inches (parameter passed directly to :class:`matplotlib.figure.Figure`) 
	'''
	def __init__(self, extent:ASTRegion=None, figsize:Tuple[float,float]=(12.0, 12.0)):
		
		self.astPlot = None # type: Ast.Plot
		self._figure = None # the matplotlib.figure.Figure object
		
		# -------------------------------------------------------
		# Create frame set that will map the position in the plot
		# (i.e. pixel coordinates) to the sky (i.e. WCS)
		#fits_chan = ASTFITSChannel()
		
		naxis1 = 100
		naxis2 = 100

		if isinstance(extent, ASTCircle):
			circle_extent = extent
		elif isinstance(extent, ASTRegion):
			circle_extent = extent.boundingCircle()
			#center = circle_extent.center
		else:
			raise ValueError("Could not determine a center point for the provided extent. Hint: provide an ASTRegion object.")
				

		# The NAXIS values are arbitrary; this object is used for mapping.
		cards = {
			"CRVAL1":circle_extent.center[0], # reference point value (image center) in sky coords (RA)
			"CRVAL2":circle_extent.center[1], # reference point value (image center) in sky coords (dec)
			"CTYPE1":"RA---TAN", #"GLON-TAN", # projection type: first coordinate RA, projection tangential
			"CTYPE2":"DEC--TAN", #"GLAT-TAN"
			"CRPIX1":naxis1/2 + 0.5, # reference point (image center) point in pixel coords
			"CRPIX2":naxis2/2 + 0.5,
			"CDELT1":2.1*circle_extent.radius.to_value(u.deg)/naxis1,
			"CDELT2":2.1*circle_extent.radius.to_value(u.deg)/naxis2,
			"NAXIS1":naxis1,
			"NAXIS2":naxis2,
			"NAXES":2,
		}
		#naxis1 = cards['NAXIS1']
		#naxis2 = cards['NAXIS2']
		pix2sky_mapping = ASTFrameSet.fromFITSHeader(fits_header=cards)
		# -------------------------------------------------------
		
		#  Create a matplotlib figure, 12x12 inches in size.
		#dx = figsize[0] # 12.0
		#dy = figsize[1] # 12.0
		dx, dy = figsize
		self._figure = plt.figure( figsize=(dx,dy) ) # -> matplotlib.figure.Figure Ref: https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.pyplot.figure.html
		fig_aspect_ratio = dy/dx
		
		#  Set up the bounding box of the image in pixel coordinates, and get
		#  the aspect ratio of the image.
		bbox = (0.5, 0.5, naxis1 + 0.5, naxis2 + 0.5)
		fits_aspect_ratio = ( bbox[3] - bbox[1] )/( bbox[2] - bbox[0] )
		
		#  Set up the bounding box of the image as fractional offsets within the
		#  figure. The hx and hy variables hold the horizontal and vertical half
		#  widths of the image, as fractions of the width and height of the figure.
		#  Shrink the image area by a factor of 0.7 to leave room for annotated axes.
		if fig_aspect_ratio > fits_aspect_ratio :
		  hx = 0.5
		  hy = 0.5*fits_aspect_ratio/fig_aspect_ratio
		else:
		  hx = 0.5*fig_aspect_ratio/fits_aspect_ratio
		  hy = 0.5
		
		hx *= 0.7
		hy *= 0.7
		gbox = ( 0.5 - hx, 0.5 - hy, 0.5 + hx, 0.5 + hy )
		
		#  Add an Axes structure to the figure and display the image within it,
		#  scaled between data values zero and 100. Suppress the axes as we will
		#  be using AST to create axes.
		ax_image = self._figure.add_axes( [ gbox[0], gbox[1], gbox[2] - gbox[0],
								  gbox[3] - gbox[1] ], zorder=1 ) # -> matplotlib.axes.Axes
		ax_image.xaxis.set_visible( False )
		ax_image.yaxis.set_visible( False )
		#ax_image.imshow( hdu_list[0].data, vmin=0, vmax=200, cmap=plt.cm.gist_heat,
		#				origin='lower', aspect='auto')
		
		#  Add another Axes structure to the figure to hold the annotated axes
		#  produced by AST. It is displayed on top of the previous Axes
		#  structure. Make it transparent so that the image will show through.
		ax_plot = self._figure.add_axes( [ 0, 0, 1, 1 ], zorder=2 ) # rect = [x0, y0, width, height], 1 = full canvas size
		ax_plot.xaxis.set_visible(False)
		ax_plot.yaxis.set_visible(False)
		ax_plot.patch.set_alpha(0.0)
		
		#  Create a drawing object that knows how to draw primitives (lines,
		#  marks and strings) into this second Axes structure.
		grf = Grf.grf_matplotlib( ax_plot )
		
		#print(f"gbox: {gbox}")
		#print(f"bbox: {bbox}")
		
		# box in graphics coordinates (area to draw on, dim of plot)
		#plot = Ast.Plot( frameset.astObject, gbox, bbox, grf )
		self.astPlot = Ast.Plot( pix2sky_mapping.astObject, gbox, bbox, grf, options="Uni1=ddd:mm:ss" )
		 #, options="Grid=1" )
		#plot.set( "Colour(border)=2, Font(textlab)=3" );
		
		self.astPlot.Grid = True # can change the line properties
		self.astPlot.Format_1 = "dms"
		
		# colors:
		# 1 = black
		# 2 = red
		# 3 = lime
		# 4 = blue
		# 5 = 
		# 6 = pink
		
		self.astPlot.grid()
		self.astPlot.Width_Border = 2
		
		self.imageAxes = ax_image
		
	def figure(self):
		'''
		Return the :class:`matplotlib.figure.Figure` object for plot customization outside of this API.
		'''
		return self._figure
	
	def addRegionOutline(self, region:Union[ASTRegion,Ast.Region], colour:str="#4a7f7b", color=None, style:int=1):
		'''
		Overlay the outline of the provided region to the plot.
		
		:param region: the region to draw
		:param colour: a color name (e.g. ``black``) or hex code (e.g. ``#4a7f7b``)
		:param color: synonym for 'colour'
		:param style: line style: 1=solid, 2=solid, 3=dashes, 4=short dashes, 5=long dashes
		'''
		if isinstance(region, ASTRegion):
			region = region.astObject
		elif isinstance(region, Ast.Object):
			pass
		else:
			raise ValueError(f"Region provided must either be an ASTRegion or starlink.Ast.Object; was given '{type(region)}'.")
		
		if color:
			colour = color
		
		original_colour = self.astPlot.Colour_Border
		original_style = self.astPlot.Style

		self.astPlot.Style = style
		self.astPlot.Colour_Border = colour

		self.astPlot.regionoutline(region)
		
		self.astPlot.Colour_Border = original_colour
		self.astPlot.Style = original_style
	
	def addPoints(self, points, style:int=1, size:float=None, colour:str=None, color:str=None):
		'''
		Draw a point onto an existing plot.
		
		.. list-table:: Marker Styles
           :widths 20 25 25
		   :header-rows: 1
		   
		   * - marker_style value
		     - style
			 - string equivalent
		   * - 1
		     - small circle
			 - circle
           * - 2
            - cross
			- cross
           * - 3
            - star
			- star
           * - 4
            - larger circle
			- circle
           * - 5
            - x
			- x
           * - 6
            - pixel dot
			- dot
           * - 7
            - triangle pointing up
			- triangle
           * - 8
            - triangle pointing down
			- triangle down
           * - 9
            - triangle pointing left
			- triangle left
           * - 10
            - triangle pointing right
			- triangle right
           * - 11
            - ...
		
		:param points: point should be in degrees (e.g. list or numpy.ndarray, or a pair (list/tuple) of astropy.units.Quantity values, or a SkyCoord, or a container of these (all in the same form)
		:param style: an integer corresponding to one of the built-in marker styles
		:param colour: marker plot colour
		:param color: synonym for 'colour'
		:param size: scale point size by this value
		'''
		if isinstance(style, str):
			try:
				marker_style = markerstr2value[style]
			except KeyError:
				raise ValueError(f"The provided marker style value '{style}' is unknown.")
		else:
			# check for integer?
			marker_style = style
			
		if len(points) == 0:
			raise Exception("No points were provided to plot.")
		
		if size:
			current_marker_size = self.astPlot.Size_Markers
			self.astPlot.Size_Markers = size
		
		if color and not colour:
			colour = color
		if colour:
			# save current colour
			current_marker_colour = self.astPlot.Colour_Markers
			self.astPlot.Colour_Markers = colour

		# use first point to determine type
		
		if isinstance(points[0], SkyCoord):
			for p in points:
				point = [p.ra.to(u.rad).value, p.dec.to(u.rad).value]
				self.astPlot.mark(point, marker_style)
		elif isinstance(points[0], (np.ndarray, list, tuple)):
			for p in points:
				point = np.deg2rad(p)
				self.astPlot.mark(point, marker_style)
		elif len(points[0]) == 2 and isinstance(point[0], Quantity):
			for p in points:
				point = [p[0].to(u.rad).value, p[1].to(u.rad).value]
				self.astPlot.mark(point, marker_style)
		else:
			raise ValueError("Unhandled point type.")
				
		#self.astPlot.mark(point, marker_style)
		
		# restore plot values
		# -------------------
		if size:
			self.astPlot.Size_Markers = current_marker_size
		
		if colour:
			self.astPlot.Colour_Markers = current_marker_colour
	
	# def addCurve(self, start, finish):
	# 	'''
	# 	
	# 	'''
	# 	self.astPlot.
	
	def show(self):
		'''
		Display the plot (passthrough for :meth:`matplotlib.pyplot.show`).
		'''
		plt.show()
