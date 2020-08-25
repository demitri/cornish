
from typing import Union
from .cornish import CornishPlot

import matplotlib.pyplot as plt
import astropy.units as u
import starlink.Grf as Grf
import starlink.Ast as Ast
import starlink.Atl as Atl

from ..region.region import ASTRegion
from ..region.circle import ASTCircle

from cornish import ASTFITSChannel, ASTFrameSet

class SkyPlot(CornishPlot):
	'''
	A convenience class to create a sky plot in Matplotlib.
	
	:param extent: an ASTRegion that encompasses the full area to plot
	'''
	def __init__(self, extent:ASTRegion=None, figsize=(12.0, 12.0)):
		
		self.astPlot = None # type: Ast.Plot
		self._figure = None # the matplotlib.figure.Figure object
		
		# -------------------------------------------------------
		# Create frame set that will map the position in the plot
		# (i.e. pixel coordinates) to the sky (i.e. WCS)
		fits_chan = ASTFITSChannel()
		
		naxis1 = 100
		naxis2 = 100

		if isinstance(extent, ASTCircle):
			circle_extent = extent
		elif isinstance(extent, ASTRegion):
			circle_extent = extent.boundingCircle()
			center = circle_extent.center
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
		return self._figure
	
	def addRegionOutline(self, region:Union[ASTRegion,Ast.Object], color:str="#4a7f7b", style:int=1):
		'''
		Overlay the outline of the provided region to the plot.
		
		:param region: the region to draw
		:param color: a color name (e.g. ``black"``) or hex code (e.g. ``#4a7f7b``)
		:param style: line style: 1=solid, 2=solid, 3=dashes, 4=short dashes, 5=long dashes
		'''
		if isinstance(region, ASTRegion):
			region = region.astObject
		elif isinstance(region, starlink.Ast.Object):
			pass
		else:
			raise ValueError(f"Region provided must either be an ASTRegion or starlink.Ast.Object; was given '{type(region)}'.")
		
		old_color = self.astPlot.Colour_Border
		old_style = style

		self.astPlot.Style = style
		self.astPlot.Colour_Border = color

		self.astPlot.regionoutline(region)
		
		self.astPlot.Style = old_style
		self.astPlot.Colour_Border = old_color
		
	def show(self):
		'''
		'''
		self.astPlot.Colour_Border = "black"
		plt.show()
