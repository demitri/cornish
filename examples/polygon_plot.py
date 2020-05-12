#!/usr/bin/env python

# This script plots a polygon created from points.

import sys
import warnings

import numpy as np
from cornish import ASTPolygon
from cornish import ASTICRSFrame, ASTFrameSet, ASTBox, ASTFITSChannel, ASTCircle, ASTCompoundRegion

from astropy.io import fits
import matplotlib.pyplot as plt
import starlink.Grf as Grf
import starlink.Ast as Ast
import starlink.Atl as Atl

points = np.array([[ 24.9220814, -2.32553877e-01],
 [ 24.8690619, -2.13198227e-01],
 [ 24.8080071, -1.56379062e-01],
 [ 24.7961841, -1.32038075e-01],
 [ 24.7603950, -3.85297093e-02],
 [ 24.7542463,  1.17204538e-02],
 [ 24.7542463,  1.17204538e-02],
 [ 24.7769168,  8.05598701e-02],
 [ 24.8216773,  1.66088007e-01],
 [ 24.8332202,  1.83192204e-01],
 [ 24.8724133,  2.11948177e-01],
 [ 24.9190898,  2.36432081e-01],
 [ 25.0443598,  2.38795506e-01],
 [ 25.0694520,  2.34736352e-01],
 [ 25.1083355,  2.26095876e-01],
 [ 25.1263480,  2.15632984e-01],
 [ 25.1730281,  1.80042404e-01],
 [ 25.2055145,  1.40829694e-01],
 [ 25.2459929,  1.54015514e-02],
 [ 25.2459929,  1.54015514e-02],
 [ 25.2426565, -3.86550958e-02],
 [ 25.2219908, -9.67569213e-02],
 [ 25.1986233, -1.49820068e-01],
 [ 25.0872686, -2.32297073e-01]])
icrs_frame = ASTICRSFrame()
polygon = ASTPolygon(frame=icrs_frame, points=points)

galex_circle = ASTCircle(frame=icrs_frame, center=[24.617269485878584,0.2727299618460874], radius=1.1312250143591236)

compound_region = ASTCompoundRegion(regions=[polygon, galex_circle], operation=Ast.AND)

# define the extend of the plot
#bounding_circle = polygon.boundingCircle()
bounding_circle = compound_region.boundingCircle()

# -------------------------------------------------------
# Create frame set that will map the position in the plot
# (i.e. pixel coordinates) to the sky (i.e. WCS)

fits_chan = ASTFITSChannel()

cards = {
	"CRVAL1":bounding_circle.center[0], # reference point (image center) in sky coords
	"CRVAL2":bounding_circle.center[1],
	"CTYPE1":"RA---TAN", #"GLON-TAN", # projection type
	"CTYPE2":"DEC--TAN", #"GLAT-TAN",
	"CRPIX1":50.5, # reference point (image center) point in pixel coords
	"CRPIX2":50.5,
	"CDELT1":2.1*bounding_circle.radius/100,
	"CDELT2":2.1*bounding_circle.radius/100,
	"NAXIS1":100,
	"NAXIS2":100,
	"NAXES":2,
}
naxis1 = cards['NAXIS1']
naxis2 = cards['NAXIS2']
pix2sky_mapping = ASTFrameSet.fromFITSHeader(fits_header=cards)
# -------------------------------------------------------
#pix2sky_mapping.system = "Galactic"

print(bounding_circle.center)

#  Create a matplotlib figure, 12x12 inches in size.
dx=12.0
dy=12.0
fig = plt.figure( figsize=(dx,dy) )
fig_aspect_ratio = dy/dx

#  Set up the bounding box of the image in pixel coordinates, and get
#  the aspect ratio of the image.
naxis1 = int(cards["NAXIS1"])
naxis2 = int(cards["NAXIS2"])
bbox = (0.5, 0.5, naxis1 + 0.5, naxis2 + 0.5)
fits_aspect_ratio = ( bbox[3] - bbox[1] )/( bbox[2] - bbox[0] )
#fits_aspect_ratio = 1

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
ax_image = fig.add_axes( [ gbox[0], gbox[1], gbox[2] - gbox[0],
						  gbox[3] - gbox[1] ], zorder=1 )
ax_image.xaxis.set_visible( False )
ax_image.yaxis.set_visible( False )
#ax_image.imshow( hdu_list[0].data, vmin=0, vmax=200, cmap=plt.cm.gist_heat,
#				origin='lower', aspect='auto')

#  Add another Axes structure to the figure to hold the annotated axes
#  produced by AST. It is displayed on top of the previous Axes
#  structure. Make it transparent so that the image will show through.
ax_plot = fig.add_axes( [ 0, 0, 1, 1 ], zorder=2 )
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
plot = Ast.Plot( pix2sky_mapping.astObject, gbox, bbox, grf, options="Uni1=ddd:mm:ss" )
 #, options="Grid=1" )
#plot.set( "Colour(border)=2, Font(textlab)=3" );

plot.Grid = True # can change the line properties
plot.Format_1 = "dms"

# colors:
# 1 = black
# 2 = red
# 3 = lime
# 4 = blue
# 5 = 
# 6 = pink

plot.grid()
plot.Width_Border = 2

#plot.Colour_Border = "#0099cc"
#plot.regionoutline(bounding_circle.astObject)

plot.Colour_Border = "#106942"
plot.regionoutline(polygon.astObject)

plot.Colour_Border = "blue"
plot.regionoutline(galex_circle.astObject)
#plt.plot(galex_circle.center[0],galex_circle.center[1],'ro')

plot.Colour_Border = "red"
plot.Style = 3
plot.regionoutline(bounding_circle.astObject)

plt.show()

