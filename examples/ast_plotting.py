#!/usr/bin/env python

import sys
import warnings

import numpy as np
from cornish import ASTFrameSet, ASTBox, ASTFITSChannel

from astropy.io import fits
import matplotlib.pyplot as plt
import starlink.Grf as Grf
import starlink.Ast as Ast
import starlink.Atl as Atl

# suppress Astropy FITS header warnings
warnings.filterwarnings(action="ignore", message="The following header keyword is invalid or follows an unrecognized non-standard convention")

# create region from FITS image
hdu_list = fits.open("frame-r-000094-5-0131.fits")
box_region = ASTBox(fits_header=hdu_list[0].header)

print(f"box_region is of type '{box_region.astObject.Class}'");
print(box_region.astObject.isabox()) # --> false, why?
print(isinstance(box_region.astObject, Ast.Box))
print(f"{type(box_region)}, {type(box_region.astObject)}")

# check that box_region is a sky system before setting "system" to a sky system
#
# box_region.getregionframe() # <- check is sky frame

print(box_region.system)

# set the region as being in ICRS
# - internally, adds a new frame with the new system (and mapping)
box_region.system = "ICRS" # "Galactic" -> will produce an error if not a sky system

# Could instead use the FITS file to determine the
# frame to draw in. This is done by getting the "system"
# from box region.
#
# Checks:
#
# * make sure the system is a sky system


bounding_circle = box_region.boundingCircle()

fits_chan = ASTFITSChannel()

cards = {
	"CRVAL1":bounding_circle.center[0],  # reference point (image center) in sky coords
	"CRVAL2":bounding_circle.center[1],
	"CTYPE1":"RA---TAN", #"GLON-TAN"
	"CTYPE2":"DEC--TAN", #"GLAT-TAN"
	"CRPIX1":50.5, # reference point (image center) point in pixel coords
	"CRPIX2":50.5,
	"CDELT1":2.1*bounding_circle.radius/100,
	"CDELT2":2.1*bounding_circle.radius/100,
	"NAXIS1":100,
	"NAXIS2":100,
	"NAXES":2,
}

frameset = ASTFrameSet.fromFITSHeader(fits_header=cards)

# artificial pixel coordinates are square

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
ax_image.imshow( np.log(hdu_list[0].data), vmin=0, vmax=200, cmap=plt.cm.gist_heat,
				origin='lower', aspect='auto')

#  Add another Axes structure to the figure to hold the annotated axes
#  produced by AST. It is displayed on top of the previous Axes
#  structure. Make it transparent so that the image will show through.
ax_plot = fig.add_axes( [ 0, 0, 1, 1 ], zorder=2 ) # rect = [x0, y0, width, height], 1 = full canvas size
ax_plot.xaxis.set_visible(False)
ax_plot.yaxis.set_visible(False)
ax_plot.patch.set_alpha(0.0)

#  Create a drawing object that knows how to draw primitives (lines,
#  marks and strings) into this second Axes structure.
grf = Grf.grf_matplotlib( ax_plot )

# box in graphics coordinates (area to draw on, dim of plot)
plot = Ast.Plot( frameset.astObject, gbox, bbox, grf )
plot.Grid = True # can change the line properties

plot.grid()
plot.regionoutline(box_region.astObject)
#plot.regionoutline(bounding_circle.astObject)
plt.show()

