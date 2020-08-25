#!/usr/bin/env python

'''
This script provides an example on how to draw regions and image data
using the cornish.SkyPlot class.
'''

import io
import bz2
import sys
import requests
from astropy.io import fits
from cornish import ASTCircle, ASTPolygon
from cornish.plot.matplotlib import SkyPlot

import matplotlib.pyplot as plt

filename = "/Users/demitri/Documents/Repositories/GitLab/Nightlight/Nightlight/Sample FITS files/frame-r-006073-4-0063.fits.bz2"

hdu_list = fits.open(filename)
hdu = hdu_list[0]

fits_bc = ASTPolygon.fromFITSHeader(hdu.header).boundingCircle()

bc = ASTCircle(frame=fits_bc.frame,
			   center=fits_bc.center,
			   radius=1.15*fits_bc.radius)

#print(fits_bc)
#print(bc)

skyplot = SkyPlot(extent=bc, figsize=(5,5))

skyplot.addRegionOutline(region=fits_bc, color="#106942", style=1)
skyplot.addRegionOutline(region=fits_bc.toPolygon(npoints=7), color="orange", style=3)

fits_file = fits.open(name=filename)
data = fits_file[0].data

#ax_image.imshow( data, vmin=0, vmax=200, cmap=plt.cm.gist_heat,
#				origin='lower', aspect='auto')


# todoL wrap this into the SkyPlot API
# imshow ref: https://matplotlib.org/3.2.1/api/_as_gen/matplotlib.axes.Axes.imshow.html?highlight=imshow#matplotlib.axes.Axes.imshow
skyplot.imageAxes.imshow( data, vmin=-50, vmax=200, cmap=plt.cm.gist_heat,origin='lower', aspect='auto')


plt.show()


