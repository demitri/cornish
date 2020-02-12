#!/usr/bin/env python

import astropy.units as u
from cornish.channel import ASTFITSChannel
#from cornish.mapping.frame import FrameSet
from cornish import ASTFrameSet, ASTCircle, ASTICRSFrame
from astropy.io import fits

hdu_list = fits.open("frame-r-000094-5-0131.fits")
hdu1 = hdu_list[0]

# create frameset that describes image array
# contains pixel frame and WCS frame
frameset = ASTFrameSet.fromFITSHeader(hdu1.header)

# need mapping between pixel coordinate to sky coords
pix2sky_mapping = frameset.mapping #.astObject
sky2pix_mapping = pix2sky_mapping.inverseMapping()

# usage:
# sky2pix_mapping.transform([[ra1,ra2],[dec1,dec2]])

# want to create a circle, but on the sky
# need to create a frame the circle will be in
# leaving epoch unspecified will leave it to be inherited 
# from FITS file later
icrs_frame = ASTICRSFrame()

# get a new frame set that contains a mapping
# from the original base frame (FITS file)
# to the template icrs_frame
#
# this new frameset defines a mapping from the pixel to sky frame
new_frameset = frameset.framesetWithMappingTo(icrs_frame)

# new_frameset is created to guarantee that
# axis 1 = RA, axis 2 = dec, regardless of what WCS
# the FITS file is - this is unnecessary if the FITS
# file frameset uses axis1=ra, axis2=dec,
# but this guarantees the mapping is correct.
# this hides the details of what coordinates the fits
# file are define in, and Just Works.

# create a circle in the sky frame
circle_region = ASTCircle(frame=new_frameset.currentFrame,
						  center=[354.47207716*u.deg, 0.75203964*u.deg],
						  radius=30*u.arcsec)

# get mapping from region to pixel
region2pixel_mapping = new_frameset.mapping.inverseMapping()

region.mask(map=region2pixel_mapping,
			mask_inside=True,
			fits_coordinates=True,
			image=hdu1.data) #note overwrites image
			

