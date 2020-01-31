#!/usr/bin/env python

# WORK IN PROGRESS - this script doesn't work yet.

import cornish
from cornish.mapping import ASTSkyFrame, ASTFrameSet
import starlink.Ast as Ast

j2000_frame = ASTSkyFrame() # default frame: ICRS
gaia_frame = ASTSkyFrame(system="ICRS")
gaia_frame.epoch = 2015.5

#frame_set = ASTFrameSet(base_frame=j2000_frame)
#frame_set.addToBaseFrame(gaia_frame)
# from, to, domainlist (optional)
#frame_set_converter = Ast.Frame.convert(j2000_frame.astObject, gaia_frame.astObject)
frame_set_converter = ASTFrameSet.fromFrames(gaia_frame, j2000_frame)
print(frame_set_converter.astObject)

#print(frame_set_converter.convert([12.345], [-32.44]))
print(frame_set_converter.convert([1.5], [-0.57]))

