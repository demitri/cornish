__version__ = '0.0.0'

import logging

from .ast_object import ASTObject

# channels
from .channel.ast_channel import ASTChannel
from .channel.fits_channel import ASTFITSChannel

# mapping
from .mapping.frame.frame import ASTFrame
from .mapping.frame.frame_set import ASTFrameSet
from .mapping.frame.sky_frame import ASTSkyFrame, ASTICRSFrame
from .mapping.frame.time_frame import ASTTimeFrame

# regions
from .region.circle import ASTCircle
from .region.box import ASTBox
from .region.polygon import ASTPolygon
from .region.region import ASTRegion


# Set up logger
# =============
try:
	cornish_logger
except NameError:
	cornish_logger = logging.getLogger('cornish')
	cornish_logger.setLevel(logging.CRITICAL)       # set log level for logger
	
	# define console logger
	console_handler = logging.StreamHandler()
	console_handler.setLevel(logging.CRITICAL)		# set log level for THIS handler
	cornish_logger.addHandler(console_handler)

