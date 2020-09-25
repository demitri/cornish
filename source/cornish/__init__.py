
import logging

from .ast_object import ASTObject

# channels
from .channel.ast_channel import ASTChannel
from .channel.fits_channel import *

# mapping
from .mapping.frame.frame import *
from .mapping.frame.frame_set import *
from .mapping.frame.sky_frame import *
from .mapping.frame.time_frame import *
from .mapping.frame.compound_frame import *

# regions
from .region.box import *
from .region.circle import *
from .region.compound_region import *
from .region.polygon import *
from .region.region import *

from .version import __version__

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

