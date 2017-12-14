__version__ = '0.0.0'

import logging

from .ast_object import ASTObject
from .mapping.frame.frame import ASTFrame

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

