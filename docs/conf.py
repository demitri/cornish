
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from datetime import datetime
# sys.path.insert(0, os.path.abspath('.'))

# import from local files rather than an installed version:
sys.path.insert(0, os.path.abspath('../source'))

import cornish

# -- Project information -----------------------------------------------------

project = 'Cornish, A Python Interface to the Starlink AST WCS Library'
author = 'Demitri Muna'
copyright = f'2015-{datetime.today().year}, {author}'

# The full version, including alpha/beta/rc tags
release = cornish.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx.ext.autodoc',
	'sphinx.ext.intersphinx',
	'sphinx.ext.graphviz',
	'sphinx.ext.inheritance_diagram',
	'sphinx.ext.viewcode',
	'sphinx.ext.todo',
	'sphinx.ext.napoleon', # must appear before 'sphinx-autodoc-typehints'
	'sphinx_autodoc_typehints'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# See: https://pypi.org/project/sphinx-autodoc-typehints/
# sphinx_autodoc_typehints_options = {
# 	'set_type_checking_flag' : True
# }

# -- Extensions settings -------------------------------------------------

# number of days to cache remotely downloaded 'inv' files (default = 5)
intersphinx_cache_limit = 5

intersphinx_mapping = {
	'astropy'    : ('https://docs.astropy.org/en/stable', None),
	'numpy'      : ('https://numpy.org/doc/stable/', None),
	'matplotlib' : ('http://matplotlib.sourceforge.net', None)
}

todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_theme_options = {
	'logo_only': True,
	'display_version': True,
	'prev_next_buttons_location': 'bottom',
	'style_external_links': True,
	# Toc options
	'collapse_navigation': False,
	'sticky_navigation': False,
	'navigation_depth': 3,
	'includehidden': True,
	'titles_only': False
}
