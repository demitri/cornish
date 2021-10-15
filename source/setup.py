import re
import setuptools
from distutils.core import setup

import numpy as np

def get_property(prop:str, project:str):
	'''
	Read the requested property (e.g. '__version__', '__author__') from the specified Python module.
	Ref: https://stackoverflow.com/a/41110107/2712652
	'''
	result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
	return result.group(1)

try:
	with open('HISTORY.rst') as history_file:
		history = history_file.read()
except:
	history = ""

sources = [] # e.g. C sources
data_files = []
include_dirs = ['cornish', np.get_include()]		# -I directories
library_dirs = []			# -L directories
libraries = []		# libraries to include
# equivalent of a bare #define in C source
define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_18_API_VERSION')] # warn if using a deprecated NumPy API, defined in numpy/numpyconfig.h
extra_link_args = [] # e.g. ['-framework', 'OpenGL', '-framework', 'GLUT'])

description = ("A pure Python interface to the excellent Starlink AST library.")
#long_description = '''<long description here>'''
with open('README.rst') as readme_file:
    readme = readme_file.read()


# list of classifiers: https://pypi.org/classifiers/
classifiers = [
    "Development Status :: 4 - Beta",
    "License :: Other/Proprietary License (TBD)",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research"
]

exec(open('cornish/version.py').read())
setup(
    name="cornish",
    version=__version__,
    author="Demitri Muna",
    author_email="demitri@trillianverse.org",
    #version=get_property('__version__', 'cornish'),
    description=description,
    long_description=f"{readme}\n\n{history}",
    #long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
    #license="GPL",
    #classifiers=classifiers,
    url="https://github.com/demitri/cornish",
    #project_urls={
    #        "Bug Tracker": "https://bugs.example.com/HelloWorld/",
    #        "Documentation": "https://docs.example.com/HelloWorld/",
    #        "Source Code": "https://code.example.com/HelloWorld/",
    #    },
    project_urls={
    	"Documentation": "https://github.com/demitri/cornish",
    	"Source Code":"https://github.com/demitri/cornish",
    },
    #setup_requires=['wheel'], # needed to package for distribution
    install_requires=[
        'starlink-pyast>=3.15.3'
    ],
    packages=['cornish'],
    #packages=setuptools.find_packages(include=["cornish", "cornish.*"]),#exclude=["tests"]), # https://setuptools.readthedocs.io/en/latest/setuptools.html#id15
    data_files=data_files,
    #ext_package="cornish", # will compile the methods from the extension to the namespace "cornish"
    #ext_modules=[c_extension], # alternative: cythonize(etc), needs "from Cython.Build import cythonize"
    #include_dirs=[],
    #py_modules=['cornish'],
    python_requires='>=3.6'
)
