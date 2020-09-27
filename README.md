[![Documentation Status](https://readthedocs.org/projects/cornish/badge/?version=latest)](https://cornish.readthedocs.io/en/latest/?badge=latest)

# Cornish

Cornish is a Python interface over the [Starlink AST](http://starlink.eao.hawaii.edu/starlink/AST) astronomical software library, part of the [Starlink Software Collection](http://starlink.eao.hawaii.edu/). Cornish is designed and written by Demitri Muna.

The Starlink AST library is a collection of tools for working with world coordinate systems and excels at coordinate system transformations, representing and working with regions on the celestial sphere (e.g. polygons, circles, boxes, etc.), plotting, and more.

A [Python interface](http://starlink.github.io/starlink-pyast/pyast.html) called [starlink-ast](https://pypi.org/project/starlink-pyast/) written by the library's authors, David Berry and Tim Jenness, is available, however it thinly exposes the C interface. A working knowledge of the C API is realistically a prerequisite to using `starlink-ast`. Similarly, the primary [documentation for the library](http://www.starlink.ac.uk/cgi-bin/htxserver/sun211.htx/sun211.html?xref_) is the that of the C version which does not refer to the Python interface.

The aim of Cornish is to be a fully Pythonic interface to the library. It doesn't replace the existing `starlink-ast` Python interface; rather, it is a wrapper around that. The current development focus is defining and working with regions on the sky. It accepts and returns [Astropy](https://www.astropy.org) objects where possible.

Cornish is currently under active development and all APIs are subject to change. It is not recommended to be used in a pipeline yet, but it is becoming increasingly mature and particularly useful for interactive use. The project is being released in this state as it is a dependency of the Trillian and SciID projects.

Many thanks to David Berry for the generous and extremely responsive help and advice in the development of the library.

## Installation

#### Install from PyPy

    pip install cornish

#### Install from source

    git clone https://github.com/demitri/cornish.git
    cd cornish/source
    python setup.py install

## Documentation
  
Code examples and API documentation can be found here: [https://cornish.readthedocs.io/en/latest/](https://cornish.readthedocs.io/en/latest/)

## Examples

A few examples demonstrate the library's capabilities.


```python
from cornish import ASTRegion
from astropy.io import fits

# read FITS HDU (assuming a 2D image with a WCS)
hdu1 = fits.open("my_file.fits")[0]

# define a polygon around the border the image
polygon = ASTRegion.fromFITSHeader(hdu1.header)

# get a circle that bounds that polygon
circle = polygon.boundingCircle()

# create a polygon from a circle region
polygon2 = circle.toPolygon(npoints=20)

# check if polygons overlap
polygon.overlaps(polygon3)

```

(If it seems too simple... that's kind of the point.)
