.. Cornish documentation master file, created by
   sphinx-quickstart on Thu May  7 16:08:35 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cornish: A Python Interface for the Starlink AST Library
========================================================

Cornish is a Python interface over the Starlink AST astronomical software library, part of the Starlink Software Collection. Cornish is designed and written by Demitri Muna.

The Starlink AST library is a collection of tools for working with world coordinate systems and excels at coordinate system transformations, representing and working with regions on the celestial sphere (e.g. polygons, circles, boxes, etc.), plotting, and more.

A Python interface called starlink-ast written by the library's authors, David Berry and Tim Jenness, is available, however it thinly exposes the C interface. A working knowledge of the C API is realistically a prerequisite to using starlink-ast. Similarly, the primary documentation for the library is the that of the C version which does not refer to the Python interface.

The aim of Cornish is to be a fully Pythonic interface to the library. It doesn't replace the existing starlink-ast Python interface; rather, it is a wrapper around that. The current development focus is defining and working with regions on the sky. It accepts and returns Astropy objects where possible.

Cornish is currently under active development and all APIs are subject to change. It is not recommended to be used in a pipeline yet, but it is becoming increasingly mature and particularly useful for interactive use. The project is being released in this state as it is a dependency of the Trillian and SciID projects.

Many thanks to David Berry for the generous and extremely responsive help and advice in the development of the library.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Cornish Regions API
-------------------

.. toctree::
   :maxdepth: 2

   regions

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
