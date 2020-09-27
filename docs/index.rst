.. Cornish documentation master file, created by
   sphinx-quickstart on Thu May  7 16:08:35 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cornish: A Python Interface for the Starlink AST WCS Library
============================================================

Cornish is a Python interface over the Starlink AST astronomical software library, part of the Starlink Software Collection. Cornish is designed and written by Demitri Muna.

The Starlink AST library is a collection of tools for working with world coordinate systems and excels at coordinate system transformations, representing and working with regions on the celestial sphere (e.g. polygons, circles, boxes, etc.), plotting, and more.

A Python interface called `starlink-pyast <http://starlink.github.io/starlink-pyast/pyast.html>`_ (``pip install starlink-pyast``) written by the library's authors, David Berry and Tim Jenness, is available, however it thinly exposes the C interface. A working knowledge of the C API is realistically a prerequisite to using ``starlink-pyast``. Similarly, the primary documentation for the library is the that of the C version which does not refer to the Python interface.

The aim of Cornish is to be a fully Pythonic interface to the library. It doesn't replace the existing ``starlink-pyast`` Python interface; rather, it is a wrapper around it. The current development focus is defining and working with regions and coordinate systems on the sky. It accepts and returns Astropy objects where possible. While care is taken to mimic the original API where it makes sense, the library places greater emphasis on making the API as Pythonic as possible and will rename methods where it will make operations and concepts more clear to the Python user.

Cornish is currently under active development and all APIs are subject to change. It is not recommended to be used in a pipeline yet, but it is becoming increasingly mature and particularly useful for interactive use. The project is being released in this state as it is a dependency of the Trillian and SciDD projects, also written by Demitri Muna.

Many thanks to David Berry for the generous and extremely responsive help and advice in the development of the library.

Examples
--------

.. toctree::
   :maxdepth: 2
   
   examples_regions
   examples_plotting
   
Cornish API
-----------

.. toctree::
   :maxdepth: 2

   api_regions
   api_mapping_and_frames
   api_plotting
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
