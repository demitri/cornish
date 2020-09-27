Working With Regions
====================

.. py:module:: cornish

Starlink AST is primarily a library for working with world coordinate systems (WCS) for astronomical data. It has extensive functionality related to coordiate system translations, platting, and more. One of the primary tools of the library (and the initial focus of the Cornish Python interface) is the handling of regions. Regions can be Cartesian, e.g. a pixel grid on a CCD, or else in a celestial sphere frame, e.g. a region on an :py:class:`~cornish.mapping.frame.ASTSkyFrame`. For the latter, all lines connecting vertices are great circles.

All region objects are all subclassed from :class:`~ASTRegion`. See the `C library reference documentation <http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211.html>`_ for the full details of the objects. Region classes include:

* :class:`~cornish.ASTCircle` [`C AST Reference <http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211ss27.html>`_]
* :class:`~cornish.ASTPolygon` [`C AST Reference <http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211ss166.html>`_]
* :class:`~cornish.ASTBox` [`C AST Reference <http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211ss22.html>`_]
* :class:`~cornish.ASTCompoundRegion` [`C AST Reference <http://starlink.eao.hawaii.edu/devdocs/sun211.htx/sun211ss35.html>`_]

When creating a region, note that the frame the region to be defined in must be specified. The class :class:`~frame.ASTICRSFrame` is provided as a convenience, defined as sky frame in the ICRS system with epoch 2000.0 and equinox 2000.0. Note that regions default to the :class:`~frame.ASTICRSFrame` if one is not provided.

Region objects are defined in the top level of the ``cornish`` namespace.

All Regions
-----------

The :class:`~cornish.ASTRegion` class is the superclass for all Cornish region objects. It contains a great deal of functionality. The API reference is a worth looking at to see what is available. A short listing:

* :attr:`~ASTRegion.points`: Return a list of points thatdescribe the region. The shape is dependent on the typr of region.
* :meth:`~ASTRegion.boundingCircle`: Returns a new :class:`~cornish.ASTCircle` object that bounds the given region. 
* :meth:`~ASTRegion.overlaps`: Check whether two regions overlap. 
* :meth:`~ASTRegion.isIdenticalTo`: Check whether two regions are identical.
* :meth:`~ASTRegion.isFullyWithin`: Check whether one region is fully within another.
* :meth:`~ASTRegion.fullyEncloses`: Check whether one region is fully encloses another.
* :meth:`~ASTRegion.boundaryPointMesh`: Returns an array of evenly distributed points that cover the boundary of the region.
* :meth:`~ASTRegion.interiorPointMesh`: Returns an array of evenly distributed points that cover the surface of the region.
* :meth:`~ASTRegion.containsPoint`: Determine whether a given point lies within a given region.

See the API reference for more methods and properties.

Circles
-------

Creating Circles
^^^^^^^^^^^^^^^^

Circles can be defined as either a center point and a radius or else a center point and another on the circumference. Coordinates can be specified an :class:`astropy.coodinates.SkyCoord` object or pairs of values in degrees. The examples below demonstrate various ways to create circle regions.

.. code-block:: python

	from cornish import ASTCircle, ASTICRSFrame, ASTSkyFrame
	from cornish.constants import SYSTEM_GALACTIC, EQUINOX_J2010
	from astropy.coordinates import SkyCoord
	import astropy.units as u
	
	# note that the default frame is ICRS, epoch=2000.0, equinox=2000.0
	
	# defined as center + radius
	# --------------------------	
	
	# using Astropy objects
	center = SkyCoord(ra="12d42m22s", dec="-32d18m58s")
	circle = ASTCircle(center=center, radius=2.0*u.deg)
	
	# using float values, defaults to degrees
	circle = ASTCircle(center=[12.7061, -31.6839], radius=2.0) # assumes degrees
	circle = ASTCircle(center=[12.7061*u.deg, -31.6839*u.deg], radius=2.0*u.deg) # Quanitites also accepted
	
	# defined as center + circumference point
	# ---------------------------------------
	circle = ASTCircle(center=center, edge_point=[12.7061, -32.6839]) # edge_point also takes SkyCoord
	
	# define the circle in another frame
	# ----------------------------------
	gal_frame = ASTSkyFrame(system=SYSTEM_GALACTIC)
	gal_frame.equinox = EQUINOX_J2010
	ASTCircle(frame=gal_frame, center=center, radius=2.0*u.deg)

Circle Properties
^^^^^^^^^^^^^^^^^

Circles have :attr:`~ASTCircle.radius` and :attr:`~ASTCircle.centre` properties as one might expect, and both can be directly modified:

.. code-block:: python

	circle.radius
	>>> <Quantity 2. deg>
	
	circle.centre # or "center" if you prefer...
	>>> array([ 12.70611111, -32.31611111]) # output in degrees
		
New circles can be created by a scale factor or increased by addition from an existing circle.

.. code-block:: python
	
	scaled_circle = circle * 2.0
	scaled_circle.radius
	>>> <Quantity 4. deg>
	
	larger_circle = circle + 6*u.deg
	larger_circle.radius
	>>> <Quantity 8. deg>
	
Converting to Polygons
^^^^^^^^^^^^^^^^^^^^^^

For code that requires a polygon region as an input the method :meth:`~cornish.ASTCircle.toPolygon` will convert a circle to an :class:`~cornish.ASTPolygon`. The default is to sample 200 points for the polygon but this can be customized by using the `npoints` parameter (often even 20 are sufficient). Note that all of the polygon points fall on the circle's circumference, so the resulting region is fully inscribed by the original circle.

.. code-block:: python

	polygon = circle.toPolygon()
	finer_polygon = circle(toPolygon(npoints=200))
	
All regions have a :py:meth:`~cornish.ASTRegion.boundingCircle` property that returns an :class:`~cornish.ASTCircle` that bounds the region. In the case of :class:`~cornish.ASTCircle` objects, this method returns the original circle.

Polygons
--------

A polygon is a collection of vertices that lie in a specific frame. The default frame :class:`~cornish.ASTICRSFrame` is used if none is specified.

.. code-block:: python

    from cornish import ASTPolygon, ASTICRSFrame
    import numpy as np
    
    points = np.array([[ 12.70611111, -30.31611111],
                       [ 13.42262189, -30.41196836],
                       [ 14.07300863, -30.69069244],
                       [ 14.59623325, -31.12642801],
                       [ 14.94134955, -31.67835614],
                       [ 15.07227821, -32.29403528],
                       [ 14.97204342, -32.91392471],
                       [ 14.6459242 , -33.47688136],
                       [ 14.12273328, -33.92626054],
                       [ 13.4533703 , -34.21603194],
                       [ 12.70611111, -34.31611111],
                       [ 11.95885193, -34.21603194],
                       [ 11.28948894, -33.92626054],
                       [ 10.76629802, -33.47688136],
                       [ 10.4401788 , -32.91392471],
                       [ 10.33994401, -32.29403528],
                       [ 10.47087267, -31.67835614],
                       [ 10.81598897, -31.12642801],
                       [ 11.3392136 , -30.69069244],
                       [ 11.98960033, -30.41196836]])
    polygon = ASTPolygon(frame=ASTICRSFrame(), points=points)

Points can be specified as an array of coordinate points (as above) or as parallel arrays of each dimension (which is just ``points.T`` from above):

.. code-block:: python

    points = np.array([[ 12.70611111,  13.42262189,  14.07300863,  14.59623325, ...],
                       [-30.31611111, -30.41196836, -30.69069244, -31.12642801, ...]])

.. todo:: Provide example of how to convert a region from one frame to another.

Boxes
-----

.. todo:: Box section coming soon! (But it's pretty straightforward from the :class:`~ASTBox` API.)

Compound Regions
----------------

.. todo:: Counpound regions section coming soon! (But it's pretty straightforward from the :class:`~ASTBox` API.)
