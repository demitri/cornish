
import pytest
import numpy as np
from numpy.testing import assert_approx_equal

from cornish import ASTPolygon, ASTICRSFrame
import astropy.units as u

# points generated from:
#center = SkyCoord(ra="12d42m22s", dec="-32d18m58s")
#circle = ASTCircle(center=center, radius=2.0*u.deg)
#polygon = circle.toPolygon(npoints=20)

# center ~ [12.70611111, -32.31611111] deg

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

# polygon_points, area (sr)
sky_polygon_area_tests = [
	("((0,0), (90,0), (90,90))", 1.5707963267948966)
]

def test_polygon_creation_coordinate_pairs():
	'''
	Test the creation of a polygon using points as coordinate pairs.
	'''
	polygon = ASTPolygon(frame=ASTICRSFrame(), points=points)

def test_polygon_creation_parallel_arrays():
	'''
	Test the creation of a polygon using parallel arrays of ra,dec.
	'''
	polygon = ASTPolygon(frame=ASTICRSFrame(), points=points.T)

def test_polygon_correct_boundedness():
	'''
	Check that a newly created polygon is not unbounded.
	'''
	polygon = ASTPolygon(frame=ASTICRSFrame(), points=points)

	assert polygon.isBounded == True, "Polygon is not bounded."

	assert polygon.containsPoint([12.70611111, -32.31611111]), "Center point of polygon not inside polygon."

@pytest.mark.parametrize("polygon_points, area", sky_polygon_area_tests)
def test_polygon_area(polygon_points, area):
	'''
	Test the area of a polygon on the sky.
	'''
	polygon = ASTPolygon(frame=ASTICRSFrame(), points=polygon_points)
	assert_approx_equal(polygon.area.to(u.sr).value, area)
