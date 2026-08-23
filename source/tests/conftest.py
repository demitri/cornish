
'''
Shared fixtures for the SPEC-09 plotting tests.

The golden TAN field (CRVAL 30/45, 2° radius) is the PoC field the design
session used throughout SPEC-09A; several test modules and the golden-stream
fixtures are all built from it, so it is defined once here.
'''

import pytest

#: The golden field: 2° radius at (30, 45) ICRS, in the §2.2 synthetic-header
#: form (note CDELT1 < 0 — RA increases to the left).
GOLDEN_CENTRE = (30.0, 45.0)
GOLDEN_RADIUS_DEG = 2.0
GOLDEN_NAXIS = 100
GOLDEN_CDELT = 2.1 * GOLDEN_RADIUS_DEG / GOLDEN_NAXIS

GOLDEN_TAN_CARDS = {
	"CRVAL1": GOLDEN_CENTRE[0], "CRVAL2": GOLDEN_CENTRE[1],
	"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
	"CRPIX1": GOLDEN_NAXIS / 2 + 0.5, "CRPIX2": GOLDEN_NAXIS / 2 + 0.5,
	"CDELT1": -GOLDEN_CDELT, "CDELT2": GOLDEN_CDELT,
	"NAXIS1": GOLDEN_NAXIS, "NAXIS2": GOLDEN_NAXIS, "NAXES": 2,
	"RADESYS": "ICRS",
}

#: The construction options SkyPlot v2 passes to Ast.Plot (§2.2).
GOLDEN_PLOT_OPTIONS = ("Grid=1,Border=1,Colour(Grid)=11,Width(Grid)=0.5,"
                       "Colour(Border)=11,Width(Border)=2,TickAll=0,"
                       "MinTickLen(1)=0,MinTickLen(2)=0")


def golden_frameset():
	''' The pix->sky frameset of the golden TAN field. '''
	from cornish import ASTFrameSet
	return ASTFrameSet.fromFITSHeader(fits_header=dict(GOLDEN_TAN_CARDS))


def golden_ast_plot(grf, options: str = GOLDEN_PLOT_OPTIONS):
	'''
	Build an ``Ast.Plot`` on the golden field drawing into *grf*.

	All arguments are positional: pyast constructors silently swallow keyword
	arguments (the options string is the 5th positional).
	'''
	import starlink.Ast as Ast
	frameset = golden_frameset()
	gbox = (0.0, 0.0, 1.0, 1.0)
	bbox = (0.5, 0.5, GOLDEN_NAXIS + 0.5, GOLDEN_NAXIS + 0.5)
	return Ast.Plot(frameset.astObject, gbox, bbox, grf, options)


@pytest.fixture
def golden_cards():
	return dict(GOLDEN_TAN_CARDS)
