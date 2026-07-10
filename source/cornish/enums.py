
'''
Enumerations for values drawn from fixed AST value sets.

Every enum here is welded to the underlying AST source of truth: integer
enums are built directly on the ``starlink.Ast`` module constants, and string
enums hold the *canonical* values AST reports via ``get("System")`` (AST
canonicalizes input aliases; the accepted alias sets are documented per enum
and verified against AST 9.3 / starlink-pyast 4.0.1).

At AST-facing boundaries pass ``enum.value`` (or the enum itself for the
IntEnums — they compare and format as their integer values).
'''

from enum import IntEnum, StrEnum

import starlink.Ast as Ast

__all__ = ['OverlapType', 'RegionOperation', 'MeshType', 'PixelComparison', 'Marker',
           'SkySystem', 'SpecSystem', 'TimeSystem',
           'SKY_SYSTEM_INPUT_ALIASES', 'SPEC_SYSTEM_INPUT_ALIASES']

class OverlapType(IntEnum):
	'''
	The return codes of ``astOverlap`` (:meth:`starlink.Ast.Region.overlap`),
	describing how two regions overlap.
	'''
	NO_FRAME_MAPPING = 0 #: the frames of the two regions could not be mapped to each other
	NONE = 1             #: no overlap
	INSIDE = 2           #: the first region is completely inside the second
	CONTAINS = 3         #: the second region is completely inside the first
	PARTIAL = 4          #: partial overlap
	IDENTICAL = 5        #: the regions are identical to within their uncertainties
	NEGATION = 6         #: the second region is the exact negation of the first

class RegionOperation(IntEnum):
	'''
	Boolean operators for combining regions (``astCmpRegion``); built on the
	``starlink.Ast`` module constants.
	'''
	AND = Ast.AND
	OR = Ast.OR
	XOR = Ast.XOR

class MeshType(IntEnum):
	'''
	The ``surface`` argument of ``astGetRegionMesh``: whether the mesh covers
	the boundary of the region or its interior.

	Note the AST convention: ``SURFACE`` (0) is the *interior* mesh ("surface"
	in the sense of a filled area), ``BOUNDARY`` (1) traces the boundary.
	'''
	SURFACE = 0
	BOUNDARY = 1

class PixelComparison(IntEnum):
	'''
	Pixel-value comparison operators for ``Ast.outline`` and
	``Moc.addpixelmask``; built on the ``starlink.Ast`` module constants.
	(Note: ``Ast.LT`` is genuinely 11, not 1 — verified.)
	'''
	LT = Ast.LT
	LE = Ast.LE
	EQ = Ast.EQ
	GE = Ast.GE
	GT = Ast.GT
	NE = Ast.NE

class Marker(IntEnum):
	'''
	Plot marker styles for :meth:`starlink.Ast.Plot.mark` (the Grf marker
	type codes).
	'''
	SMALL_CIRCLE = 1
	CROSS = 2
	STAR = 3
	CIRCLE = 4
	X = 5
	DOT = 6
	TRIANGLE = 7
	TRIANGLE_DOWN = 8
	TRIANGLE_LEFT = 9
	TRIANGLE_RIGHT = 10

class SkySystem(StrEnum):
	'''
	Canonical AST SkyFrame ``System`` values — the strings ``get("System")``
	returns. AST additionally accepts input aliases which it canonicalizes on
	set: see :data:`SKY_SYSTEM_INPUT_ALIASES`.
	'''
	ICRS = "ICRS"
	FK4 = "FK4"
	FK4_NO_E = "FK4-NO-E"
	FK5 = "FK5"
	J2000 = "J2000"
	GALACTIC = "GALACTIC"
	SUPERGALACTIC = "SUPERGALACTIC"
	ECLIPTIC = "ECLIPTIC"
	HELIOECLIPTIC = "HELIOECLIPTIC"
	AZEL = "AZEL"
	GAPPT = "GAPPT"
	UNKNOWN = "Unknown"

#: Input aliases AST accepts for SkyFrame ``System`` and the canonical value
#: each one becomes (verified against AST 9.3).
SKY_SYSTEM_INPUT_ALIASES = {
	"FK4_NO_E": SkySystem.FK4_NO_E,
	"EQUATORIAL": SkySystem.FK5,
	"GEOCENTRIC": SkySystem.GAPPT,
	"APPARENT": SkySystem.GAPPT,
}

class SpecSystem(StrEnum):
	'''
	Canonical AST SpecFrame ``System`` values (verified against AST 9.3;
	input aliases in :data:`SPEC_SYSTEM_INPUT_ALIASES`).
	'''
	FREQUENCY = "FREQ"
	ENERGY = "ENER"
	WAVENUMBER = "WAVN"
	WAVELENGTH = "WAVE"
	AIR_WAVELENGTH = "AWAV"
	RADIO_VELOCITY = "VRAD"
	OPTICAL_VELOCITY = "VOPT"
	REDSHIFT = "ZOPT"
	BETA_FACTOR = "BETA"
	APPARENT_RADIAL_VELOCITY = "VELO"

#: Input aliases AST accepts for SpecFrame ``System`` and the canonical value
#: each one becomes (verified against AST 9.3). Note that the historical
#: cornish constant value "WAVENUMBER" is NOT accepted by AST (it raises);
#: the accepted input is "WAVN".
SPEC_SYSTEM_INPUT_ALIASES = {
	"ENERGY": SpecSystem.ENERGY,
	"WAVELEN": SpecSystem.WAVELENGTH,
	"AIRWAVE": SpecSystem.AIR_WAVELENGTH,
	"VRADIO": SpecSystem.RADIO_VELOCITY,
	"VOPTICAL": SpecSystem.OPTICAL_VELOCITY,
	"REDSHIFT": SpecSystem.REDSHIFT,
	"VREL": SpecSystem.APPARENT_RADIAL_VELOCITY,
}

class TimeSystem(StrEnum):
	'''
	Canonical AST TimeFrame ``System`` values (verified against AST 9.3).
	The historical cornish constant value "BEOPCH" was a typo AST rejects;
	the correct value is "BEPOCH".
	'''
	MJD = "MJD"
	JULIAN_DATE = "JD"
	JULIAN_EPOCH = "JEPOCH"
	BESSELIAN_EPOCH = "BEPOCH"
