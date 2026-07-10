
'''
Deprecated constants layer.

The bare string/int constants that used to live here are deprecated in favor
of the enums in :mod:`cornish.enums` (:class:`SkySystem`, :class:`SpecSystem`,
:class:`TimeSystem`, :class:`MeshType`, :class:`PixelComparison`, …).
Accessing a deprecated name emits a :class:`DeprecationWarning` and returns
the corresponding enum member's value (the canonical AST string), so existing
code keeps working. Two historical values were typos AST rejects outright and
now map to the corrected canonical values: ``SYSTEM_WAVENUMBER``
("WAVENUMBER" → "WAVN") and ``SYSTEM_BESSELIAN`` ("BEOPCH" → "BEPOCH").
'''

import warnings

from .enums import SkySystem, SpecSystem, TimeSystem, MeshType, PixelComparison

# __all__ lists the deprecated names too so `from cornish.constants import *`
# keeps working (each name resolves through __getattr__ and warns).
__all__ = ['CENTER_CORNER', 'CORNER_CORNER', 'EQUINOX_J2000', 'EQUINOX_J2010', 'EQUINOX_2001']

# Not deprecated: AST Box "form" argument values (no enum counterpart yet).
CENTER_CORNER = 0
CORNER_CORNER = 1

# Not deprecated: equinox strings (free-form AST attribute values, not a fixed set).
EQUINOX_J2000 = "J2000"
EQUINOX_J2010 = "J2010"
EQUINOX_2001 = "2001"

_DEPRECATED_CONSTANTS = {
	# SkyFrame systems -> cornish.enums.SkySystem
	'SYSTEM_ECLIPTIC': SkySystem.ECLIPTIC.value,
	'SYSTEM_EQUATORIAL': SkySystem.FK5.value,       # "EQUATORIAL" is an AST input alias for FK5
	'SYSTEM_FK5': SkySystem.FK5.value,
	'SYSTEM_GALACTIC': SkySystem.GALACTIC.value,
	'SYSTEM_HORIZON': SkySystem.AZEL.value,
	'SYSTEM_FK4_BARYOCENTRIC': SkySystem.FK4_NO_E.value,
	'SYSTEM_GEOCENTRIC_APPARENT_EQUATORIAL': SkySystem.GAPPT.value,
	'SYSTEM_HELIOECLIPTIC': SkySystem.HELIOECLIPTIC.value,
	'SYSTEM_ICRS': SkySystem.ICRS.value,
	'SYSTEM_J2000': SkySystem.J2000.value,
	'SYSTEM_SUPERGALACTIC': SkySystem.SUPERGALACTIC.value,
	'SYSTEM_UNKNOWN': SkySystem.UNKNOWN.value,
	# SpecFrame systems -> cornish.enums.SpecSystem
	'SYSTEM_FREQUENCY_GHZ': SpecSystem.FREQUENCY.value,
	'SYSTEM_ENERGY_J': SpecSystem.ENERGY.value,
	'SYSTEM_WAVENUMBER': SpecSystem.WAVENUMBER.value,     # old value "WAVENUMBER" was rejected by AST
	'SYSTEM_VACUUM_WAVELENGTH': SpecSystem.WAVELENGTH.value,
	'SYSTEM_AIR_WAVELENGTH': SpecSystem.AIR_WAVELENGTH.value,
	'SYSTEM_RADIO_VELOCITY': SpecSystem.RADIO_VELOCITY.value,
	'SYSTEM_OPTICAL_VELOCITY': SpecSystem.OPTICAL_VELOCITY.value,
	'SYSTEM_REDSHIFT': SpecSystem.REDSHIFT.value,
	'SYSTEM_BETA_FACTOR': SpecSystem.BETA_FACTOR.value,
	'SYSTEM_APPARENT_RADIAL': SpecSystem.APPARENT_RADIAL_VELOCITY.value,
	# TimeFrame systems -> cornish.enums.TimeSystem
	'SYSTEM_MJD': TimeSystem.MJD.value,
	'SYSTEM_JULIAN_DATE': TimeSystem.JULIAN_DATE.value,
	'SYSTEM_JULIAN_EPOCH': TimeSystem.JULIAN_EPOCH.value,
	'SYSTEM_BESSELIAN': TimeSystem.BESSELIAN_EPOCH.value, # old value "BEOPCH" was a typo AST rejects
	# mesh types -> cornish.enums.MeshType
	'AST_SURFACE_MESH': MeshType.SURFACE.value,
	'AST_BOUNDARY_MESH': MeshType.BOUNDARY.value,
	# astoutline / addpixelmask comparisons -> cornish.enums.PixelComparison
	'AST_OUTLINE_LT': PixelComparison.LT.value,
	'AST_OUTLINE_LE': PixelComparison.LE.value,
	'AST_OUTLINE_EQ': PixelComparison.EQ.value,
	'AST_OUTLINE_GE': PixelComparison.GE.value,
	'AST_OUTLINE_GT': PixelComparison.GT.value,
	'AST_OUTLINE_NE': PixelComparison.NE.value,
}

__all__ += list(_DEPRECATED_CONSTANTS)

def __getattr__(name):
	if name in _DEPRECATED_CONSTANTS:
		warnings.warn(f"cornish.constants.{name} is deprecated; use the corresponding enum in cornish.enums.",
		              DeprecationWarning, stacklevel=2)
		return _DEPRECATED_CONSTANTS[name]
	raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
