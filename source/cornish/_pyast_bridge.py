
'''
The single seam between cornish and pyast idioms (SPEC-04 / SPEC-04A).

This private module owns every unit and array-shape conversion between
user-facing coordinate forms (SkyCoord, Quantity, bare arrays in degrees,
strings) and the ``(naxes, npoints)`` float64 radians-or-native arrays that
pyast wants. After the SPEC-04 migration, the functions here are the ONLY
``deg2rad``/``rad2deg`` call sites in cornish (enforced by the
``test_no_conversions_outside_bridge`` gate).

It also owns the effective-frame dispatch (`is_sky`/`current_frame_of`) that
kills the verified pyast footguns::

    Ast.FrameSet.isaskyframe()  -> False even when the current frame is a SkyFrame
    Ast.Region.isaskyframe()    -> False even for a region built on a SkyFrame

plus uncertainty materialization (`as_uncertainty_region`), the sanctioned
dump-reading escape hatch (`dump_value`), and the pickle helper
(`reconstruct`).

This is a PRIVATE module: its API may move only in step with all internal
callers, and it is exported nowhere. It returns bare ``np.ndarray`` only —
rich returns (``asSkyCoord=True`` etc.) are the API layer's job.

The behavior contract (decision tables, error taxonomy, verified pyast
behaviors V1–V18) is frozen in SPEC-04A; do not change behavior here without
updating that document.
'''

from __future__ import annotations
from typing import Optional, Sequence, Union, Literal, TYPE_CHECKING

import re
import ast as _python_ast

import numpy as np
import astropy.units as u
import astropy.coordinates
from astropy.coordinates import SkyCoord
from astropy.time import Time
import starlink.Ast as Ast

from .ast_object import ASTObject

if TYPE_CHECKING:
    from .mapping.frame.frame import ASTFrame

#: Anything that carries a frame the bridge can resolve: a pyast Frame
#: (including SkyFrame, FrameSet, Plot, Region — all Frame subclasses) or a
#: cornish wrapper around one. Plain non-Frame mappings are NOT FrameLike.
FrameLike = Union[Ast.Frame, "ASTFrame"]

#: Accepted point-input forms; the exhaustive contract is the SPEC-04A §2 decision table.
PointsInput = Union[SkyCoord, u.Quantity, np.ndarray, Sequence, str]


# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

def _unwrap(obj) -> Ast.Object:
    ''' cornish wrapper -> its ``.astObject``; pyast object -> itself; anything else -> TypeError. '''
    if isinstance(obj, ASTObject):
        obj = obj.astObject
    if isinstance(obj, Ast.Object):
        return obj
    raise TypeError(f"Expected an AST object or a cornish wrapper around one; got '{type(obj).__name__}'.")


def _check_finite(values: np.ndarray, context: str):
    '''
    Raise ValueError if any value is NaN, ±inf, or Ast.BAD (P18/D13).

    This must run BEFORE any unit conversion: deg2rad(Ast.BAD) silently
    corrupts the sentinel (verified V11), and NaN input is always an
    upstream bug (fail-loud rule).
    '''
    values = np.asarray(values)
    offending = ~np.isfinite(values) | (values == Ast.BAD)
    if np.any(offending):
        indices = np.argwhere(offending).tolist()
        raise ValueError(
            f"{context} contains non-finite or Ast.BAD value(s) at index/indices {indices}; "
            f"coordinate values must be finite."
        )


def _coerce_array(arr: np.ndarray, naxes: int, parallel_axes: Optional[bool]) -> tuple:
    '''
    Apply the SPEC-04A §2.2 shape rules to a plain numeric array.

    :returns: (array of shape (naxes, npoints), is_single_point_form)
    '''
    if arr.size == 0:
        raise ValueError("No points provided (empty input).")  # P17
    if arr.ndim == 1:
        if parallel_axes is not None:
            raise ValueError("'parallel_axes' applies only to rank-2 array input; "
                             "a rank-1 single point has no orientation to declare.")
        if arr.shape[0] != naxes:
            raise ValueError(f"A rank-1 point must have exactly naxes={naxes} values "
                             f"(got shape {arr.shape}).")  # P19
        return arr.reshape(naxes, 1), True
    if arr.ndim != 2:
        raise ValueError(f"Accepted point-array shapes are ({naxes},), (n, {naxes}), and "
                         f"({naxes}, n); got shape {arr.shape}.")  # P19
    a, b = arr.shape
    if a == naxes and b == naxes:
        # the square case: DECISIONS D3 — interpreted as PAIRS (rows are points)
        # unless the caller explicitly declares parallel_axes=True
        if parallel_axes is True:
            return arr, False
        return np.array(arr.T), False
    if a == naxes:
        # rows are axes (parallel form)
        if parallel_axes is False:
            raise ValueError(f"parallel_axes=False contradicts the input shape {arr.shape}: "
                             f"only the first dimension matches naxes={naxes}.")
        return arr, False
    if b == naxes:
        # rows are points (pairs form)
        if parallel_axes is True:
            raise ValueError(f"parallel_axes=True contradicts the input shape {arr.shape}: "
                             f"only the second dimension matches naxes={naxes}.")
        return np.array(arr.T), False
    raise ValueError(f"Accepted point-array shapes are ({naxes},), (n, {naxes}), and "
                     f"({naxes}, n); got shape {arr.shape}.")  # P19


def _skycoord_to_rad(sc: SkyCoord, eff_frame: Ast.Frame) -> np.ndarray:
    '''
    Convert a scalar or 1-D array SkyCoord into a (2, n) radians array in the
    effective frame's system (§2.3). Always transforms — even a SkyCoord that
    appears to already be in the target frame — to remove the
    "matched by name but differed by equinox" bug class.
    '''
    if sc.ndim > 1:
        raise ValueError(f"Multi-dimensional SkyCoord arrays are not accepted "
                         f"(got shape {sc.shape}); flatten yourself.")  # P2 rule
    target = astropy_frame_for(eff_frame)
    converted = sc.transform_to(target)
    # read frame-neutrally: .ra raises AttributeError on e.g. galactic SkyCoords (V6)
    spherical = converted.spherical
    lon = np.atleast_1d(spherical.lon.rad)
    lat = np.atleast_1d(spherical.lat.rad)
    return np.stack([lon, lat])


def _quantity_to_rad_values(q: u.Quantity) -> np.ndarray:
    '''
    Angular Quantity -> plain radians ndarray, with the P15/P16/P18 checks.
    Callers have already established the frame is a sky frame (D12).
    '''
    if q.unit.physical_type == 'dimensionless':
        # Decided loud because .to(u.rad) SUCCEEDS on dimensionless input, which
        # would create a units trap inside the module that exists to kill units traps.
        raise ValueError("Ambiguous input: a dimensionless Quantity would silently read as "
                         "radians; attach an angular unit or pass a bare array (degrees).")  # P16
    raw = np.asarray(q.value, dtype=np.float64)
    _check_finite(raw, "the provided Quantity")  # P18, before conversion
    try:
        return np.asarray(q.to_value(u.rad), dtype=np.float64)
    except u.UnitConversionError as e:
        raise ValueError(f"Cannot interpret a Quantity with unit '{q.unit}' as sky coordinates; "
                         f"an angular unit is required.") from e  # P15


_NONSKY_REFUSAL = ("{what} cannot be interpreted in a non-sky frame: sky coordinates are "
                   "uninterpretable there, and a Quantity's unit cannot be reconciled with "
                   "unknowable native frame units — pass bare numbers in the frame's native "
                   "units instead (DECISIONS D12).")


def _convert(points: PointsInput, frame: FrameLike, *,
             parallel_axes: Optional[bool] = None) -> tuple:
    '''
    The shared input classifier/converter behind `to_frame_units` and
    `is_single_point` (they must raise identically — same classifier).

    :returns: (float64 array of shape (naxes, npoints), is_single_point_form)
    '''
    # frame resolution happens FIRST, independently of the points (§2.1)
    eff = current_frame_of(frame)
    sky = eff.isaskyframe()
    naxes = int(eff.get("Naxes"))

    if parallel_axes not in (None, True, False):
        raise TypeError(f"'parallel_axes' must be True, False, or None "
                        f"(got '{type(parallel_axes).__name__}').")

    # Input dispatch order is frozen (§2.2): each step matches by an explicit
    # isinstance test, NEVER by whether a coercion succeeds — coercion-success
    # dispatch is what created the P16 trap (u.Quantity([10.0, 20.0]) succeeds
    # on a bare list as dimensionless, and u.Quantity("(10,20)") raises
    # TypeError, which would kill string input).

    # step 1: str (checked first — a str is also a Python sequence)
    if isinstance(points, str):
        if parallel_axes is not None:
            raise ValueError("'parallel_axes' applies only to array input (got a points string).")
        arr = parse_points_string(points)
        _check_finite(arr, "the parsed points string")  # P18
        arr, single = _coerce_array(arr, naxes, None)
        if sky:
            arr = np.deg2rad(arr)  # string grammar carries degrees on sky frames (P14)
        return np.ascontiguousarray(arr, dtype=np.float64), single

    # P17 emptiness — a guarded probe, not a bare len() call: scalar SkyCoord,
    # scalar Quantity, and rank-0 ndarray all raise TypeError under len() and
    # are by definition non-empty single/scalar forms. Runs before steps 3/5 so
    # their vacuously-true all(...) never sees an empty sequence.
    try:
        length = len(points)
    except TypeError:
        length = None
    if length == 0:
        raise ValueError("No points provided (empty input).")  # P17

    # step 2: SkyCoord
    if isinstance(points, SkyCoord):
        if parallel_axes is not None:
            raise ValueError("'parallel_axes' applies only to array input (got a SkyCoord).")
        if not sky:
            raise ValueError(_NONSKY_REFUSAL.format(what="A SkyCoord"))  # P1/P2 basic-frame rule
        rad = _skycoord_to_rad(points, eff)
        if rad.size == 0:
            raise ValueError("No points provided (empty input).")  # P17 re-applied
        _check_finite(rad, "the converted SkyCoord coordinates")
        return np.ascontiguousarray(rad, dtype=np.float64), points.isscalar

    # step 3: non-str sequence whose elements are all SkyCoord
    if isinstance(points, (list, tuple)) and all(isinstance(el, SkyCoord) for el in points):
        if parallel_axes is not None:
            raise ValueError("'parallel_axes' applies only to array input (got SkyCoord elements).")
        if not sky:
            raise ValueError(_NONSKY_REFUSAL.format(what="A SkyCoord"))  # P3 basic-frame rule
        # each element converts independently — no shared-frame requirement (P3)
        rad = np.concatenate([_skycoord_to_rad(el, eff) for el in points], axis=1)
        if rad.size == 0:
            raise ValueError("No points provided (empty input).")  # P17 re-applied
        _check_finite(rad, "the converted SkyCoord coordinates")
        # a one-element SkyCoord sequence is a sequence form, not a single-point form
        return np.ascontiguousarray(rad, dtype=np.float64), False

    # step 4: Quantity
    if isinstance(points, u.Quantity):
        if not sky:
            raise ValueError(_NONSKY_REFUSAL.format(what="A Quantity"))  # P6–P9 basic-frame rule (D12)
        rad_values = _quantity_to_rad_values(points)
        arr, single = _coerce_array(rad_values, naxes, parallel_axes)
        return np.ascontiguousarray(arr, dtype=np.float64), single

    # step 5: non-str sequence whose elements are all Quantity
    if isinstance(points, (list, tuple)) and all(isinstance(el, u.Quantity) for el in points):
        if not sky:
            raise ValueError(_NONSKY_REFUSAL.format(what="A Quantity"))  # P4/P5 basic-frame rule (D12)
        try:
            stacked = u.Quantity(points)  # unifies units (verified V7)
        except u.UnitConversionError as e:
            # a stacking failure is a units VALUE problem (§2.5), not P20
            raise ValueError(f"The elements of the Quantity sequence have incompatible units "
                             f"({[str(el.unit) for el in points]}).") from e
        if stacked.size == 0:
            raise ValueError("No points provided (empty input).")  # P17 re-applied
        rad_values = _quantity_to_rad_values(stacked)
        arr, single = _coerce_array(rad_values, naxes, parallel_axes)
        return np.ascontiguousarray(arr, dtype=np.float64), single

    # step 6: everything else through float coercion
    try:
        raw = np.asarray(points)
    except Exception as e:
        raise TypeError(f"Cannot interpret an object of type '{type(points).__name__}' "
                        f"as coordinate points.") from e  # P20
    if raw.dtype == np.bool_:
        # True silently becoming 1.0 degrees is a bug, not a coordinate
        raise TypeError("Boolean values cannot be interpreted as coordinate points.")  # P20
    if not (np.issubdtype(raw.dtype, np.floating) or np.issubdtype(raw.dtype, np.integer)):
        raise TypeError(f"Cannot interpret an object of type '{type(points).__name__}' "
                        f"(array dtype '{raw.dtype}') as coordinate points.")  # P20
    arr = raw.astype(np.float64)
    if arr.size == 0:
        raise ValueError("No points provided (empty input).")  # P17 re-applied
    _check_finite(arr, "the provided points")  # P18 — native pass-through does not exempt BAD/NaN
    arr, single = _coerce_array(arr, naxes, parallel_axes)
    if sky:
        arr = np.deg2rad(arr)  # bare numbers mean degrees on sky frames (P10–P13)
    return np.ascontiguousarray(arr, dtype=np.float64), single


# ---------------------------------------------------------------------------
# public surface (signatures frozen in SPEC-04A §1.2)
# ---------------------------------------------------------------------------

def to_frame_units(
    points: PointsInput,
    frame: FrameLike,
    *,
    parallel_axes: Optional[bool] = None,
    squeeze: bool = False,
) -> np.ndarray:
    """User-facing coordinates -> the (naxes, npoints) float64 array pyast wants.

    Sky frames (per is_sky(frame)): output is radians. Non-sky frames: values
    pass through in the frame's native units. THE ONLY deg2rad call site in
    cornish, together with to_frame_distance().

    :param points: any form in the §2 decision table
    :param frame: the frame the points are destined for (Frame, SkyFrame,
        FrameSet, Region, or cornish wrapper); for FrameSets the CURRENT frame
        governs, for Regions the encapsulated current frame governs
    :param parallel_axes: orientation override for bare/Quantity *array* input:
        True = rows are axes (naxes, npoints); False = rows are points
        (npoints, naxes); None (default) = infer from shape, with the square
        (naxes, naxes) case resolving to *pairs* per DECISIONS D3.
        Passing a non-None value with any non-array input form raises ValueError.
    :param squeeze: if True, require exactly one point and return shape
        (naxes,) instead of (naxes, 1); more than one point raises ValueError
    :returns: float64 C-contiguous array, shape (naxes, npoints), or (naxes,)
        when squeeze=True
    :raises TypeError: unsupported input type / frame kind (§2.5)
    :raises ValueError: bad value: wrong shape, wrong units, non-finite,
        Ast.BAD, empty, contradiction with parallel_axes, unconvertible
        SkyCoord system (§2.5)
    """
    arr, _ = _convert(points, frame, parallel_axes=parallel_axes)
    if squeeze:
        if arr.shape[1] != 1:
            raise ValueError(f"squeeze=True requires exactly one point "
                             f"(got {arr.shape[1]} points).")
        return np.ascontiguousarray(arr[:, 0])
    return arr


def from_frame_units(
    points: np.ndarray,
    frame: FrameLike,
    *,
    normalize: bool = True,
    bad: Literal["raise", "nan"] = "raise",
) -> np.ndarray:
    """pyast output -> user-facing degrees (sky) / native (non-sky).

    Input must be the pyast layout: shape (naxes, npoints), or (naxes,) for a
    single point. Output is transposed to (npoints, naxes) — or kept (naxes,)
    for 1-D input — as float64. For sky frames values are rad2deg'ed; for
    non-sky frames values pass through. THE ONLY rad2deg call site in cornish,
    together with from_frame_distance().

    :param points: (naxes, npoints) or (naxes,) array from pyast
    :param frame: same resolution rules as to_frame_units
    :param normalize: apply frame.norm() before unit conversion (norm is a
        verified identity for basic frames, so this is safe to leave True
        everywhere; set False only where AST documentation says the values are
        already normalized and the extra call is measured to matter)
    :param bad: policy for values equal to Ast.BAD in the input:
        'raise' (default) raises ValueError naming the offending indices;
        'nan' replaces them with np.nan BEFORE any unit conversion (explicit
        opt-in for masking workflows, e.g. off-sky pixels in all-sky
        projections)
    :raises TypeError: wrong input type / frame kind
    :raises ValueError: wrong shape (first dimension != naxes), or Ast.BAD
        present with bad='raise', or NaN/inf present that did not come from
        bad='nan' handling
    """
    if bad not in ("raise", "nan"):
        raise ValueError(f"'bad' must be 'raise' or 'nan' (got {bad!r}).")

    eff = current_frame_of(frame)  # TypeError for non-FrameLike
    sky = eff.isaskyframe()
    naxes = int(eff.get("Naxes"))

    try:
        raw = np.asarray(points)
    except Exception as e:
        raise TypeError(f"Cannot interpret an object of type '{type(points).__name__}' "
                        f"as a points array.") from e
    if raw.dtype == np.bool_ or not (np.issubdtype(raw.dtype, np.floating)
                                     or np.issubdtype(raw.dtype, np.integer)):
        raise TypeError(f"Expected a numeric array (got dtype '{raw.dtype}').")
    arr = raw.astype(np.float64)

    if arr.ndim not in (1, 2) or arr.shape[0] != naxes:
        # accepting (npoints, naxes) too would re-import the ambiguity the bridge kills
        raise ValueError(f"Expected the pyast layout ({naxes}, npoints) or ({naxes},); "
                         f"got shape {arr.shape}.")

    # order of operations is frozen: BAD handling -> norm -> rad2deg -> transpose.
    # BAD must be resolved before norm()/rad2deg (both corrupt the sentinel, V11).
    arr = arr.copy()  # never mutate the caller's array
    bad_mask = (arr == Ast.BAD)

    # NaN/inf that did NOT come from bad='nan' handling is always an error
    pre_existing_nonfinite = ~np.isfinite(arr)
    if np.any(pre_existing_nonfinite):
        raise ValueError(f"The points array contains non-finite value(s) at index/indices "
                         f"{np.argwhere(pre_existing_nonfinite).tolist()}.")

    if np.any(bad_mask):
        if bad == "raise":
            raise ValueError(f"The points array contains Ast.BAD value(s) at index/indices "
                             f"{np.argwhere(bad_mask).tolist()}; pass bad='nan' to convert "
                             f"them to NaN for masking workflows.")
        arr[bad_mask] = np.nan  # BEFORE any unit conversion (D13)

    if normalize:
        # norm on the object AS PASSED: a Region normalizes via its own
        # encapsulated frame, a FrameSet via its current frame (V3), so
        # normalization can never disagree with provenance. norm propagates
        # NaN unmodified without raising (V17).
        arr = np.asarray(_unwrap(frame).norm(arr), dtype=np.float64)

    if sky:
        arr = np.rad2deg(arr)

    if arr.ndim == 2:
        return np.ascontiguousarray(arr.T, dtype=np.float64)
    return np.ascontiguousarray(arr, dtype=np.float64)


def to_frame_distance(value: Union[float, int, u.Quantity], frame: FrameLike) -> float:
    """Scalar geodesic distance/size in user units -> frame units (radians on sky).

    Sky: Quantity -> .to(u.rad).value (angular units required); bare number ->
    treated as DEGREES and converted. Non-sky: bare number passes through;
    Quantity raises ValueError (native units unknowable — pass .value yourself).

    :raises TypeError: not a real number or Quantity (bool is rejected)
    :raises ValueError: non-scalar Quantity; non-angular or dimensionless
        Quantity on a sky frame; any Quantity on a non-sky frame; non-finite
        value; value == Ast.BAD
    """
    eff = current_frame_of(frame)
    sky = eff.isaskyframe()

    if isinstance(value, bool):
        raise TypeError("A bool cannot be interpreted as a distance.")
    if isinstance(value, u.Quantity):
        if not value.isscalar:
            raise ValueError(f"Expected a scalar distance (got a Quantity of shape {value.shape}).")
        if not sky:
            raise ValueError(_NONSKY_REFUSAL.format(what="A Quantity distance"))
        if value.unit.physical_type == 'dimensionless':
            raise ValueError("Ambiguous input: a dimensionless Quantity would silently read as "
                             "radians; attach an angular unit or pass a bare number (degrees).")
        raw = float(value.value)
        _check_finite(raw, "the provided distance")
        try:
            return float(value.to_value(u.rad))
        except u.UnitConversionError as e:
            raise ValueError(f"Cannot interpret a Quantity with unit '{value.unit}' as an "
                             f"angular distance on a sky frame.") from e
    if isinstance(value, (int, float, np.integer, np.floating)):
        v = float(value)
        _check_finite(v, "the provided distance")
        return float(np.deg2rad(v)) if sky else v
    raise TypeError(f"Expected a real number or Quantity for the distance "
                    f"(got '{type(value).__name__}').")


def from_frame_distance(value: float, frame: FrameLike) -> float:
    """Inverse of to_frame_distance: frame units -> degrees (sky) / native (non-sky).

    :raises ValueError: non-finite or Ast.BAD input (no 'nan' escape here;
        a scalar distance of BAD is always an upstream error — see
        boundingCircle's existing guard for the pattern)
    """
    eff = current_frame_of(frame)
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TypeError(f"Expected a real number for the distance (got '{type(value).__name__}').")
    v = float(value)
    _check_finite(v, "the provided distance")
    return float(np.rad2deg(v)) if eff.isaskyframe() else v


def is_sky(obj: Union[Ast.Object, "ASTObject"]) -> bool:
    """True iff the EFFECTIVE frame of obj is a sky frame. Dispatch table in §3.

    Fixes the verified footgun: Ast.FrameSet.isaskyframe() returns False even
    when the current frame is a SkyFrame, and Ast.Region.isaskyframe() returns
    False even when the region is defined on a SkyFrame (V1).

    :raises TypeError: obj carries no frame (e.g. a plain Mapping, a Channel)
    """
    return current_frame_of(obj).isaskyframe()


def current_frame_of(obj: Union[Ast.Object, "ASTObject"], *, copy: bool = False) -> Ast.Frame:
    """The effective frame of obj as a pyast Frame. Dispatch table in §3.

    LIVENESS WARNING (verified, V2): with copy=False the returned object is
    whatever AST naturally hands back — a LIVE pointer for FrameSet.getframe
    (mutations propagate into the FrameSet), a deep copy for
    Region.getregionframe, and the object itself for a plain Frame. Treat the
    copy=False return as READ-ONLY. Pass copy=True to get a uniformly
    independent deep copy.

    :raises TypeError: obj carries no frame
    """
    ast_obj = _unwrap(obj)
    # dispatch order matters: Region before FrameSet before Frame, because an
    # Ast.Region answers True to isaframe() (V1) and an Ast.FrameSet is also a
    # Frame subclass
    if isinstance(ast_obj, Ast.Region):
        frame = ast_obj.getregionframe()  # the region's CURRENT frame; a deep copy (V2)
    elif isinstance(ast_obj, Ast.FrameSet):  # includes Ast.Plot (V13)
        frame = ast_obj.getframe(Ast.CURRENT)  # LIVE pointer (V2) — read-only doctrine
    elif isinstance(ast_obj, Ast.Frame):
        frame = ast_obj
    else:
        raise TypeError(f"An object of AST class '{ast_obj.Class}' carries no frame.")
    if copy:
        frame = frame.copy()  # one rule, uniform
    return frame


def parse_points_string(s: str) -> np.ndarray:
    """Parse the polygon string grammar into a float ndarray (degrees implied
    by the caller's frame — this function does NO unit handling).

    Accepted: e.g. "((131.7,5.3),(131.8,3.7),(132.5,3.8))" and the
    square-bracket equivalent. Implementation: strip all characters not in
    [0-9 [ ] ( ) space comma . e E + -] (the existing polygon.py scrub), then
    ast.literal_eval, then np.asarray(dtype=float). Result must be rank 1 or 2.

    :raises ValueError: anything that does not survive that pipeline
        (chained from the underlying SyntaxError/ValueError/TypeError)
    """
    if not isinstance(s, str):
        raise TypeError(f"Expected a string of points (got '{type(s).__name__}').")
    scrubbed = re.sub(r'[^\d\[\]\(\) ,\.eE+-]', '', s)
    try:
        parsed = _python_ast.literal_eval(scrubbed)
        result = np.asarray(parsed, dtype=float)
    except (SyntaxError, ValueError, TypeError, MemoryError) as e:
        raise ValueError(f"Could not parse the provided string into an array of coordinate "
                         f"points: {s!r}") from e
    if result.ndim not in (1, 2):
        raise ValueError(f"The points string {s!r} parsed to rank {result.ndim}; "
                         f"expected rank 1 or 2.")
    return result


def astropy_frame_for(frame: FrameLike) -> "astropy.coordinates.BaseCoordinateFrame":
    """The astropy frame instance corresponding to an effective-sky frame's
    AST System/Equinox/Epoch, per the §2.3 map. Used in BOTH directions: as
    the transform_to target for SkyCoord INPUT, and as the frame in which
    SkyCoord OUTPUT conveniences (M24; SPEC-08 asSkyCoord) build their
    results — so returned SkyCoords are in the destination's actual system,
    never ICRS-by-fiat.

    :raises ValueError: is_sky(frame) is False, or the system has no exact
        astropy equivalent (the §2.3 refusal list; same message both directions)
    :raises TypeError: frame is not FrameLike
    """
    eff = current_frame_of(frame)
    if not eff.isaskyframe():
        raise ValueError(f"The effective frame (AST class '{eff.Class}', system "
                         f"'{eff.get('System')}') is not a sky frame; there is no "
                         f"corresponding astropy coordinate frame.")

    def _ast_year_to_time(value: str) -> Time:
        # AST returns bare-number epoch/equinox strings (verified V4); per the
        # SUN/211 convention, < 1984.0 is Besselian, >= 1984.0 is Julian
        year = float(value)
        return Time(year, format='byear') if year < 1984.0 else Time(year, format='jyear')

    system = eff.get("System")  # AST canonicalizes input aliases on set (V4)
    if system == "ICRS":
        return astropy.coordinates.ICRS()
    if system == "FK5":
        return astropy.coordinates.FK5(equinox=_ast_year_to_time(eff.get("Equinox")))
    if system == "FK4":
        return astropy.coordinates.FK4(equinox=_ast_year_to_time(eff.get("Equinox")),
                                       obstime=_ast_year_to_time(eff.get("Epoch")))
    if system == "FK4-NO-E":
        return astropy.coordinates.FK4NoETerms(equinox=_ast_year_to_time(eff.get("Equinox")),
                                               obstime=_ast_year_to_time(eff.get("Epoch")))
    if system == "GALACTIC":
        return astropy.coordinates.Galactic()
    if system == "SUPERGALACTIC":
        return astropy.coordinates.Supergalactic()
    if system == "ECLIPTIC":
        return astropy.coordinates.GeocentricMeanEcliptic(
            equinox=_ast_year_to_time(eff.get("Equinox")))
    # J2000 is a DISTINCT AST system (mean dynamical equatorial), not an FK5
    # alias — deliberately refused rather than approximated (V4)
    raise ValueError(f"cannot convert a SkyCoord into AST system '{system}' (no exact "
                     f"astropy equivalent); transform yourself and pass bare degrees")


def is_single_point(points: PointsInput, frame: FrameLike) -> bool:
    """True iff ``points`` is in a SINGLE-POINT FORM: a scalar SkyCoord, a
    rank-1 bare array/list/tuple or Quantity of shape (naxes,), a sequence of
    exactly naxes scalar Quantities, or a string that parses to rank 1.

    Used by wrapper methods (pix2world/world2pix/convertPoints, pointInRegion
    callers) to preserve single-point-in -> single-point-out return shapes.
    Explicit ARRAY forms containing one point — (1, naxes) or (naxes, 1) — are
    NOT single-point forms (array in -> array out). This predicate is decided
    by FORM, not by npoints.

    Raises exactly what to_frame_units would raise for the same input (it is
    the same classifier); in the wrapper-method pattern it is called first, so
    invalid input fails here with the §2.5 taxonomy.
    """
    _, single = _convert(points, frame)
    return single


def as_uncertainty_region(
    unc: Union[None, float, u.Quantity, Ast.Region, "ASTObject"],
    frame: FrameLike,
    centre: np.ndarray,
) -> Optional[Ast.Region]:
    """Normalize a user 'uncertainty' argument into the Ast.Region that region
    constructors require — or None (AST default uncertainty).

    None -> None. Ast.Region -> itself. cornish ASTRegion -> its .astObject.
    float/Quantity -> Ast.Circle(<unwrapped frame>, 1, centre,
    [to_frame_distance(unc, frame)]): an uncertainty circle of that radius.

    ``centre`` is a (naxes,) point ALREADY in frame units (callers invoke this
    immediately after to_frame_units, passing one of the converted points).
    AST only uses the uncertainty region's extent (it re-centres the region on
    each boundary point), so any in-region point is acceptable.

    CRITICAL implementation rules (finding N1, §8): (a) the returned region
    must be passed to pyast constructors POSITIONALLY — pyast silently ignores
    the ``unc=`` keyword on Circle, Box, and Polygon (verified, V10); (b) when
    this returns None, callers must OMIT the positional argument entirely —
    pyast rejects a positional None (``TypeError: argument 5 must be
    starlink.Ast.Region, not None``; verified, V18). The call-site pattern is
    therefore: build the argument tuple, append the uncertainty region only if
    non-None, splat into the constructor.

    :raises TypeError: unc of any other type
    :raises ValueError: numeric unc rejected by to_frame_distance (incl. any
        Quantity on a non-sky frame), or <= 0
    """
    if unc is None:
        return None
    if isinstance(unc, ASTObject):
        unwrapped = unc.astObject
        if not isinstance(unwrapped, Ast.Region):
            raise TypeError(f"A cornish wrapper passed as an uncertainty must wrap an "
                            f"Ast.Region (got '{type(unwrapped).__name__}').")
        return unwrapped
    if isinstance(unc, Ast.Region):
        return unc
    if isinstance(unc, bool):
        raise TypeError("A bool cannot be interpreted as an uncertainty.")
    if isinstance(unc, (int, float, np.integer, np.floating, u.Quantity)):
        radius = to_frame_distance(unc, frame)  # frame-units policy (D15)
        if radius <= 0:
            raise ValueError(f"An uncertainty must be positive (got {unc!r}).")
        centre = np.ascontiguousarray(centre, dtype=np.float64)
        return Ast.Circle(_unwrap(frame), 1, centre, [radius])
    raise TypeError(f"Cannot interpret an object of type '{type(unc).__name__}' as an "
                    f"uncertainty; expected None, a number, a Quantity, or a Region.")


def dump_value(
    ast_object: Union[Ast.Object, "ASTObject"],
    component: str,
    *,
    include_defaults: bool = True,
) -> str:
    """Read a named top-level component from an AST object's text dump.

    The sanctioned pure-Python escape hatch (02 §3) for values AST computes
    but pyast does not expose. Grammar and worked examples in §7. Component
    names are the labels appearing in the dump — e.g. a CmpRegion's operator
    is ``"Operator"`` (NOT the C attribute name "Oper"; verified V8, SPEC-04
    corrected accordingly).

    :param include_defaults: dump lines prefixed with '#' are values AST
        reports but were not explicitly set; True (default) reads them too,
        False restricts to explicitly-set components
    :returns: for scalar components, the value token with the trailing
        comment stripped and surrounding double quotes removed; for
        object-valued components (e.g. "Unc"), the nested Begin..End block
        text, newline-joined, exactly as dumped.

        CAUTION for object-valued components: the block parses through an
        ``Ast.Channel`` without error, but AST elides any frame the component
        shares with its enclosing object (its normal dump behavior), so the
        read-back object silently carries a generic Cartesian ``Ast.Frame``
        instead of the original frame. Verified for both ``"Unc"`` and
        ``CmpRegion``'s ``"RegionA"``/``"RegionB"``: a sky circle read back
        this way reports ``isaskyframe() == False``, degrades ``overlap()``
        results, AND returns distorted geometry — ``circlepars()`` radius is
        inflated by exactly 1/cos(dec), since Cartesian distance is applied
        to the cos(dec)-scaled stored RA offset. Only the raw stored axis
        values (``getregionpoints()``) survive unchanged. Anything else
        requires re-attaching the enclosing object's frame rather than
        trusting the isolated block. Pinned by
        ``test_dump_value_object_component_loses_shared_frame``.
    :raises KeyError: component not present at the top level of the dump
        (message must name the component and the object class)
    :raises TypeError: not an AST object / cornish wrapper
    """
    ast_obj = _unwrap(ast_object)
    if not isinstance(component, str):
        raise TypeError(f"'component' must be a string (got '{type(component).__name__}').")

    lines = str(ast_obj).splitlines()
    depth = 0
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        i += 1
        if not stripped:
            continue
        is_default = stripped.startswith('#')
        content = stripped.lstrip('#').strip() if is_default else stripped
        first_word = content.split(None, 1)[0] if content else ''
        if first_word == 'Begin':
            depth += 1
            continue
        if first_word == 'End':
            depth -= 1
            continue
        if first_word == 'IsA':
            # facet separator WITHIN a record — depth unchanged (D14)
            continue
        if depth != 1:
            continue  # "top level" = depth 1: inside the outermost record only
        if '=' not in content:
            continue
        name, _, remainder = content.partition('=')
        if name.strip() != component:
            continue
        if is_default and not include_defaults:
            continue  # component names are unique within a record; keep scanning to the KeyError
        # strip the trailing comment: value runs to the tab (or line end)
        value = remainder.split('\t', 1)[0].strip()
        if value:
            # scalar component: remove surrounding double quotes, no type coercion
            if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            return value
        # object-valued component: the object follows as a nested record —
        # return its lines Begin..End inclusive, indentation preserved
        block_lines = []
        block_depth = 0
        while i < len(lines):
            block_line = lines[i]
            block_word = block_line.strip().split(None, 1)[0] if block_line.strip() else ''
            i += 1
            if block_word == 'Begin':
                block_depth += 1
            block_lines.append(block_line)
            if block_word == 'End':
                block_depth -= 1
                if block_depth == 0:
                    return '\n'.join(block_lines)
        raise ValueError(f"Malformed AST dump: the object-valued component '{component}' "
                         f"has no closing 'End' line.")
    raise KeyError(f"Component '{component}' is not present at the top level of the dump "
                   f"of this {ast_obj.Class}.")


def reconstruct(cls: type, ast_string: str) -> "ASTObject":
    """Pickle helper (module-level so it is picklable by qualified name); see §6.

    Reads ast_string back into a pyast object via Ast.Channel +
    channel_io.ListSource, allocates cls via cls.__new__(cls) (bypassing
    __init__ and its validation variance across subclasses), assigns
    .astObject, and returns the instance. Python-side attributes are restored
    afterwards by pickle via the standard __reduce__ state dict.

    :raises ValueError: the string does not read back as an AST object
        (chained from the underlying Ast error). Per the §6 frozen rule and
        D16, NO class-compatibility check is performed beyond a successful
        read — the dump came from a live instance of cls, and tampered
        pickles are out of scope.
    """
    from .channel.channel_io import ListSource  # local import to avoid any import-order issues

    try:
        ast_object = Ast.Channel(ListSource(ast_string), None).read()
    except Ast.AstError as e:
        raise ValueError(f"The provided string could not be read back as an AST object.") from e
    if ast_object is None:
        raise ValueError("The provided string could not be read back as an AST object "
                         "(the channel read returned nothing).")
    obj = cls.__new__(cls)  # deliberately bypasses __init__ (D16)
    obj.astObject = ast_object
    return obj
