# START HERE

Cornish is a Pythonic OO wrapper around Starlink AST (via `starlink-pyast`), currently focused on regions and actively growing broader AST function coverage.

## What this is

Cornish wraps the C-flavored `starlink-pyast` interface to the Starlink AST library (world coordinate systems, coordinate transforms, regions on the celestial sphere, plotting) in a proper Pythonic object model. It accepts and returns Astropy objects where possible. The initial impetus was region support (polygons, circles, boxes, compound regions), but the broader goal is a complete OO wrapper around AST's C API surface — extending use cases and function coverage over time. It's a dependency of the Trillian and SciDD projects; not yet recommended for pipeline use, but increasingly mature for interactive use.

## Current state

- Regions (box, circle, polygon) are fairly mature and the primary use case so far.
- Compound regions are actively in progress — not recommended to use yet (see recent commits and `AST+Cornish Wish List.md`). Note: several wish-list blockers are fixed upstream in AST 9.3 (`getregiondisc` and boundary meshes now work on compound regions — verified 2026-07-09); converting to a MOC (`region.toMoc()`) is the recommended way to get bounding circles/areas for compound geometry.
- `ASTMoc` (IVOA Multi-Order Coverage maps) and STC-S serialization (`region.toSTCS()` / `ASTRegion.fromSTCS()`) added 2026-07-09, with tests in `source/tests/test_moc_stcs.py`.
- **Branch `spec-04-bridge` (2026-07-10) implements SPEC-04 + SPEC-10**: the private `cornish/_pyast_bridge.py` module is now the ONLY place unit/shape conversions happen (machine-enforced by `test_no_conversions_outside_bridge`); `cornish/enums.py` holds AST-welded enums; all exceptions descend from `cornish.exc.CornishError` (no bare `Exception` raised — gate-enforced); regions pickle via `__reduce__`; `fromFITSHeader` has a fast mapregion path and carries `.wcs`. Suite 320/320 on the branch. Awaiting final review rounds + merge — read the review-state note in the local `ORCHESTRATION.md` before touching it.
- `addPoints` accepts the full bridge input contract on the branch (SkyCoord/Quantity/bare degrees).
- FITS header handling (`ASTRegion.fromFITSHeader`, `channel/fits_channel.py`) is in active use. On the branch, `ASTRegion.fromFITSHeader` works again as a dispatcher (the README's first example is valid there); on master it still raises "deprecated".
- A full project evaluation + implementation-ready spec corpus (API gaps, C-vs-pyast verdict, recipes, plotting/`cornish-view` design) lives in `local_development/2026-07-09 Claude Fable project review/` (gitignored; machine-local). **Resuming that effort? Read its `ORCHESTRATION.md` first** — it tracks what's done and the next action.
- The pyast `unc=`-keyword-silently-ignored bug (SPEC-04A §8): FIXED on the branch — uncertainties are delivered positionally and a 1-arcsec default is materialized on sky frames; on master no cornish-built region has ever carried its stated uncertainty to AST. Upstream bug report still to be filed (see local ORCHESTRATION loose ends).
- Version: see `source/cornish/version.py` (currently 1.1).

## Where things are

- `source/cornish/region/` — region classes (`region.py`, `box.py`, `circle.py`, `polygon.py`, `compound_region.py`)
- `source/cornish/channel/` — AST channel I/O, including FITS header handling (`fits_channel.py`)
- `source/cornish/mapping/` — coordinate mapping/frame code
- `source/cornish/plot/` — plotting support (matplotlib)
- `source/tests/` — pytest suite (frames, polygons, fitschan)
- `docs/` — Sphinx docs source (published at readthedocs.io); `api_regions.rst`, `api_mapping_and_frames.rst`, `api_plotting.rst`
- `AST+Cornish Wish List.md` (repo root) — open work items / known gaps (NaN handling in `boundingCircle`, compound region support, exposing `getregionframeset`); treat this as the TODO list until it's folded into `AI/`

## Conventions a fresh session would otherwise violate

- **All coordinate unit/shape conversions go through `cornish/_pyast_bridge.py`** (`to_frame_units`/`from_frame_units` and the distance pair). Never write `deg2rad`/`rad2deg`/`.to(u.rad)` in package code — a test gate fails the suite if you do (angle-of-result lines are the only allowlisted exception).
- **Pass ALL arguments to pyast constructors POSITIONALLY** — pyast silently swallows unknown keyword arguments (verified victims: region `unc=` dropped, CmpRegion `oper=` always-OR, FrameSet `frame=` discarded). A test gate (`test_no_kwargs_into_pyast_constructors`) fails the suite on any kwarg into an `Ast.<Class>(...)` call. Uncertainty additionally is appended only when non-None (pyast rejects a positional `None`; D15) — the call sites carry comments saying so; do not "clean them up".
- **Values from fixed sets use the enums in `cornish/enums.py`** (OverlapType, RegionOperation, MeshType, …), not bare strings/ints; `cornish/constants.py` is a deprecated alias layer.
- **No bare `raise Exception`** — TypeError for wrong types, ValueError for bad values, `cornish.exc` subclasses for domain errors (SPEC-10 policy; AST-walk gate enforces the unraised-exception pattern too).
- Test channel/`FitsChan` objects with `is not None`, never truthiness — an empty `Ast.FitsChan` is falsy.
