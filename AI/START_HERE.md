# START HERE

Cornish is a Pythonic OO wrapper around Starlink AST (via `starlink-pyast`), currently focused on regions and actively growing broader AST function coverage.

## What this is

Cornish wraps the C-flavored `starlink-pyast` interface to the Starlink AST library (world coordinate systems, coordinate transforms, regions on the celestial sphere, plotting) in a proper Pythonic object model. It accepts and returns Astropy objects where possible. The initial impetus was region support (polygons, circles, boxes, compound regions), but the broader goal is a complete OO wrapper around AST's C API surface — extending use cases and function coverage over time. It's a dependency of the Trillian and SciDD projects; not yet recommended for pipeline use, but increasingly mature for interactive use.

## Current state

- Regions (box, circle, polygon) are fairly mature and the primary use case so far.
- Compound regions are actively in progress — not recommended to use yet (see recent commits and `AST+Cornish Wish List.md`). Note: several wish-list blockers are fixed upstream in AST 9.3 (`getregiondisc` and boundary meshes now work on compound regions — verified 2026-07-09); converting to a MOC (`region.toMoc()`) is the recommended way to get bounding circles/areas for compound geometry.
- `ASTMoc` (IVOA Multi-Order Coverage maps) and STC-S serialization (`region.toSTCS()` / `ASTRegion.fromSTCS()`) added 2026-07-09, with tests in `source/tests/test_moc_stcs.py`.
- **SPEC-04 + SPEC-10 are MERGED to master (2026-07-10)**: the private `cornish/_pyast_bridge.py` module is the ONLY place unit/shape conversions happen (machine-enforced by `test_no_conversions_outside_bridge`); `cornish/enums.py` holds AST-welded enums; all exceptions descend from `cornish.exc.CornishError` (no bare `Exception` raised — gate-enforced); regions pickle via `__reduce__`; `fromFITSHeader` has a fast mapregion path and carries `.wcs`. Suite 379/379; review cycle closed dry across codex/sonnet/opus.
- `addPoints` accepts the full bridge input contract (SkyCoord/Quantity/bare degrees).
- FITS header handling (`ASTRegion.fromFITSHeader`, `channel/fits_channel.py`) is in active use; `ASTRegion.fromFITSHeader` works again as a dispatcher (the README's first example is valid).
- A full project evaluation + implementation-ready spec corpus (API gaps, C-vs-pyast verdict, recipes, plotting/`cornish-view` design) lives in `local_development/2026-07-09 Claude Fable project review/` (gitignored; machine-local; has its own local git repo). **Resuming that effort? Read its `ORCHESTRATION.md` first** — it tracks what's done and the next action. Both design freezes are done (SPEC-04A 2026-07-09, SPEC-09A plotting/cornish-view 2026-07-10, each reviewed to dry) — implementation phases are what remain.
- The pyast kwarg-swallowing bug class (SPEC-04A §8; `unc=`/`oper=`/`frame=` silently dropped): FIXED here — all pyast constructor arguments are passed positionally (gate-enforced), uncertainties are delivered and a 1-arcsec default is materialized on sky frames. Before 2026-07-10 no cornish-built region ever carried its stated uncertainty to AST. Upstream bug report still to be filed (see local ORCHESTRATION loose ends).
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
- Test channel/`FitsChan` objects (and parameters generally) with `is not None`, never truthiness — an empty `Ast.FitsChan`/astropy `Header` is falsy, and falsy-but-provided parameters must conflict or fail, never be silently ignored.
- **Wrapper setters route user values to AST through `ASTObject._setAttribute`** (quotes the value — astSet parses comma-separated settings lists — and chains `Ast.AstError` into `ValueError`). Never call `astObject.set()` with a user-provided value directly.
- **Integer-valued parameters route through `cornish/_validation.as_integer`** (NumPy ints accepted, bools/floats rejected loudly). AST index parameters are 1-based; `Ast.BASE` (0) and `Ast.CURRENT` (-1) are AST sentinels, not Python-style indices.
