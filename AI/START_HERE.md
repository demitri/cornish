# START HERE

Cornish is a Pythonic OO wrapper around Starlink AST (via `starlink-pyast`), currently focused on regions and actively growing broader AST function coverage.

## What this is

Cornish wraps the C-flavored `starlink-pyast` interface to the Starlink AST library (world coordinate systems, coordinate transforms, regions on the celestial sphere, plotting) in a proper Pythonic object model. It accepts and returns Astropy objects where possible. The initial impetus was region support (polygons, circles, boxes, compound regions), but the broader goal is a complete OO wrapper around AST's C API surface — extending use cases and function coverage over time. It's a dependency of the Trillian and SciDD projects; not yet recommended for pipeline use, but increasingly mature for interactive use.

## Current state

- Regions (box, circle, polygon) are fairly mature and the primary use case so far.
- Compound regions are actively in progress — not recommended to use yet (see recent commits and `AST+Cornish Wish List.md`). Note: several wish-list blockers are fixed upstream in AST 9.3 (`getregiondisc` and boundary meshes now work on compound regions — verified 2026-07-09); converting to a MOC (`region.toMoc()`) is the recommended way to get bounding circles/areas for compound geometry.
- `ASTMoc` (IVOA Multi-Order Coverage maps) and STC-S serialization (`region.toSTCS()` / `ASTRegion.fromSTCS()`) added 2026-07-09, with tests in `source/tests/test_moc_stcs.py`.
- `addPoints` API is in flux — still being reworked, not settled.
- FITS header handling (`ASTRegion.fromFITSHeader`, `channel/fits_channel.py`) is in active use. (`ASTRegion.fromFITSHeader` itself raises "deprecated" — use `ASTPolygon.fromFITSHeader`; the README's first example needs updating.)
- A full project evaluation + implementation-ready spec corpus (API gaps, C-vs-pyast verdict, recipes, plotting/`cornish-view` design) lives in `local_development/2026-07-09 Claude Fable project review/` (gitignored; machine-local). **Resuming that effort? Read its `ORCHESTRATION.md` first** — it tracks what's done and the next action. The bridge-module design (`specs/SPEC-04A_bridge_design.md`) is frozen and double-reviewed to dry (2026-07-09); SPEC-04 implementation is unblocked.
- Known latent issue recorded there (SPEC-04A §8): pyast silently ignores the `unc=` keyword on region constructors, so no cornish-built region has ever carried its stated uncertainty to AST — the fix ships with the SPEC-04 migration; until then treat documented "1 arcsec uncertainty" claims as aspirational.
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

- None identified yet — add here as they come up.
