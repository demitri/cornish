# START HERE

Cornish is a Pythonic OO wrapper around Starlink AST (via `starlink-pyast`), currently focused on regions and actively growing broader AST function coverage.

## What this is

Cornish wraps the C-flavored `starlink-pyast` interface to the Starlink AST library (world coordinate systems, coordinate transforms, regions on the celestial sphere, plotting) in a proper Pythonic object model. It accepts and returns Astropy objects where possible. The initial impetus was region support (polygons, circles, boxes, compound regions), but the broader goal is a complete OO wrapper around AST's C API surface — extending use cases and function coverage over time. It's a dependency of the Trillian and SciDD projects; not yet recommended for pipeline use, but increasingly mature for interactive use.

## Current state

- Regions (box, circle, polygon) are fairly mature and the primary use case so far.
- Compound regions are actively in progress — not recommended to use yet (see recent commits and `AST+Cornish Wish List.md`).
- `addPoints` API is in flux — still being reworked, not settled.
- FITS header handling (`ASTRegion.fromFITSHeader`, `channel/fits_channel.py`) is in active use.
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
