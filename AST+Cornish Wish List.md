# AST/Cornish Wish List

*Status review 2026-07-09: several items below were re-tested against AST 9.3.0 (starlink-pyast 4.0.1); results noted inline.*

### ~~Error Handling of NaN returned from `boundingCircle`~~ — resolved 2026-07-09

* The actual mechanism: AST signals "no disc" with `Ast.BAD` (a large *finite* float, −1.8e308), not NaN. `ASTRegion.boundingCircle()` now checks for both and raises `ValueError` for empty/unbounded regions (e.g. a MOC of the AND of disjoint regions) instead of letting bad values propagate. Regression test: `test_moc_stcs.py::test_empty_moc_bounding_circle_raises`.

### Compound Region Support

* ~~`getregiondisc` returns NaN values for a compound region~~ — **fixed upstream**: returns finite centre/radius on AST 9.3 (verified with a two-circle union).
* ~~`compund_region.boundaryPointMesh` returns an empty array~~ — **fixed upstream**: returns a real mesh on AST 9.3 (268 points for a two-circle union).
* How are component regions accessed? — `ASTCompoundRegion.componentRegions()` (repaired 2026-07-09; was returning `regions()`, a call on a list).
* Bounding polygon? / Bounding circle? — recommended path: `compound.toMoc().boundingCircle()` (robust; exact for the MOC's cell set, which approximates the compound geometry at the chosen order); `boundingCircle()` directly on the compound also works on AST 9.3.
* Plotting a compound region? — `regionoutline` accepts any `Ast.Region`; should work now that meshes do (untested).

### Support for `region.frameset`

* This requires the C function `getregionframeset` to be exposed to the Python API.
* 2026-07-09 finding: on AST 9.3 the FrameSet encapsulated in a region produced by `mapregion` is simplified to identity (verified via object dump), so the C function would mostly return a trivial mapping. The useful fix is cornish-side: retain the originating pixel↔sky `ASTFrameSet` as an attribute (e.g. `region.wcs`) when a region is created from a FITS header.

