# Region algebra, the type registry, and compound repair (SPEC-03)

Implemented on branch `claude/spec-03a` (2026-09-03), against the frozen
`specs/SPEC-03A_algebra_design.md` in the machine-local design corpus
(`local_development/2026-07-09 Claude Fable project review/`). Not merged —
that is the owner's call. Suite 379 → 557 passing, also under
`-W error::DeprecationWarning`.

## What you can now write

```python
a | b        a.union(b, c, d)            # union; n-ary folds to a balanced tree
a & b        a.intersection(b, c)        # intersection
a - b        a.difference(b)             # a with b removed (a & ~b)
a ^ b        a.symmetricDifference(b)    # symmetric difference
~a                                       # negation: same class, not a compound
a.simplified()                           # AST's simplify(), registry-wrapped
circle.dilatedBy(30*u.arcmin)            # also circle + 0.5
circle.scaledBy(1.2)                     # also circle * 1.2, circle / 1.2
```

`+` between two regions **raises** (D4): its old base-class meaning was
intersection, which contradicted `ASTCircle.__add__`'s "grow the radius". Use
`&`. Raw `starlink.Ast` regions define no Python operators, so a raw left
operand reaches the reflected form; every reflected form preserves the order you
wrote, because the result's frame is the first-written operand's and `-` is not
commutative.

`ASTCompoundRegion` gained `components` (the two children, as COPIES),
`flattened()` (all leaves, iteratively, left to right), `operation`, and a real
`toPolygon()`. `componentRegions()` still works as a deprecated delegate.

## Three things that will bite you if you forget them

**1. AST does not refuse a sky × basic region pair.** At equal axis counts it
identity-maps the basic frame's native values *as radians* into the sky frame —
a pixel-frame circle at (50, 50) silently becomes a sky region near (50 rad,
50 rad), and containment, overlap, disc and mesh are all garbage with no
exception anywhere. `overlap()` even answers `NONE` rather than
`NO_FRAME_MAPPING`. Cornish gates this itself, on **every** construction path
including the plain `ASTCompoundRegion(regions=..., operation=...)` constructor,
and `overlapType` answers `NO_FRAME_MAPPING` for a mixed pair without consulting
AST. Operand coercion is region-only: a frame satisfies the gate's own checks,
so it must be rejected as a TYPE error at the boundary.

**2. AST's bounds machinery ignores negation; its membership machinery does
not.** `getregiondisc` and `getregionbounds` answer about the un-negated shape,
so `(~a).bounds` used to hand back `a`'s own box while `~a` covers the rest of
the sky. `.bounds`, `boundingBox()` and `boundingCircle()` now probe just
outside the bound they are about to return and raise if the region has members
there. `ASTCircle.boundingCircle()`'s `return self` shortcut carried the same
hole and now refuses a negated circle.

**3. Tree shape is not cosmetic.** Every `Ast.CmpRegion` deep-copies its
operands, so a left fold re-copies the growing tree. Measured at n = 1000:
1.34 s to build left-folded versus 37 ms balanced, and `overlap`/`getregiondisc`
blow a two-minute budget versus ~8 s. `pointinregion` is shape-INsensitive, so a
regression test must watch build and overlap, not containment.

## AST defects found while building this

All verified live against AST 9.3.0 / starlink-pyast 4.0.1 / Python 3.13. The
first three are recorded in the frozen spec's Appendix B as corrections to it.

| # | defect | how it shows |
|---|---|---|
| 1 | `Ast.CmpRegion` **recognizes** the composition `(a & ~b) \| (~a & b)` and stores it as a native XOR node | so cornish cannot avoid building native XOR, `(a ^ b).operation` reads XOR, and the spec's plan for conformant STC-S output does not work |
| 2 | an XOR node **nested** inside another compound serializes to STC-S as a plain `Union` of its operands | `(a ^ b) \| c` writes as `Union ( a b c )`, a different region — membership differs at the a∩b overlap. A post-write text guard misses it; the guard here is structural |
| 3 | `getregionmesh(SURFACE)` is **negation-blind for `Ast.Box`** — and only for Box | a negated box's interior mesh is the un-negated box's, 0 of 195 points of which are members; Circle, Polygon and CmpRegion all honour negation |
| 4 | AST's **spherical-polygon `pointinregion` has isolated, measure-zero True answers outside the polygon** | one sits on the ordinary FITS-field polygon of a 512×512 TAN header: (31.357°, 44.016°) reads True for a 0.36°-radius polygon while every neighbour 0.02° away reads False. This is why every escape probe is corroborated by four neighbours |
| 5 | AST reads and writes STC-S `Difference` as the **symmetric** difference | IVOA STC 1.33 defines it as asymmetric (`R1 AND NOT R2`) and says explicitly it is not symmetrical. `toSTCS()` refuses to emit it; `fromSTCS()` warns (`STCSDifferenceWarning`) and still parses. Ready-to-post report: `local_development/ast_stcs_difference_bug_report.md` |

Pre-existing cornish bugs the spec had already identified and this work fixed:
leaf `toPolygon()` was negation-broken in **both** leaf classes (`ASTCircle` had
no orientation logic at all; `ASTBox`'s `containsPoint(centre)` guard tested the
freshly built, never-negated polygon against the box's own structural centre and
so was trivially always True), and `ASTFrame.isSkyFrame` answered `False` for
every `ASTFrameSet` on a sky frame and every `ASTRegion` defined on one.

## `toPolygon()` on a compound: what it refuses, and why

The boundary mesh is split at its seams into per-component arcs, the arcs are
chained into loops, and the single loop becomes the polygon. Every step is a
heuristic over a finite sample, so the output is CHECKED, not trusted: the
compound's own interior mesh must lie inside the polygon and its complement's
outside it. Without that post-condition an ordinary two-circle union at
`npoints=8` returns a polygon that disagrees with the compound almost
everywhere.

It raises rather than guessing when:

- the compound is empty. Two conduct classes: the disjoint-AND class returns a
  BAD bounding disc, but `a - a` / `a ^ a` return a **finite** disc of about the
  operand's scale and a full, ordinary-looking boundary mesh — only the SURFACE
  mesh's `Ast.NPTIN` gives it away;
- the boundary is more than one loop, or is pinched at a crossing point
  (`RegionNotConnected`) — but only when the arcs are resolved well enough to
  support that claim; under-resolved sampling says "increase npoints" instead of
  asserting a topological fact about a region that may well be connected;
- the sample cannot support verification at all.

Known limitations, loud rather than silent: a union of two boxes sharing a
corner cannot be chained (their arcs do not meet endpoint-to-endpoint), and a
very thin lens can read as pinched at a given `npoints`.

Note `MeshSize` applies **per component**, so an n-leaf union at `npoints=200`
yields roughly n × 200 boundary points, and setting it on a parent propagates
into its children (they are live).

## Out of scope here

`ASTMoc.add` can silently contribute nothing when the region is small relative
to the MOC cell scale (MaxOrder-dependent, not negation-specific), and
`MaxOrder=0` + `addregion` segfaults the interpreter. Both are real and both
belong to SPEC-01, not this work — which is why `boundingCircle()`'s error text
deliberately does **not** point users at `toMoc()` for negated geometry.

## Review state

Silent-skip audit clean. Sonnet 3 rounds (round 1 found a dead safety-net branch
whose test then exposed a real pooling bug; rounds 2 and 3 dry). Opus 4 rounds,
each finding something real about the verification post-condition, closing with
sonnet's independent dry round over the whole chain. Codex was unavailable
(usage limit) for the whole run and was substituted by sonnet per the
orchestrator's ruling — a codex reading of this branch is still owed.
