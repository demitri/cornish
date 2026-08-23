# Plotting v2 and `cornish-view` (SPEC-09 stages 0–1)

Implemented on branch `claude/spec-09a` against the frozen design
`local_development/2026-07-09 Claude Fable project review/specs/SPEC-09A_plotting_design.md`.
That document remains the source of truth; this file records what a fresh
session needs in order not to break it.

## The one architectural idea

`SkyPlot` v2 drives `Ast.Plot` through **`GrfSequence`**, a recording Grf
backend, and captures the primitive stream. SVG, JSON, the browser and
matplotlib are all *replays* of that stream — there is no backend option (D17).

One backend means one metrics source and one geometry truth, so a plot is
identical in a file, a notebook, a browser and a matplotlib window. AST does all
layout server-side; no consumer measures text or places labels.

    SkyPlot ──▶ Ast.Plot ──▶ GrfSequence ──▶ commands ──▶ envelope/document
                                                            ├─▶ render_svg  (pure)
                                                            ├─▶ cornish-view.js (browser)
                                                            └─▶ to_matplotlib

## Module map

| module | role |
|---|---|
| `plot/serial.py` | `GrfSequence` — the recording Grf backend |
| `plot/envelope.py` | schema v1.0 codec: validators, RFC 8785 canonicalization, hashes |
| `plot/svg.py` | `render_svg` (pure) + the stdlib greyscale PNG encoder |
| `plot/sky.py` | `SkyPlot` v2, the synthetic-header builder, `show()`/live mode |
| `plot/palette.py`, `plot/metrics.py` | colour grammars; renderer-free font metrics |
| `view/registry.py` | `PlotRegistry` — transport-free plot state + subscriptions |
| `view/server.py` | the stdlib daemon (`python -m cornish.view`), SSE, tokens |
| `view/blueprint.py` | the Flask blueprint — the only module importing Flask |
| `view/static/` | `cornish-view.js`, `viewer.html`, `envelope.schema.json`, the node suite |

## Conventions a fresh session would otherwise violate

- **No transport imports under `cornish/plot/`**, module-level or function-local.
  Gate: `test_plot_package_imports_no_transport`. Exactly two function bodies are
  allowlisted (`SkyPlot.show`, `SkyPlot._live_post`); a third fails the gate.
  `cornish.view` imports `cornish.plot`, never the reverse.
- **`render_svg` takes a DOCUMENT dict**, never a bare envelope. It is pure and
  byte-stable — the ETag pattern depends on that. Do not introduce
  nondeterminism (dict iteration order, timestamps, `id()`).
- **Canonical JSON is RFC 8785 (JCS), not `json.dumps(sort_keys=True)`.**
  Python and JS disagree on the formatting of nearly every envelope number
  (`1.0` vs `1`, `1e-07` vs `1e-7`), which would make pinned hashes
  irreproducible in the JS consumers. `_jcs_number` implements ECMAScript
  `Number::toString`; it is verified against a real JS engine.
- **Two distinct colour grammars.** `normalize_colour` for user parameters;
  the wider `normalize_palette_entry` for palette entries, which must also admit
  `transparent`/`rgb()` because `DEFAULT_PALETTE` slot 0 is `rgba(0,0,0,0)`.
  Slot 0 deliberately fails the user grammar. Do not merge them.
- **`astSet` accepts only an INTEGER for a `Colour(...)` attribute.**
  `Colour(Border)=red` and `=#4a7f7b` are both rejected (AST falls through to the
  SkyFrame). CSS colours must go through `GrfSequence.ColToInt` first — which is
  also what makes them serialize as literal CSS strings per D18.
- **`Ast.BAD` is a finite float** (−1.797e308), so it passes any ordinary
  finiteness check. `GrfSequence` refuses it explicitly rather than recording it.
- **The §3.3 fallbacks are behaviours, not rejections** — out-of-range colour →
  foreground, style → solid, font → normal, marker → dot, unknown op → skipped.
  Every one emits BOTH §3.7 channels (SVG comment + logged warning). Never make
  one silent, and never turn one into a hard error: they are what let a 1.0
  consumer survive a 1.1 producer.
- **Golden fixtures are never refreshed silently.** `--regenerate-goldens` only,
  and read the printed diff — a changed label set means the layout changed.
  Regenerating invalidates astplot-js's vendored copy.
- **`envelope.schema.json` is hand-maintained**, held in sync by
  `test_envelope_schema.py` (schema and codec must reach the same verdict on the
  whole invalid-body battery). If you change the codec, change the schema.

## Findings from implementation worth keeping

- **An auto-generated `plot_id` is NOT emitted** in `document()`/`envelopes()`;
  an explicit one is. An embedded auto id made `curl -d @saved.json
  .../plot/demo` fail with 409 against any URL but its own — the §9 exit
  criterion — and made two structurally identical plots hash differently, so
  ETags never matched. Golden fixtures likewise pin no plot id.
- **pyast defect**: for some invalid projection codes (`ZZZ` reproduces it,
  `QQQ` does not) `FitsChan.read()` raises `UnicodeDecodeError` from AST's own
  byte string. It is a `ValueError` subclass, so it would otherwise pass a
  `ValueError` test while telling the user nothing about projections. Mapped in
  `SkyPlot._frameset_from_cards`. **Candidate for the upstream pyast report**
  alongside the `unc=` kwarg bug.
- **`astropy.coordinates` imports `http`/`socket`/`urllib.request`** at import
  time, so the runtime half of the purity gate cannot assert their absence; it
  checks `cornish.view`/flask/aiohttp instead. Pinned by its own test so the
  narrowing stays justified.
- **A dead SSE viewer is not detected by a failed write** — writing to a
  half-closed socket succeeds (verified on macOS), so the subscription would
  linger until the next 15 s ping. The server polls the read side for EOF.
- **The default compass placement is a knife-edge by construction** (anchor 12%
  in, arms 12% long → the tip lands ON the boundary; curvature pushes it 0.012%
  past). Tips may touch the edge within 0.5%; the refusal path still fires for
  all-sky corners and over-long arms.
- **The `/plot/{id}` route regex must match the id lazily** — the legal id
  charset contains `.`, so a greedy group swallows `.json`/`.svg`.

## Testing

`source/tests/`: `test_grf_conformance.py` (§7.1), `test_envelope.py` (§7.2),
`test_envelope_schema.py`, `test_golden_streams.py` (§7.3),
`test_render_svg.py` (§7.4 + the purity gate), `test_view_daemon.py` (§7.5),
`test_skyplot_v2.py` (§7.7), `test_flask_example.py`, `test_js_conformance.py`.

The JS suite (`cornish/view/static/cornish-view.test.mjs`, `node --test`) is run
BY the Python suite, so a JS regression fails the build a Python developer
actually runs. It is skipped, not passed, when node is absent.

## Not yet done

- **Stage 2**: `view/plotspec.py` (the declarative `plot_spec` v1.0 endpoint)
  and the `POST /plotspec/{id}` route. SPEC-09A §5.5 freezes the schema; the
  fable-side schema freeze for it is still an open board item.
- **Stage 3**: image regridding (`addImage(wcs=...)` via `Mapping.resample`),
  colour LUTs beyond `greys`/`greys_r`.
- The §8 astplot-js change proposals (marker geometry in `grfMark`, fixture
  vendoring, the recording bridge that would make "two producers" real).
