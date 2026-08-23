# Fable orchestration session for the cornish modernization effort: design corpus, review cycles, and phase coordination.

Parked: 2026-07-30 America/New_York

STATUS: This was the original overnight-review session (2026-07-09) that grew into the effort's standing orchestrator. All of its own work is complete and committed: the design corpus (00–05, DECISIONS D1–D30+, HANDOFF, ORCHESTRATION, SPEC-01…11 + 07B/07C + 04A/09A/03A), Phase 0/0.5/1 closed, Phase 2 row 1 (SPEC-04+10) merged to master. NOT done by this session and possibly done by others in the 17-day gap: Phase 2 rows 2 (SPEC-03) and 5 (SPEC-09 stages 0–1) implementation, the recipe waves, and resolution of the D9-amendment proposal. This session held the fable-side duties: review adjudication, spec authorship, and board maintenance.

## Last session output (verbatim)

> SPEC-03A is on the board with an honest caveat worth pausing on. Current state:
>
> **Landed**: `SPEC-03A_algebra_design.md` frozen (session 4), 12 probe scripts / 32 verified pyast behaviors preserved, D24–D30 appended — and the cycle found **four live master bugs** along the way (leaf `toPolygon()` negation-broken on Circle and Box, `boundingCircle()` silently bogus on `~a`/`a^a`, an `ASTMoc.add` silent case), frozen-for-repair in the spec. Rows 2 (SPEC-03) and 5 (SPEC-09 stages 0–1) are both cleared for opus, parallel-safe.
>
> **The caveat the board records honestly**: you terminated the cycle at sonnet round 23 with the codex phase waived — and rounds 20–21 were *still yielding P1s*. So SPEC-03A is very-heavily-reviewed but not certified dry, and the one reviewer family that historically finds what sonnet passes (codex found 14 post-sonnet-dry findings on SPEC-04A, blockers included) never read it. I won't countermand your explicit waiver by launching codex myself — but I'd be failing the adjudication job if I didn't say: a *single* codex round on SPEC-03A is the cheapest risk reduction available before opus builds on it. Your call; it costs codex quota, zero fable.
>
> **The process signal underneath your frustration is real, and I'd rather fix the protocol than defend round 23.** Three cycles of data now: SPEC-04A converged in 6+4 rounds *because* it switched reviewer families; SPEC-09A took 18 sonnet rounds; SPEC-03A hit 23 without drying. Same-family full-read rounds on 100KB documents have an asymptotic P3 tail — the dry-round rule, written for code diffs, over-fires there. I've put a concrete proposal on the board for your sign-off (it changes D9, so it's yours to ratify):
>
> - **Family-switch is mandatory after 6 rounds** of one reviewer on a document (empirically where same-family yield degrades but cross-family yield doesn't);
> - **Document cycles terminate on severity-trajectory**: two consecutive rounds with nothing ≥ P2 = effective-dry for prose artifacts (code keeps the strict dry rule);
> - **Hard cap 10 rounds per family**, then escalate to you with the findings log instead of grinding.
>
> **Next step**: launch opus on row 2 (`Implement SPEC-03 from …/specs/SPEC-03A_algebra_design.md`) and row 5 (`Implement SPEC-09 stages 0–1 from …/specs/SPEC-09A_plotting_design.md`) — separate branches, both prompts are one-liners now, and the D9 implementation reviews (which stay strict-dry, being code) will catch anything the truncated document cycle left behind. Fable then has nothing queued until adjudication is needed or wave A lands.

*(Note added at park time: the codex-round item above resolved in the gap — ORCHESTRATION row 2 records a codex second-reader cycle run 2026-07-26, 3 rounds, closed dry.)*

## Resume context

**This session's role** — fable-side orchestration for the cornish effort: author/adjudicate design specs, run/adjudicate review cycles, keep the board current. Implementation is opus work; routine review rounds are sonnet/codex; this session acts only where judgment or taste is needed.

**Boot documents (read in this order; the board outranks this file for project position):**
1. `local_development/2026-07-09 Claude Fable project review/ORCHESTRATION.md` — the live status board; execute its first unticked item in the fable-shaped list. (The corpus directory is gitignored from cornish but is its own git repo — check `git -C local_development log` for changes by other sessions.)
2. `local_development/2026-07-09 Claude Fable project review/05_review_log.md` — full review history.
3. `HANDOFF.md` + `DECISIONS.md` (same folder) — binding build order, trap list, D-decisions.
4. `AI/START_HERE.md` — repo-level state notes other sessions maintain.

**Binding rules and gates:**
- D9 review protocol (per-implementation: silent-skip audit → codex to dry → sonnet to dry → opus if real findings), board-tick convention (tick + timestamp on every completed step), corpus commits carry the Claude-Session trailer, cornish master commits/pushes are Demitri's call.
- **The user's global review policy was updated (~2026-07-24) after this session's practices**: review prompts must now be the minimal generated form (`Review <hash>` + optional one-line) — this session's seeded sonnet prompts are cited in the global CLAUDE.md as the failure mode; and **one session per reviewer** (codex-reply into the same thread across rounds, same for sonnet). Follow the new policy, not this session's older habits.
- Fable is now subscription-standard with **no end date** — the "fable window / allocation scarcity" framing throughout this session and some board notes is obsolete. Global rule: remove any verbatim fable end-date claims from project docs when encountered.

**Environment facts:**
- `[stable]` Repo: `$GH/cornish`, master, AI/ tracked. Corpus: `local_development/2026-07-09 Claude Fable project review/` (gitignored; own git repo, branch main). Siblings: `$GH/astplot-js` (Phase B layer system landed), `$GH/fitsjs` (early). pyast 4.0.1 / AST 9.3 at `/usr/local/anaconda/bin/python3`.
- `[re-verify]` Phase 2 rows 2 (SPEC-03) and 5 (SPEC-09 stages 0–1) implementation status — cleared as of 2026-07-13, codex-dry on 03A as of 07-26, but 17 days of other-session activity is unknown to this session. Check the board + `git log` both repos.
- `[re-verify]` Whether the D9-amendment proposal (verbatim in "Last session output" above) is still pending on the board, was ratified, or was superseded by the 07-24 global policy update.
- `[re-verify]` Whether the pyast upstream bug report (`local_development/pyast_unc_kwarg_bug_report.md`) was filed (Demitri's action item since 07-10).

**Done — do not redo:** the entire design corpus and its review cycles (see board Phase 0–1, all ticked, and 05_review_log); the SPEC-07 recipe catalogs R1–R43 (ideation is complete; do not re-brainstorm); the C-vs-pyast evaluation (verdict: stay on pyast); the Grf/plotting architecture proofs (`appendix_data/`).

**Fable-shaped queue (from the board, as of park):** adjudication on-call for Phase 2 review disagreements · recipe-surface coherence pass AFTER SPEC-07 wave A lands · SPEC-09 stage-2 plot-spec schema freeze when that row comes up · cross-repo arbitration if the SPEC-09A §8 astplot-js proposals hit resistance · mid-course corpus adjudication after ~3 more Phase 2 rows land (specs-vs-reality drift, strained D-decisions, rule-of-two check).

**Pending user decisions (re-ask only if the board still shows them open):**
1. The D9 amendment for document review cycles — three bullets, verbatim in "Last session output" above; awaiting ratification.
2. (Resolved in gap: codex round on SPEC-03A — ran 2026-07-26, dry.)

## Working tree state at park

- `$GH/cornish` (master): clean.
- `$GH/cornish/local_development` (corpus repo, main): clean.
- Possibly touched / not verified: `$GH/astplot-js`, `$GH/fitsjs` — read-only in this session (orientation docs and source reads; no edits made by this session).

## First actions on resume

1. Read `ORCHESTRATION.md` (boot doc 1) — it is authoritative over this file for project position; note its "Last updated" line.
2. Re-verify the three `[re-verify]` items (Phase 2 row status via board + `git log` on both repos; D9-proposal status; pyast bug-report filing).
3. If the D9 proposal is still open, re-present it to Demitri (verbatim from "Last session output") noting the 07-24 global policy update may partially supersede it.
4. Then take the first applicable item from the fable-shaped queue; work per the CURRENT global review policy (minimal generated prompts, one session per reviewer), not this session's older prompt style.
5. Cleanup rule: when the resumed session completes the orchestration effort (all Phase 2 rows closed or the effort is formally wound down), remove ONLY the registry entry whose Details link exactly matches this file's name, then delete this file. If the effort is superseded or abandoned instead, do NOT silently delete — update STATUS here and the registry entry to match, or ask Demitri before deleting.
