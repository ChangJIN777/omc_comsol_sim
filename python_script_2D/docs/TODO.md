# 2D boomerang pipeline — to-do list

Status as of 2026-08-10. Ordered so each step is independently verifiable and
nothing downstream is trusted before its dependencies land.

**The template is now parametric, and that makes the face-index bug live.**

`comsol/trusty_boomerang.mph` has been rebuilt from the parameterized export
script — verified by inspecting the archive directly: it carries all 13
parameters and the `leg1/leg2/leg3` tags, and is saved solution-free at 2.7 MB
(down from 1.1 GB). So the round-2 work *is* in effect and the COMSOL Java API
in `trusty_boomerang_script.m` is validated by the fact that the rebuild
succeeded. Consequences:

1. **`_require_parameters` should now pass.** `_require_studies` will still
   raise — the template's studies are `Study_symmetric` and an unlabelled one
   (COMSOL reports `Study 1`), not `mech evenz` / `mech oddz`. Relabel in the
   GUI, or set `OMC2D_STUDY_EVENZ="Study_symmetric" OMC2D_STUDY_ODDZ="Study 1"`
   for a smoke test.
2. **The face-index bug is no longer dormant.** It was harmless only because
   the geometry never rebuilt. The geometry now rebuilds on every candidate, so
   any `(a, w, r)` that changes topology can silently reattach the Floquet pairs
   to the wrong walls — and COMSOL will solve it without complaining. **This is
   the single highest risk in the project. Do not run a sweep until section 2
   lands.**
3. **Do not vary geometry for the first smoke test.** Run at the reference
   design only, confirm the eigenfrequencies match the pre-parameterization
   model, and *then* perturb `th` — it changes no in-plane topology, so it is
   the cheapest test that face selections survive a rebuild.

Still unverified: nothing has been *solved* through the new template, so the
BC selections and eigenvalue behaviour after a rebuild are untested.

**Regression baseline.** The pre-parameterization eigenfrequencies are preserved
in [`tests/data/baseline_oddz_reference.csv`](../tests/data/baseline_oddz_reference.csv)
— 27 k-points × 10 odd-z bands, extracted from the `std1EvgFrq` evaluation group
of the solved 1.1 GB `.mph` before it was deleted. That CSV is now the **only**
surviving record; `.mph` files are git-ignored and both 1.1 GB copies are gone.
`tests/test_pipeline_2d.py` guards it (`test_baseline_fixture_is_intact`) and
pins what it says about the reference design
(`test_baseline_reproduces_the_known_gaps`).

Caveat on its strength: **odd-z only** — only `std1` was ever solved, so there
is no even-z baseline. Given that the 24 expression conversions already match
the original literals bit-exactly and the rebuilt topology is identical, this is
a confirmation rather than the primary check.

**Target frequency is wrong for this design.** The baseline's largest odd-z gap
is **12.29 → 16.64 GHz (30% normalized, centre 14.46 GHz)**. But
`targets_2d.yaml` has `target_frequency_GHz: 8.0`, and `gap_near_frequency`'s
`rel_tol=0.5` gives a [4, 12] GHz window — which *does* contain a second,
smaller gap at 5.29 → 6.42 GHz (19%). So the scorer does not report "no gap"; it
reports a plausible-looking `G_m = 0.193`, score `+0.120`, from the wrong
spectral feature, and against `min_gap: 0.20` that reads as a near-miss rather
than the 30% success it is. At `target_frequency_GHz: 14.5` the same data scores
`G_m = 0.300`, score `+0.300`. The YAML flags 8.0 as a placeholder — this is the
data to replace it from, once you decide what frequency the device needs.

---

## ✅ Done (tranche 1 — pure Python, no COMSOL)

- [x] Fixed the `r1`/`r2` swap. `r1` = **junction** fillet (inner corners at
      `w/2`), `r2` = **tip** fillet (outer corners at `hypot(r, w/2)`), per
      `addFillet()` at `buildBoomerangUnitCell.m:119-149`. Corrected in
      `geometry2d.py`, `README_template_2d.md`, and `bounds_2d.yaml`'s header
      — **but not `bounds_2d.yaml:42-43`, whose inline comments still say the
      opposite; see section 1.**
- [x] Fixed the cell in-radius: `a*sqrt(3)/4`, not `a*sqrt(3)/2` (the cell is
      a *rhombus* of side `a`, not a hexagon). Was 2x too permissive.
- [x] Fixed the web formula to use the tip **corner** distance
      `hypot(r, w/2)`, not `r + w/2` (`w` is transverse to the leg).
- [x] Added fillet-selection validity checks — the MATLAB annulus scheme
      mis-selects over part of the bounds box. (Incomplete: they ignore that
      `h_fil1` moves the junction vertices before `h_disksel2` is evaluated —
      see section 1.)
- [x] Added `clearances()` and `REFERENCE` so the formulas are testable apart
      from thresholds.
- [x] Pointed `_TEMPLATE` at `trusty_boomerang.mph` (+ `OMC2D_TEMPLATE` env).
- [x] Added `_require_parameters` / `_require_studies` guards, and
      `OMC2D_STUDY_EVENZ` / `OMC2D_STUDY_ODDZ` escape hatches — converts the
      two worst silent failures into loud ones.
- [x] Made a short eigenvalue return raise instead of NaN-padding (a NaN
      column silently zeroed the gap for an otherwise-good design).
- [x] `_study_dataset` now prefers the parametric dataset and raises instead
      of falling back to the model default (which could be the *other*
      parity's solution — a silent parity mix-up).
- [x] Corrected `materials.yaml` to `rho=3500`, `C44=578` GPa (matches
      `LoadMaterialParams.m:17,26` and the `.mph`), and flagged it dead config.
- [x] Re-anchored `run_one_2d.yaml`'s default `u` on the MATLAB reference
      design (verified: it maps to exactly 480/140/177/10/10/220 nm); the old
      midpoint of the box is infeasible (14.0 nm web).
- [x] `.gitignore`: `*.mph.lock`, `*.mph.recovery*`, `.DS_Store`. The stale
      lock file is gone.
- [x] Tests: 16 passing.
- [~] `_combined_bands` truncation for `gap_mode: "complete"`. The *idea* is
      right and the design (sort across families, compare consecutive columns)
      is right — **but the ceiling is computed with the wrong reduction and
      still admits false-positive complete gaps.** See section 1.

## ✅ Done (tranche 2 — `comsol/trusty_boomerang_script.m` parameterized)

Applies to the export script only; the `.mph` has not been regenerated.

- [x] Added `w, r, r1, r2, th` parameters (values from `test_Boomerang.m:9-16`)
      plus derived `hx, hy, selw, dsel`. All 24 geometry expressions checked to
      reproduce the original literals at the defaults.
- [x] Renamed Rectangle tags `r1,r2,r3` → `leg1,leg2,leg3` and updated the
      Compose formula. Not strictly required (COMSOL parses object names and
      numeric expressions in separate namespaces) but correct and free.
- [x] Fixed `ZsymSel`'s box to `x ∈ [-dsel, 3a/2+dsel]`,
      `y ∈ [-dsel, a√3/2+dsel]`, `z ∈ [-dsel, dsel]`. It previously selected
      nothing. `condition='allvertices'` is the right condition — do **not**
      switch to `somevertex`, which would also grab all four side walls
      (their lower edges sit at z=0) and put a Symmetry BC on the Floquet
      faces. Invariant this relies on: **`dsel < th/2`** (10 nm vs ≥75 nm
      across the bounds — safe, but state it if `dsel` ever changes).
- [x] Header block documenting what is still not done.
- [x] API forms cross-checked against existing builders: cell-array-of-char
      `set('table', {...})` (`photonic_bandStructure_code.m:46`,
      `buildBoomerangUnitCellStrip_v2.m:57`), `set('size'/'pos', {...})`
      (`buildRibUnitCell_LN.m:67-79`), `setIndex('distance','th',0)`
      (`buildLowerBoomerangUnitCell.m:76`), `set('quickz','-th/2')`
      (`buildRibUnitCell_LN.m:37`), string `set('xmin', ...)`
      (`buildLowerBoomerangUnitCell.m:103`). All have direct precedent.

---

## ✅ Done (tranche 3 — cheap correctness fixes + free/fixed variable split)

- [x] **`objective2d._combined_bands`: `nanmax` → per-k `nanmin`.** The ceiling
      is now `min over k of min(evenz[k,-1], oddz[k,-1])`. The old max-over-k
      form was the frequency below which a family *sometimes* had all its
      modes; a gap needs *always*. Fuzzed over 40 000 random spectra: old form
      admitted **419 false-positive complete gaps**, new form **0**.
      `_combined_bands` now returns `(bands, ceiling)`, and `_mechanical_gap`
      returns `(Gap, info)` carrying `truncation_ceiling_hz` and
      `n_bands_usable`, both written into the result record — so "no gap
      exists" is now distinguishable from "not enough bands were solved".
- [x] The `np.sort` NaN-placement invariant is documented in the docstring.
- [x] `bounds_2d.yaml` inline `r1`/`r2` comments corrected (they contradicted
      the header block twelve lines above).
- [x] Third fillet-selection check tightened to `w/2 + sqrt(3)*r1 < r - sw`,
      accounting for `h_fil1` displacing the junction vertices before
      `h_disksel2` is evaluated. Still assumes a 60° interior angle — confirm
      against `mphgeom` on the first rebuild.
- [x] `REQUIRED_PARAMS` extended to the full template contract
      (`a,w,r,r1,r2,th,hx,hy,selw,dsel,k,kx,ky`), so a template missing a
      derived helper fails with the driver's clear message rather than deep
      inside a COMSOL geometry rebuild.
- [x] `test_selw_matches_the_comsol_template` parses `selw` out of the `.m` and
      asserts it equals `fillet_select_halfwidth`. Still two sources of truth —
      see section 1.
- [x] `test_ceiling_uses_min_over_k_not_max` added; the older truncation test
      could not see the bug because its top bands barely varied.
- [x] **Design space reduced to 3 free variables: `a`, `w`, `r`.** `r1`, `r2`
      and `th` are now fixed at the `test_Boomerang.m` values (10, 10, 220 nm)
      in `bounds_2d.yaml:fixed`. `geometry2d.VARS`, `FIXED_VARS`,
      `optimizer.N_DIM` (asserted equal at import), `run_one_2d.yaml`'s `u`,
      and `run_one_2d.py`'s length check all follow. `phys_to_u` raises rather
      than silently dropping a conflicting fixed value; `load_bounds` validates
      that every free/fixed name is present and that none appears in both.
- [x] Tests: 23 passing.

## ✅ Done (tranche 4 — band persistence + the 2D band plotter)

- [x] **Band persistence wired in.** `evaluate_candidate` now takes
      `save_bands=True` (default) / `bands_dir`, dumps every successful
      mechanical solve to `results/bands/<id>.npz` via
      `acoustic_comsol_2d.save_bands`, and records the absolute path as
      `bands_npz`. A write failure is recorded as `bands_error` and never fails
      the candidate — a side effect must not cost a 54-eigensolve run.
      Both additions are new keys only, per the CLAUDE.md schema rule. Threaded
      through `run_one_2d.py` / `run_loop_2d.py` from `save_bands:` /
      `bands_dir:` in their configs; `$OMC2D_BANDS_DIR` overrides the default,
      mirroring `database.py`'s `$OMC2D_DB_PATH`. This also makes
      `objective2d.py`'s "re-score `symmetry` vs `complete` without re-solving"
      claim true rather than aspirational.
- [x] **`save_bands` no longer drops scalars.** It filtered on
      `isinstance(v, np.ndarray)`, so `a` — a plain Python float in the dict
      `run_mechanical_comsol_2d` returns — vanished from every `.npz`, which is
      exactly the key the plotter needs. Scalars are now stored as 0-d arrays;
      non-numeric values are still skipped so `np.load` never needs
      `allow_pickle=True`. `os.makedirs` uses `abspath` so a bare filename
      works. Returns the path. Guarded by
      `test_save_bands_keeps_scalar_lattice_constant`.
- [x] **`scripts/plot_bands_2d.py` + `configs/plot_bands_2d.yaml`.** Two input
      kinds: `baseline-csv` (the real 27×10 odd-z fixture — runs today, no
      COMSOL) and `npz`. Ticks Γ/M/K/Γ at k = 0/1/2/3 with separators at the
      segment joins; both parities in distinguishable colours; gap shading per
      `targets_2d.yaml`'s `gap_mode`; the truncation ceiling drawn as a dotted
      line in `complete` mode with an "any gap above this is an artifact"
      annotation. **`complete` mode errors out on single-parity data** instead
      of quietly reporting a symmetry gap under the wrong name. No optical
      kind. Not a copy of the 1D plotter — see its module docstring for the
      four things that differ and why each would have been wrong.
- [x] **The target-frequency trap is now a figure, not a paragraph.** The plot
      draws `target_frequency_GHz` and its ±`rel_tol` window, and shades BOTH
      the gap the scorer picks and the largest gap anywhere whenever they
      differ. On the baseline that is immediately legible: a 19.3% gap at
      5.86 GHz shaded amber inside the [4, 12] GHz window, and the real 30%
      gap at 12.29→16.64 GHz hatched green above it, unscored. That figure is
      the argument for section 6's `target_frequency_GHz` decision.
- [x] Tests: **30 passing** (+7). New: `save_bands` scalar retention, plotter
      smoke render of the baseline CSV to PNG, `complete`-mode render including
      the ceiling invariant, `complete`-on-one-parity refusal, and BZ-loop
      closure being gap-neutral and idempotent.

## ✅ Done (tranche 5 — per-parity results in the record)

- [x] **Both z-parity families are now reported, not just the winner.** The
      record gained `mech_parity` (`"evenz"` / `"oddz"`, or `"complete"`),
      `mechanical_gap_{evenz,oddz}` and
      `mechanical_center_frequency_{evenz,oddz}`. `mechanical_gap` /
      `mechanical_center_frequency` keep their existing meaning as the *scored*
      values, so this is additive per the CLAUDE.md schema rule. Both parities
      were always solved and both were already in the `.npz`; the record simply
      threw the loser away and never said which family had won — which makes a
      scored candidate impossible to judge, since `gap_mode: symmetry` is only
      acceptable *because a single parity couples to the transducer*.
- [x] Per-family gaps are populated in **`complete` mode too** (they answer
      "what would each family have given alone?", the first question when a
      complete gap comes back empty). They are computed on the *untruncated*
      family arrays, which is correct: within one family every mode below its
      top band is known, so a single-family gap needs no ceiling — the
      `_combined_bands` ceiling exists only because *stacking* leaves the
      interleaving unknown.
- [x] **Tie-break made explicit.** `objective2d.best_parity` (public) picks the
      winner and breaks exact ties toward `PARITY_ORDER[0]` = `evenz`,
      deliberately, so the label cannot flip between runs of a near-z-symmetric
      design. `scripts/plot_bands_2d.py:select_gaps` now calls the same helper
      instead of re-implementing the comparison, so a figure legend and a
      record can no longer disagree about which family won — previously they
      agreed only by coincidence (`sorted()` + a strict `>`).
- [x] `mech_parity` rides in `database.py`'s `record` JSON blob — **no new
      column, deliberately.** SQLite's JSON1 is built in here (3.53.3), so it is
      already queryable:
      `SELECT id, score, json_extract(record,'$.mech_parity') AS parity,
      json_extract(record,'$.mechanical_gap_oddz') FROM candidates
      WHERE parity='oddz' ORDER BY score DESC;`
      A real column would need an `ALTER TABLE` guarded by `PRAGMA table_info`,
      since `CREATE TABLE IF NOT EXISTS` will not migrate an existing
      `runs2d.sqlite`. Not worth it for a field that is queryable as-is;
      revisit only if parity ends up in a hot query path.
- [x] Progress output surfaces it: `run_one_2d.py`'s mechanical "done" line now
      carries `[oddz: 24.8% @ 7.99 GHz | evenz 2.3% oddz 24.8%]`, and
      `run_loop_2d.py`'s per-iteration line prints `parity=`. Both docstrings
      updated to describe the record's mechanical fields.
- [x] Tests: **35 passing** (+5). Winning-family label; losing family with no
      gap reporting `0.0` rather than NaN or a missing key; `complete` mode
      still carrying both families; the evenz tie-break pinned (plus
      `best_parity` on a single-family subset); and plotter-vs-record agreement
      checked on data where *evenz* wins, the reverse of the other fixture.

## ✅ Done (tranche 6 — "is a complete gap the overlap of the two gaps?")

- [x] **Answered and written down at the code.** `_combined_bands`'s docstring
      gained an EQUIVALENT FORMULATIONS section with the derivation, and
      `objective2d`'s module docstring points at it, because this question will
      be asked again. Short answer: stack-then-search returns a **strict
      superset** of `union over (i,j) of (evenz_gap_i ∩ oddz_gap_j)`. A window
      can be free of a family because it lies *below that family's lowest band*,
      not only because it is inside one of that family's named gaps — the
      overlap recipe cannot express those. 84% of random band pairs contain at
      least one such window, and one can be the **largest** complete gap
      (fixture in the suite: overlap 2.3%, correct answer 11.8%, a 5× miss).
- [x] The containment direction is exact and now enforced:
      `test_complete_gap_contains_every_pairwise_intersection` fuzzes random
      flat/dispersive families with unequal band counts and asserts every
      pairwise overlap below `min_k c_N` lands in exactly one combined gap.
      Verified over 30 000 pairs / 48 283 overlaps offline: 0 misses, 0
      multi-coverage.
- [x] **The qualifier is `min_k c_N`, not the ceiling** — a band straddling the
      ceiling is dropped by the `keep` mask, taking any gap that ended at its
      bottom with it. Pinned by `test_ceiling_can_hide_a_real_overlap`. One-sided
      conservatism: it can hide a real complete gap, never invent one, and the
      fix is more `neigs` rather than a different reduction. **This is a live
      argument for raising `neigs` above 10 in both studies** (already flagged
      in tranche 3 for a different reason).
- [x] **Gap EDGES added to the record** (Hz): `mechanical_gap_{lower,upper}_frequency`
      for the scored gap and `..._{lower,upper}_frequency_{evenz,oddz}` per
      family. Additive. They were algebraically recoverable as
      `centre*(1 ∓ G/2)`, but only under *this* project's normalization
      (`G = Δf/mean`, see `bandgap.Gap`) — `Δf/f_lower` is a common enough
      alternative that a reader reconstructing edges could be quietly wrong.
      Flat keys, not a nested dict, to match the existing
      `mechanical_gap_{evenz,oddz}` naming; `json_extract` reaches either, so
      nesting would only break the symmetry. Now
      `max(lower_*) .. min(upper_*)` reads the overlap straight off a record.
- [x] Tests: **40 passing** (+5). Fuzzeed containment; strict-superset via a
      family with no gap at all; an extra beating the overlap 5×; the ceiling
      hiding an overlap and reappearing when both families get another band;
      and gap edges present for the scored gap and both families, including the
      `centre*(1 ∓ G/2)` identity and all-zeros for a gapless family.
- [x] **Decided against an "intersection view" in the plotter** — see the report;
      in `complete` mode the shaded span already *is* the correct answer, and
      drawing the overlap alongside it would advertise a quantity the code
      deliberately does not score.

## ✅ Done (tranche 7 — the COMPLETE gap is now what this pipeline scores)

- [x] `configs/targets_2d.yaml`: **`gap_mode: "complete"`**. `mechanical_gap`,
      `mechanical_center_frequency`, `mechanical_gap_{lower,upper}_frequency`
      and `score` now all refer to the two-family gap. `"symmetry"` remains
      fully supported and the YAML now argues *why* it exists as an option
      (only one parity couples to the transducer) rather than presenting it as
      the sane default.
- [x] **Per-family diagnostics kept, and now load-bearing.**
      `mechanical_gap_{evenz,oddz}`, their edges, and `mech_parity` are recorded
      in complete mode too — they are the explanation for a narrow complete gap.
      `test_targets_yaml_scores_the_complete_gap` asserts they survive, so a
      future tidy-up cannot drop them.
- [x] In-code fallbacks flipped to match the shipped config
      (`objective2d.py:344`, `plot_bands_2d.py:330`), so a targets file missing
      the key cannot silently switch modes back. Pinned by a test that greps for
      the fallback string.
- [x] **The no-COMSOL demo still works.** `configs/plot_bands_2d.yaml` now pins
      `gap_mode: symmetry` explicitly, with a comment saying why: its default
      input is the odd-z-only baseline, and complete mode correctly *raises* on
      single-parity data. Verified both ways — the default invocation renders,
      and removing the pin reproduces the guard's actionable message.
- [x] Docs made coherent: `objective2d`'s module docstring, both run scripts'
      docstrings (including the expected `parity=complete` and the plateau
      warning), `plot_bands_2d.py`'s docstring, and `CLAUDE.md` — which now says
      plainly that a complete gap **is** required here, and flags that
      `../python-scripts/CLAUDE.md` says the opposite for the 1D project so the
      language is not carried over. `README.md` needed no change (it never
      described `gap_mode`).
- [x] Tests: **42 passing** (+2) — the shipped `gap_mode` matches what the code
      paths expect and takes the complete branch with diagnostics intact; and
      the default plot config both pins `symmetry` and actually renders.

---

## 1. Residual single-source-of-truth issue

- [ ] **`selw` is still stored twice** — as the `selw` parameter in the `.m`
      and as `fillet_select_halfwidth` in `bounds_2d.yaml`. A test now asserts
      they match, which catches drift but does not prevent it. Once the driver
      can write template parameters, push `selw` (and `dsel`, `hx`, `hy`) from
      Python so the YAML is authoritative.
- [ ] **The fixed values are stored twice too** — `r1/r2/th` in
      `bounds_2d.yaml:fixed` and as defaults in the `.m`. Same fix: have the
      driver push all six geometry parameters, which it already intends to do.
      Until then a template edit can silently diverge from what Python thinks
      it is solving. Consider extending the `selw` test to cover these three.
- [ ] **`load_bounds()` re-reads and re-parses the YAML on every call**, and
      both `u_to_geometry` and `check_feasibility` call it when `bounds=None`
      — two file reads per candidate. Irrelevant next to a COMSOL solve, but it
      makes feasibility Monte Carlos over the bounds box slow enough to notice
      (~40 s for 40 k samples). Cache it if that becomes annoying.

## 2. Replace the index-based BC selections — **now the top risk**

`trusty_boomerang_script.m:220-231` still pins `pbcX` to faces `[1 22]`,
`pbcY` to `[2 9]`, and both parity BCs to `[3]`. These are absolute boundary
numbers valid only for the default geometry. Growing `r` past the cell wall
splits a side face; changing `r1`/`r2` adds or removes fillet faces. COMSOL
renumbers, the BCs reattach to the wrong walls, and the solve **succeeds** —
returning a band structure for a problem nobody posed.

This was previously last-ish in the template work because the geometry never
rebuilt. Parameterization removes that protection. **Do not run a sweep until
this lands.**

- [ ] Port `bndindex.m` to Python (~30 lines: `geom.java.getVertexCoord()` for
      the in-plane test, `geom.java.getAdj(2,0)` for face→vertex adjacency,
      keep faces whose vertex set is entirely in-plane). Test that it
      reproduces MATLAB's `P.xEnd1/xEnd2/yEnd1/yEnd2/zEnd` at the reference
      design. **Highest-value single artifact in this project.**
      Call sites mirror `buildBoomerangUnitCell.m:106-112`.
- [ ] For the z=0 parity plane, `ZsymSel` is already correct and can take over
      immediately — no `bndindex` needed:

      ```matlab
      model.component('comp1').physics('smech').feature('symBCs').selection.named('geom1_ZsymSel');
      model.component('comp1').physics('smech').feature('asymBCs').selection.named('geom1_ZsymSel');
      ```

      (`buildBoomerangUnitCell.m:100` confirms the `geom1_ZsymSel` naming.)
      Do this first — it is one line per BC and removes a third of the risk.
- [ ] For the four slanted side walls, `bndindex`-style normals are
      `[sqrt(3)*a/2, -a/2, 0]` (the `pbcX` pair) and `[0, 1, 0]` (the `pbcY`
      pair), per `buildBoomerangUnitCell.m:106-110`. Note the stale comments in
      `runBands_2D.m:212,225` ("yz planes at x = ±a/2") describe the 1D beam,
      not this cell — trust the `bndindex` calls.
- [ ] **Hard invariant for this section:** at the reference design the
      resolved indices must equal `[1 22]`, `[2 9]`, `[3]`. If they don't, the
      port is wrong, not the old model.

## 3. Regenerate the `.mph` from the parameterized `.m`

The parameterization is inert until this happens.

- [ ] Comment out the solve tail (`model.study('std1').runNoGen;` at :330 and
      the two `evaluationGroup(...).run` calls at :339 and :350), run
      `model = trusty_boomerang_script(); mphsave(model, 'trusty_boomerang.mph')`.
      Note the script as written solves **only `std1` (odd-z)** — the even-z
      study is never computed — so leaving the tail in gives a half-solved,
      still-enormous file.
- [ ] Relabel the two studies to `mech evenz` / `mech oddz` (currently
      `Study_symmetric` and, unlabelled, `Study 1`). Parity mapping is already
      correct: `symBCs` active → z-even, `asymBCs` active → z-odd. Until then,
      `OMC2D_STUDY_EVENZ="Study_symmetric" OMC2D_STUDY_ODDZ="Study 1"` works
      as a stopgap.
- [x] Saved **solution-free** — 1.1 GB → 2.7 MB, 0 stored solutions. Note it is
      **not mesh-free**: `mesh1.mphbin` is 2.07 MB of the 2.7 MB (77%) and
      `nodeType="meshed"`. That is arguably better for the driver (no remesh per
      process), but if you ever want a committable template the mesh is what
      to drop next.
- [x] Add `hx, hy, selw, dsel` to `REQUIRED_PARAMS`. Done, and verified: the
      `.mph`'s parameter set and `REQUIRED_PARAMS` now match exactly in both
      directions, so `_require_parameters` passes.
- [x] The 3-arg `model.param.set(name, expr, descr)` overload **works** — all 13
      descriptions are stored as `descr=` attributes in the rebuilt `.mph`. No
      fallback needed. (Also retired by the successful rebuild: the
      cell-array-of-char polygon `table` form and string-valued `DiskSelection`
      properties. All 42 geometry features report `BUILT`, 0 errors, 0
      warnings, and topology is identical to the pre-parameterization model at
      24 faces / 66 edges / 44 vertices.)
- [ ] ~~Low-risk but unverified~~ — superseded, see above. Retained for
      context: if a future rebuild ever errors on the 3-arg form, split into
      `model.param.set(name, expr); model.param.descr(name, descr);` — the
      form COMSOL's own GUI export emits. `'480[nm]'` in place of
      `'4.8e-07[m]'` is fine; COMSOL parses the unit bracket.
- [ ] **Hard invariant for this section: eigenfrequencies at the reference
      design must not move.** That is the only check that parameterization
      preserved the geometry. Two solves.

## 4. Decide on `geomRep('cadps')` — now affects *which designs build*, not just licensing

`trusty_boomerang_script.m:104-105` selects the Parasolid kernel and
`designBooleans(false)`. Both need the CAD Import Module / Design Module.
Previously this was a headless-licence question. Now that the geometry
rebuilds per candidate it is also a **robustness** question: Parasolid and the
COMSOL kernel differ in fillet behaviour on tight corners, so the choice
changes which `(w, r, r1, r2)` combinations produce a buildable hole — i.e. it
moves the feasible region that section 6 is trying to map.

- [ ] Revert to `geomRep('comsol')` unless Parasolid is needed.
      `buildBoomerangUnitCell.m` never sets it.
- [ ] Either way, gate it at driver startup: set the reference parameters and
      call `(model/'geometries'/'Unit Cell').java.run()` once, so a licence or
      kernel failure surfaces immediately rather than mid-sweep.

## 5. Rewrite the k-sweep — 46× fewer eigensolves

`acoustic_comsol_2d.py:246-247` writes `kx`/`ky` as *literals* per point, which
destroys the `if(k<1,...)` expressions the `.m` now explicitly warns about at
`:81-84`; the study's own 27-point `k` sweep then solves the same wavevector 27
times and `evaluate("freq")` returns the fundamental repeated. `comsol_client`
caches the model process-wide, so the damage persists across every candidate.
The module docstring at `:205-214` already describes the fix.

| approach | eigensolves / candidate |
|---|---|
| current code as written | 2,484 (degenerate garbage) |
| per-point loop, if `kx`/`ky` were left alone | 92 |
| **COMSOL sweep at `kpts=9`** | **54** |

- [ ] Set only `k`; let the Parametric step sweep once per parity; read back
      `[n_k, n_bands]` via `model.outer(dset)` + per-outer `model.evaluate`.
      Do **not** rely on COMSOL's flat outer-major ordering — undocumented and
      version-dependent.
- [ ] Also write the mirrored `model.batch('pbatch')` parameter list
      (`trusty_boomerang_script.m:302-311` holds a frozen 27-entry copy that
      desyncs silently if only the study step is updated).
- [ ] Default `n_per_segment` 15 → **9** to match `test_Boomerang.m:22`.
- [ ] Copy row 0 for the closing Γ point instead of re-solving it (matches
      `runBands_2D.m:556`).
- [ ] Call `comsol_client.release_model()` / restart before the first run.

## 6. Resolve `min_web`, then retune `bounds_2d.yaml`

`min_web` is **calibrated, not derived**: 15 nm so the reference design
(17.5 nm corner clearance) passes. Two unresolved readings — either 30 nm is
the wrong threshold for a *corner* measure (the corner is rounded by `r2`, so
real material there is wider than the sharp-corner estimate), or the reference
really is that tight.

Monte Carlo over 40k uniform samples of the **3-D** `(a, w, r)` box: only
**36.0%** is feasible (51.0% rejected on `web`, 41.0% on fillet selection —
the fillet share rose from 36.2% when the post-fillet junction displacement
was added to the check). Dropping from 6 free variables to 3 did *not* help
the feasible fraction, because both dominant constraints couple only `a`, `w`
and `r` — the three that remain free. The optimizer would still spend nearly
two-thirds of its budget on rejections.

- [ ] Measure the true minimum solid width at the reference design with
      `mphgeom` once section 3 lands, and replace the calibrated 15 nm.
      Update `test_reference_design_clearances_are_pinned` only if the geometry
      *model* changed — not to make a failing test pass.
- [ ] Tighten `r`'s upper bound relative to `a`, or make `r` a ratio variable
      `r/a` so the feasible region is closer to axis-aligned in `u`-space.
- [ ] Widen the fillet-selection-valid region at the source: give the builder a
      geometry-derived annulus half-width, e.g.
      `selw = min(r - w/2, hypot(r, w/2) - r)/2` instead of the fixed 25 nm,
      then drop the corresponding checks from `check_feasibility` *and* the
      duplicated constant in `bounds_2d.yaml`. This departs from
      `buildBoomerangUnitCell.m` — worth it, but note it in the header so the
      divergence from the MATLAB reference is deliberate and visible.

## 7. Persist what the solves cost

Band data is now saved (tranche 4). What remains:

- [ ] Nothing prunes `results/bands/`. One `.npz` per candidate at
      2 parities × 28 k × 10 bands ≈ 5 KB, so a 60-candidate loop is ~300 KB —
      not urgent, but there is no retention policy and `_hash_u` reuses the id,
      so a re-run silently overwrites the previous bands for that `u`. Decide
      whether that is wanted before the provenance change below lands (it will
      change the ids and orphan every existing file).
- [ ] Stamp config provenance into the record **and** into `_hash_u`:
      `gap_mode`, `n_bands`, `min_gap`, `target_frequency_GHz`, bounds file.
      `database.py` is `INSERT OR REPLACE` on that id, so re-running a `u`
      after flipping `gap_mode` silently overwrites, and
      `load_completed_results()` then feeds incommensurate scores to the
      optimizer. Back-compatible addition, per the CLAUDE.md rule.
- [ ] Expose `k_points_per_segment` in `targets_2d.yaml` and thread it through
      `objective2d.py:evaluate_candidate` (currently pinned at the default).
- [ ] Write `None`, not `f_o_target`, for `optical_center_frequency` when
      `require_opt=False` — a 193.4 THz value next to `optical_gap=0.0` is
      indistinguishable from a real optical solve that found no gap.

## 8. Wire `materials.yaml` in, or delete it

Nothing in `src/` reads it; the material comes from the `.mph`'s `diamond`
node. Editing the YAML changes nothing about a solve. It is now *correct* dead
config rather than *wrong* dead config, which is an improvement but not a fix.

- [ ] Preferred: add `C11,C12,C44,rho,rxtal` as COMSOL parameters and push
      them from the YAML, computing the rotation via a rotated coordinate
      system instead of the pre-baked `DVo` at
      `trusty_boomerang_script.m:239-242`.
- [ ] Otherwise: add a test that recomputes the 45° Bond rotation from
      `materials.yaml` and asserts it matches the six independent `DVo`
      entries (verified by hand: 1178.5 / 22.5 / 125 / 1076 / 578 / 475.5 GPa
      — an exact match only for `C44 = 578`, which is why 577 was wrong).

## 9. Validation — the gate on trusting anything

**The first end-to-end run must solve BOTH studies, not one.** Now that
`gap_mode` is `"complete"`, a single-parity solve cannot produce a score at all.
And there is no even-z band structure anywhere in this repo: `std1` (odd-z) was
the only study ever solved, and `tests/data/baseline_oddz_reference.csv` is
odd-z only. **The reference design's complete gap is therefore unknown**, and no
config or doc should be read as claiming otherwise. The `.m`'s solve tail
(`model.study('std1').runNoGen`) also only computes odd-z — see section 3.

- [ ] Run one candidate at the reference design with `kpts=9`, **both
      parities**, and diff the odd-z bands against
      `tests/data/baseline_oddz_reference.csv` (27×10, and readable without
      COMSOL). **This single comparison validates sections 2, 3 and 5
      simultaneously. Nothing downstream should be trusted until it passes.**
- [ ] Record the reference design's first-ever **even-z** spectrum and its
      **complete** gap. Two numbers to write down from that record:
      `truncation_ceiling_hz` and `n_bands_usable` — they decide whether
      `n_bands`/`neigs` = 10 is enough (section 11), and there is currently no
      way to know from one parity.
- [ ] Then perturb one variable at a time (start with `th`, which changes no
      in-plane topology) and confirm the resolved face indices track. This is
      the cheapest test that section 2 actually worked.
- [ ] Only then let the optimizer vary anything.

## 10. Consequences of scoring the COMPLETE gap (new, from tranche 7)

Ordered by how much they can distort a run.

- [ ] **The gapless plateau — the real optimizer risk, and `min_gap` does not
      fix it.** Every candidate with no complete gap in the target window gets
      `G_m = 0` *and* `f_center = 0`; the frequency penalty
      `lambda_mech_freq*((f_c - f_t)/f_t)^2` then saturates at exactly 1.0, so
      all of them score an identical **−1.40**. Complete mode makes that
      plateau much larger than symmetry mode did. Note what is *not* wrong: the
      quadratic hinge gives a *steeper* `dS/dG` below threshold (5.0 at `G=0`
      for `min_gap=0.20`, vs 1.0 above threshold), so candidates that *do* have
      a gap are well separated. The problems are the flat plateau and a **+1.05
      score cliff** between "no gap" and "any gap at all", which is bad for any
      surrogate that assumes smoothness (Optuna/TPE copes; a GP would not).
      Fix by giving gapless candidates a weak ordering signal, e.g. score the
      best single-family gap at a small weight, or use the largest complete gap
      *anywhere* so `f_center` is not identically 0. This is a scoring-design
      change; do it before switching the optimizer backend to anything
      model-based.
- [ ] **Decide `min_gap` deliberately.** 0.20 was inherited from symmetry mode,
      not chosen for complete mode, and is flagged as such in
      `configs/targets_2d.yaml`. **Recommendation: 0.10**, applied only once
      someone owns the decision — see the report and the YAML commentary. It is
      a low-stakes knob for optimizer behaviour (it shifts the plateau and
      changes the sub-threshold slope) and a high-stakes one for
      interpretability, so set it to the gap you would actually accept.
- [ ] **Re-target `target_frequency_GHz` using two-parity data.** The 8 GHz
      placeholder was already suspect from the odd-z baseline (whose real 30%
      gap sits at 14.46 GHz, outside the [4, 12] GHz window — see
      `results/figures/bands_2d.png`). A complete gap will sit at a different
      frequency again, so re-derive it from the first two-parity solve rather
      than from the odd-z figure.

## 11. `neigs` / `n_bands` — verify, don't pre-emptively raise

Complete mode is truncation-sensitive in a way symmetry mode is not: the usable
range ends at `min_k c_N`, which sits below the ceiling, and a band straddling
the ceiling silently removes any gap whose upper edge was that band's bottom
(`test_ceiling_can_hide_a_real_overlap`). One-sided: it hides real gaps, never
invents them.

**Recommendation: keep 10 for now, and verify empirically rather than guessing.**
The odd-z baseline's 10th band bottoms at 30.03 GHz, ~2.5× above the [4, 12] GHz
scoring window, so the target region has ample headroom from the one parity we
can see. Raising `neigs` on that basis alone would buy unmeasured safety at a
real cost.

- [ ] After the first two-parity solve, read `truncation_ceiling_hz` and
      `n_bands_usable` out of the record (stored for exactly this purpose). If
      the ceiling is within ~1.5× of the scoring window's upper edge, raise.
- [ ] If raising: **change all four together or the change does nothing.**
      `n_bands` in `configs/targets_2d.yaml`, plus `neigs` on `std_eigv` and
      `std_eigv1` (the two study steps) and on `solv_eigv` (the solver feature)
      in `comsol/trusty_boomerang_script.m` — then rebuild the `.mph`.
      `n_bands` > `neigs` makes `run_mechanical_comsol_2d` raise; `neigs` >
      `n_bands` computes modes that are then discarded and leaves the ceiling
      where it was.
- [ ] Cost if it comes to it: 10 → 16 is roughly +60% on the dominant
      eigensolve cost, i.e. ~54 → ~86 eigensolve-equivalents per candidate at
      `kpts=9`. Worth it to stop hiding gaps; not worth it speculatively.

## 12. Cleanup, last

- [ ] Rewrite `README_template_2d.md` from the *rebuilt* model rather than from
      intent. It has been patched (r1/r2, C44=578/rho=3500, "keep the `k`
      parameter") but still describes a template that does not exist on disk.
- [ ] Fix `runBands_2D.m:73`'s `ky_norm` third-segment formula:
      `(3*ki - kpts)*4/(3*kpts)` gives `32/3` at `ki = 3*kpts` instead of `0`.
      Affects only MATLAB's bookkeeping array, not its solve — the COMSOL
      expression at `:64` is correct and Python ported the correct one. But the
      two `ky` arrays will disagree on the K→Γ leg if compared.
- [ ] `_kxky_square` in `acoustic_comsol_2d.py` is unreachable (nothing
      overrides `unitcell="hexagonal"`). Keep it for MATLAB parity, but the
      `bz_path` docstring conflates Γ→M with Γ→X.
- [ ] `trusty_boomerang_script.m:53` hard-codes `model.modelPath('/Users/...')`
      — harmless for MPh, but it makes the script non-portable if anyone else
      replays it.
