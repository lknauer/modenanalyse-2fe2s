# Changelog

All notable changes to `modenanalyse_2fe2s` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [1.1.0] — 2026-06-01

Renames the secondary-structure abbreviation **SS → SSE** project-wide and
ships the v1.0.5 secondary-structure work as a stable release.

### Changed (breaking: output sheet and file names)

- **`SS` → `SSE` throughout.** The abbreviation for *secondary structure
  element* was renamed across code, configuration, output and docs. This
  removes the clash with disulfide (S–S / cystine) — a real ambiguity in
  sulfur-rich iron–sulfur systems — and matches the conventional term
  *secondary structure element*. Concretely:
  - **Excel sheets**: `SS_amplitude_mean` → `SSE_amplitude_mean`,
    `SS_com_amplitude` → `SSE_com_amplitude`, `SS_UMAP_clusters` →
    `SSE_UMAP_clusters`, `SS_UMAP_profile` → `SSE_UMAP_profile`, and the
    other `SS_*` metric sheets likewise.
  - **Output filename**: `*_analysis_SS.xlsx` → `*_analysis_SSE.xlsx`
    (and `*_analysis_SS_interp*` → `*_analysis_SSE_interp*`).
  - **Public API**: `analyze_ss_element` → `analyze_sse_element`,
    `analyze_all_ss` → `analyze_all_sse`, `build_ss_center_map` →
    `build_sse_center_map`, `build_ss_ca_center_map` →
    `build_sse_ca_center_map`, `compute_ss_umap_cluster` →
    `compute_sse_umap_cluster`, `detect_ss_dssp` / `detect_ss_phipsi` →
    `detect_sse_*`, `export_ss_excel` / `export_ss_interp_excel` →
    `export_sse_*`, and the constant `SS_UMAP_METRICS` → `SSE_UMAP_METRICS`.
  - **Per-mode result key**: `result["ss"]` → `result["sse"]`.
  - **Documentation** (Manual, Anleitung, Supplement, tutorials) updated to
    match; the in-repo doc-reference tests enforce code/doc consistency.
- **Configuration keys** `analyze_ss` → `analyze_sse` and `ss_chain` →
  `sse_chain`. The legacy names remain accepted in TOML with a deprecation
  warning, so existing configurations keep working unchanged.

### Migration

Re-run affected systems to regenerate `*_analysis_SSE.xlsx`. Downstream
code that reads the old `SS_*` sheet names or the `*_analysis_SS.xlsx`
filename must be updated to the `SSE_*` names. TOML configs need no change
(legacy keys are aliased), though renaming `analyze_ss` / `ss_chain` to the
new keys is recommended.

## [1.0.5] — 2026-06-01

Physics-correctness release for the secondary-structure (SS) analysis,
triggered by a code review of the per-SS-element descriptor computation
(`analyze_ss_element`) and the SS-UMAP clustering. One confirmed bug and
several physical/statistical inaccuracies of the same root classes —
silently dropped features, a mislabelled geometric quantity, unweighted
centroids used as centres of mass, and non-propagated uncertainties —
are fixed here, with eight new validation tests.

**IMPORTANT — outputs change; re-run required.** The Excel **sheet and
column names are unchanged** (downstream tooling keeps working), but the
numeric contents of `SS_com_amplitude`, `SS_tilting_angle`,
`SS_internal_amplitude`, the `SS_UMAP_*` sheets, and the cluster `kern`
COM differ. `tilting_angle` in particular changed meaning (see below).
Existing `_analysis_SS.xlsx` files must be regenerated.

### Fixed

- **SS-UMAP dropped all lateral/bending motion (silent feature bug).**
  `embedding.compute_ss_umap_cluster` listed two features as
  `bending_std`/`bending_mean`, which `analyze_ss_element` never
  produces (it returns `lateral_std`/`lateral_amplitude`). Pulled via
  `.get(metric, 0.)`, both columns were identically zero for every mode
  and every element, so lateral/bending motion never entered the
  embedding or the cluster Z-score / Fisher-F profiles (the 36
  `bending_*` features had Fisher-F ≡ 0). Now uses the real keys. A new
  regression test parses the feature list and asserts every metric is a
  key actually produced by `analyze_ss_element`.

- **`tilting_angle` was not a tilt.** It was the polar angle of the
  centroid *translation* relative to the element axis — a property of
  the net translation, not a rotation — which sat near a constant
  ~57.3° (the mean angle of an isotropic vector to a fixed axis) for
  delocalised modes and carried essentially no structural signal. It is
  now a genuine rigid-body tilt: the rocking rotation perpendicular to
  the helix axis, obtained from the mass-weighted least-squares
  infinitesimal rotation `ω` (`J ω = b`), reported as a small-angle
  amplitude in degrees.

- **`com_amplitude` was not mass-weighted.** Both the per-SS-element COM
  and the cluster-core COM (`kern`) used the unweighted atom-mean of the
  displacement vectors. Because C/N/O/S (and Fe/S in the core) differ in
  mass, the unweighted centroid is not the centre of mass. Both are now
  true mass-weighted COM displacements (`Σ mᵢ uᵢ / Σ mᵢ`).

- **`internal_amplitude` is now a rigid-body residual.** It was the
  lateral part of the centroid-relative motion (translation removed
  only). It is now the TLS-style residual after removing both
  translation *and* the rigid rotation `ω`, i.e. a true non-rigid
  internal-strain measure.

- **Principal axis is now backbone-based.** The SVD axis was computed
  over all heavy atoms (including side chains), which can tilt the axis
  away from the helix direction for short elements. A new
  `geometry.build_ss_ca_center_map` provides per-element Cα centres; the
  axis now prefers the backbone trace and falls back to the all-atom
  axis when fewer than two Cα are resolved (backward-compatible default).

- **SS uncertainties are now honest first-order estimates.** They were
  all set to a single `s_amp = u_rms·σ_eigvec·√n` with arbitrary ×1.4 /
  ×1.5 factors, and a displacement-unit σ was assigned to an angle in
  degrees. The `√n` scaling was that of a *sum*, not a mean. The SS
  sigmas now use the standard error of a mean (`σ_eigvec/√n`), drop the
  magic factors, and propagate the tilt uncertainty as a small angle
  over the radius of gyration (`σ/Rg`).

### Changed

- **SS-UMAP feature set reduced and decorrelated (changes clustering output).**
  The per-element clustering features were cut from nine to five. Dropped:
  the two overall-magnitude features (`amplitude_mean`, `amplitude_max`),
  which are collinear with the directional components and over-weighted the
  "overall amplitude" axis under z-scoring, and the two second-moment
  features (`lateral_std`, `stretching`), which largely duplicate their
  first moments. Retained: the orthogonal rigid-body partition
  (`com_amplitude`, `tilting_angle`, `internal_amplitude`) plus the
  axial/lateral directional split. All nine descriptors are still written
  to the per-metric SS sheets; only the clustering input changed. This
  affects `SS_UMAP_clusters` and `SS_UMAP_profile`.
- **Static cleanup:** removed unused imports across `config`, `core`,
  `geometry`, `logio`, `pcet_et`, `reorganization` and `embedding` (no
  behavioural change). The `runner` orchestrator was left untouched by
  design; the deliberate `hdbscan` availability probe is retained.

### Notes

The SS-UMAP metric list was the one feature set that hard-coded its names
locally instead of sharing the producer's constant (contrast
`core.SCORE_KEYS`), which is how it drifted. It is now promoted to a
shared single source of truth, `core.SS_UMAP_METRICS`, imported by
`embedding`; a regression test additionally asserts every metric is a key
produced by `analyze_ss_element`. Defaults are backward-compatible
(`axis_centers=None`, `ss_axis_center_map=None` reproduce the old
all-atom axis), so the pre-existing test suite is unaffected. The
reorganization-energy, PCET/ET, Debye-Waller / Fe-pDOS thermodynamics
and SCSD physics were reviewed alongside this release and are unchanged;
the deprecated NIS-spectrum parameters (`nis_*`) remain accepted but
unused, as before.

## [1.0.4] — 2026-05-12 (revised 2026-05-13)

Audit-driven bugfix release. A defense-in-depth audit triggered by a
user observation ("the NH-stretch sigmas seem too large relative to the
values") found four numerical bugs of the same root class — silent
inconsistent scaling of uncertainties — plus four additional defensive
gaps. All eight findings are fixed in this release, with twelve new
regression tests.

The 2026-05-13 revision adds four further reporting-only fixes
(FUND 14–17) discovered during analysis of a real-world Apd1 Rieske
batch run. Test count: 150 → 154.

### Post-release audit patches (2026-05-13)

A separate audit of a real-world [2Fe-2S]-Cys$_2$His$_2$ batch run
(Rieske-type system in both protonation states, 4 systems × 4
frequency windows, ~5900 modes each) identified four additional
issues. None affect numerical output; all are warning- or
reporting-only. Detailed in Supplement.pdf §FUND 14–17.

- **FUND 14** (`geometry._add_his_hn_info`). The
  `"His_HN NICHT erkannt"` warning fired twice per deprotonated His
  ligand because the emission was inside the outer
  `for use_pdb_only in (True, False)` loop. Hoisted out: now exactly
  one warning per ligand.
- **FUND 15** (`runner._run_analysis_single`). The PCET-disabled
  message stated `"no His ligands at cluster"` for systems with
  all-deprotonated histidines (e.g.\ Rieske-type Cys₂His₂ in the
  deprotonated state) even though `Lambda_FeN > 0` and `Fe_N_His`
  sheets confirmed His ligands were present. Root cause:
  `pcet_info.n_his` counts only protonated His. Fix: count total
  His via `coord_info.ligands` and distinguish "no His at all"
  from "His present but all deprotonated".
- **FUND 16** (`export.export_interpolated_excel`). The warning
  `"interp_boundary_mode='context' but no context modes available"`
  fired even for full-spectrum runs (`freq_min`/`freq_max` both
  `None`), where context-loading is correctly skipped and no
  anchoring is needed. Fix: only warn when the user actually set
  `freq_min`, `freq_max`, or `freq_windows`.
- **FUND 17** (`logio.RunLog.write_befund`). REPORT.txt showed
  `freq_filter: - - - cm-1` when only `freq_windows` was set,
  hiding the active filter from the audit trail. Fix: emit an
  additional `freq_windows: [...]` line whenever
  `cfg.freq_windows` is non-empty.

Four new regression tests added; total test count now 154.

### Fixed — numerical (high impact)

- **NH-stretch sigma scaling** (`core.analyze_his_hn`). The eigenvector
  rows passed in were thermally scaled by `u_rms`, but the sigma was
  computed as `sqrt(2) * cfg.sigma_eigvec` without that factor. For
  typical low-frequency modes (`u_rms ≈ 0.05 Å`), this made
  `s_hn_stretch` about 20× too large, hiding genuine PCET signals
  behind apparent noise. **The values were correct; only the sigmas
  were inflated.** The function now takes a `u_rms` argument and uses
  `cfg.sigma_eigvec * u_rms * sqrt(2)`, matching `analyze_fe_ligand`.
- **SS-element sigma hardcoded** (`core.analyze_ss_element`). Used a
  hardcoded `1e-4` regardless of `cfg.sigma_eigvec` (default `5e-4`),
  silently overriding the user's TOML configuration for SS-element
  uncertainties only. Out-of-the-box this made `s_amplitude_mean` etc.
  about 5× smaller than the rest of the pipeline. Now uses
  `cfg.sigma_eigvec` consistently. Threaded through `analyze_all_ss`
  and the runner call site.
- **SCSD sigmas hardcoded** (`core.compute_scsd_for_mode_full`). The
  geometry and eigenvector uncertainties were hardcoded as `5e-7` and
  `5e-6`, about three orders of magnitude smaller than `cfg.sigma_coord`
  (`1e-3`) and `cfg.sigma_eigvec` (`5e-4`). SCSD uncertainties now
  respond to TOML configuration.

### Fixed — silent data loss (defense-in-depth)

- **`analyze_fe_ligand` silent zero fallback** (`core.py`). If a Fe or
  ligand center fails the `c2l` lookup, the function previously
  returned `_zero_lig()` silently — the same root cause class as the
  v1.0.2 `_build_group_map` bug. Now emits a deduplicated
  `UserWarning` per `(label, failure_reason)` pair, routed through
  `_WARNED_EMPTY_LIG` (analogous to `_WARNED_EMPTY_GROUP`).
- **`analyze_his_hn` silent skip** (`core.py`). Distinguishes between
  legitimate deprotonation (`his_protonated=False` or `h_center=None`
  → silent skip, no warning) and lookup failure for protonated
  histidines (now emits a deduplicated `UserWarning`). Routed through
  `_WARNED_HN_SKIP`.
- **`_ws_fe_bindung` all-zero-row detector** (`export.py`). Mirrors the
  v1.0.2 detector in `_ws_gruppen`. If a ligand has zero stretch/bend
  across every mode, warn at the export layer — downstream check that
  catches any future silent failure that slips past
  `analyze_fe_ligand`.
- **`_ws_his_hn` all-zero-row detector** (`export.py`). Same idea, but
  only warns for ligands flagged `his_protonated=True` (deprotonated
  His legitimately produces an empty H-N entry).

### Fixed — index drift (defensive)

- **Torsion loop in `analyze_mode`** (`core.py`). The original loop
  used `ai = enumerate(gctr)` to address both `gctr` AND the filtered
  sub-eigenvector `evg_g`. If `c2l` ever missed any center earlier in
  `gctr` (`evg_g` is built via `_evg_sub` and drops missing centers),
  every subsequent `evg_g[ai]` access referenced the wrong residue's
  row. In normal operation `c2l` and `idx_map` cover the same atoms
  and the bug never triggered, but if they diverged the torsion sum
  would silently mix rows across residues. The loop now uses an
  explicit `ai_local` counter that advances only for centers actually
  consumed from `evg_g`.

### Added — tests

`tests/test_v104_audit_fixes.py` (12 new tests, all passing):

- `test_his_hn_sigma_scales_with_u_rms` — sigma is exactly proportional
  to `u_rms`.
- `test_his_hn_value_and_sigma_have_same_scaling` — significance ratio
  (value/sigma) is invariant under `u_rms` rescaling.
- `test_ss_element_sigma_responds_to_sigma_eigvec` — doubling
  `sigma_eigvec` doubles the reported `s_amplitude_mean`.
- `test_scsd_sigmas_use_cfg` — SCSD reference sigmas double when
  `sigma_coord` is doubled.
- `test_analyze_fe_ligand_warns_on_missing_fe_center` — UserWarning is
  emitted for `c2l` lookup failures.
- `test_analyze_fe_ligand_warning_deduplicated` — five identical
  failures produce exactly one warning per session.
- `test_analyze_his_hn_warns_on_protonated_lookup_failure` — protonated
  His with missing H/N emits a warning.
- `test_analyze_his_hn_does_not_warn_for_deprotonated` — deprot His
  emits no warning (legitimate silence).
- `test_ws_fe_bindung_warns_on_all_zero_ligand` — Excel-layer check
  fires when stretch/bend identically zero.
- `test_ws_fe_bindung_no_warning_when_data_present` — no false
  positive when any mode has nonzero data.
- `test_ws_his_hn_warns_only_for_protonated_zero` — His_HN detector
  ignores deprot ligands.
- `test_torsion_loop_robust_to_c2l_gaps` — index drift fix in
  isolation.

Total test count: 118 unit + 4 E2E smoke = 122 tests, all green
(previously 110).

### Practical impact on existing data

If you have v1.0.2 or v1.0.3 results in production (e.g. for a thesis
or paper), the underlying numerical **values** (`hn_stretch`,
`amplitude_mean`, `SCSD_*`, etc.) are correct. What changes in v1.0.4
is the reported **uncertainty** of those values:

- `s_hn_stretch` is ~20× smaller (was inflated)
- `s_amplitude_mean` for SS-elements is ~5× larger
- SCSD sigmas are ~10³× larger

The classification of which features cross a given Wert/Sigma > 1, > 3
threshold may therefore shift. Re-running with v1.0.4 is recommended
before drawing significance conclusions from previous outputs.

### Why these specific bugs

All eight findings share a root cause: **silent inconsistent scaling**
of an uncertainty against the value it qualifies. The v1.0.2 Apd1 bug
was the canonical example (groups sheet rows silently zero). The v1.0.4
audit systematically searched for the same pattern in: every
`s_*`-variable assignment, every `.get(..., 0.)` chained fallback,
every hardcoded constant that should have been a cfg parameter, and
every loop that indexed into a filtered sub-array using the unfiltered
index. The audit found exactly four additional sites; this release
fixes all of them and adds detection (warnings + tests) so the same
class cannot recur silently.

### Fixed — deep-sweep findings (FUND 9, 10, 11)

A second-round audit covering SCSD pipeline, PCET acceptor search,
multi-cluster aggregation, Kabsch alignment, eigenvector indexing,
B-factor accumulation, interpolation, and window-boundary handling
found three additional issues:

- **`build_pcet_info` silent skip on missing H index** (`pcet_et.py`).
  Same bug class as the v1.0.4 `analyze_his_hn` fix: a protonated His
  whose H center is not in `idx_map_h` was skipped without warning.
  Now emits a UserWarning so missing PCET HA channels become visible.
- **PDB-Gaussian matching ambiguity warning** (`geometry.py`). The
  matcher reported `n_ambiguous` (PDB atoms with multiple Gaussian
  candidates within tolerance) only as an info line. If >5% of atoms
  are ambiguous, the assignments are unstable and `cfg.coord_match_tol`
  should be tightened — this now triggers a UserWarning.
- **Window boundaries half-open** (`runner.py`). Multi-window mode
  used `[lo, hi]` closed-closed intervals, so a mode with frequency
  exactly equal to a boundary (e.g. `freq = 100.0` for windows
  `(0,100)` and `(100,300)`) was counted in BOTH windows. This caused
  silent double-counting of B-factor and Lambda_cumulative
  contributions for boundary modes. Intervals are now `[lo, hi)`
  half-open except for the very last window, which stays `[lo, hi]`
  so its upper edge mode is not lost. With realistic float
  frequencies the bug rarely triggered in practice, but the fix
  guarantees binning correctness.

Three more regression tests in `tests/test_v104_audit_fixes.py`
(`test_build_pcet_info_warns_on_missing_h_index`,
`test_pdb_matching_warns_on_high_ambiguity` /
`_no_warning_for_low_ambiguity`,
`test_window_boundaries_are_half_open`).
Total v1.0.4 test count: **16 audit tests** (was 12), **122 unit
tests** (was 118), **126 with E2E**.

### Verified clean by deep sweep

The following pipeline regions were systematically checked and are
free of the bug classes the v1.0.4 audit targeted:

- SCSD canonical D2h model (`_SCSD_MODEL_COORDS`): hardcoded
  `Fe-Fe=2.73 Å, Fe-S=2.20 Å`, but used as a *reference* for the
  difference method `SCSD_d<Irr> = SCSD_dist - SCSD_ref`, so the
  difference is robust against real cluster geometry.
- Multi-cluster pipeline (`runner._run_analysis_single` invoked per
  cluster with `dataclasses.replace(cfg, ...)`): no state contamination
  between clusters; `reset_warning_state()` resets all three
  deduplication sets per cluster.
- `get_eigvec_orca` / `parseresult_to_atoms`: heavy-only vs all-H
  filtering preserves center-IDs as 1-based original-index, so
  `ev_3d[center - 1]` always points to the correct eigenvector row
  regardless of which filter the atoms list was built with.
- Kabsch alignment (`geometry.kabsch_align`): standard SVD with
  cluster-anchor-Atom-ordering, robust to multi-cluster systems via
  best-fit Fe-pair selection.
- Mode numbering (1-based `bi.mode_nums[col]` vs 0-based `col`):
  consistently respected at every callsite.
- Debye-Waller accumulation: `_evg` is the u_rms-scaled eigenvector;
  `b_accum += sum(_evg**2)` correctly gives `u_rms² * |e|²` per atom,
  matching the textbook `B = (8π²/3) Σ <u²>`.
- `np.interp` boundary modes (`"context"`, `"zero"`, `"nearest"`):
  configurable via TOML, sortiert-monoton-Garantie auf der x-Achse,
  defensiv gegen doppelte Frequenzen.

### Fixed — FUND 12: interpolation context guarantees a boundary anchor

User-reported follow-up to FUND 11: when `freq_max` is set and the
next mode above it is farther than `interp_context_cm1` (default
5 cm⁻¹), no context mode was loaded and the interpolated pDOS lost
its right-edge anchor. Symptom: the interpolated curves drop sharply
to zero just below `freq_max` instead of decaying naturally toward
the next real mode.

Fix (`runner.py`): if no context modes exist in the
`[freq_max, freq_max + interp_context_cm1]` window, fall back to
loading exactly **one** mode -- the single closest mode above
`freq_max`. This guarantees `np.interp` has a right-edge anchor
regardless of mode density. The same fallback is applied symmetrically
at the lower edge (one mode below `freq_min` if the buffer zone is
empty). The run-log records when this fallback is used so users can
verify the anchor frequency.

### Fixed — FUND 13: synthetic zero anchor when spectrum ends inside the range

Follow-up to FUND 12: even the "one mode above" fallback can be
empty -- when the DFT spectrum genuinely ends inside the analysis
range (e.g.\ `freq_max=800` but the highest real mode is at 750).
Pre-v1.0.4 fell through to `np.interp`'s `right=0.0`, which produces
a sharp step at the upper edge rather than a smooth decay.

Fix (`runner._make_synthetic_zero_mode`): when no real mode exists
above `freq_max`, insert a **synthetic null-mode** at
`freq_max + interp_context_cm1` with all observable fields = 0. This
gives `np.interp` an explicit decay anchor across the buffer zone and
the interpolation slopes smoothly to zero. The synthetic mode carries
sentinel values (`number = -1`, `mode_type = "synthetic_zero"`,
`precision = "synthetic"`) and omits the `_evg`/`_centers`/`_c2l`
keys so the B-factor accumulator and SS-analysis loop skip it
automatically. Equivalent in effect to the pre-v1.0.4 right=0
boundary, but explicit, logged, and integratable into the
modulation-spectrum and cumulative-Lambda sums (where it contributes
zero by construction).

Three regression tests in `tests/test_v104_audit_fixes.py`:
`test_synthetic_zero_mode_has_zero_observables`,
`test_synthetic_zero_mode_marked_as_synthetic`,
`test_synthetic_zero_mode_safe_in_interp`.

### Added — physics sanity tests (Phase 3 of v1.0.4 audit)

`tests/test_sanity_physics.py` adds **19 numerical sanity tests** that
verify the implementation matches the manual formulas via inputs with
analytically known answers. These complement the regression tests
(which protect against known bugs) by protecting against the future
class of bug where someone refactors a formula and the implementation
silently diverges from the physics. Coverage:

- **OOP fraction extremals**: pure-OOP → 100%, pure-INP → 0%,
  Pythagorean closure α_OOP + α_INP = 1.
- **Bond modulation invariants**: translation → 0, rotation → 0,
  pure stretching → |Δr| = sum of axial components.
- **λ_X at T→0**: pure stretching gives ℏω/4 = ω/4 in cm⁻¹.
- **u_rms limits**: zero-point at T→0, classical equipartition at
  T ≫ ℏω/k_B, monotonic in T.
- **HA reaction-coord**: acceptor-only motion → 0,
  H-toward-A → +d, N-toward-A → -d (the canonical PT discriminator).
- **RSS aggregation**: single sub-channel → |dr|; multi sub-channel
  λ-aggregation consistency.
- **λ general properties**: zero at dr=0, sign-independent, quadratic
  in dr, quadratic in ω.

Total test count after v1.0.4 audit:
- 16 regression tests (test_v104_audit_fixes.py)
- 19 sanity tests (test_sanity_physics.py)
- 3 manual-code-ref tests (test_manual_code_refs.py)
- + 110 pre-existing tests
= **148 unit tests + 4 E2E = 152 tests total**, all green.

### Documentation — formula-to-code references in Manual_EN

Manual_EN.tex now annotates every non-trivial formula in the
"Theoretical foundations" chapter with a `\codref{module.function}`
macro that points readers at the implementing function in the
source tree. The companion test `tests/test_manual_code_refs.py`
verifies that every referenced function exists in the code,
turning manual-vs-code drift into a CI failure. Coverage:
`core.compute_thermal_amplitude` (u_rms),
`core.classify_oop_inp` (OOP fraction),
`reorganization.signed_dr_along_axis` (classical Δr),
`reorganization.compute_mode_modulations` (HA reaction coord),
`reorganization.compute_modulation_spectra` (M_X(ω)),
`reorganization.aggregate_by_parent` (RSS aggregation),
`reorganization.compute_total_reorganization` (Λ_X^total),
`core.analyze_fe_ligand` (bend split),
`core.analyze_his_hn` / `core.analyze_fe_ligand` (σ propagation),
`geometry.cluster_normal` (cluster-normal SVD),
`core.compute_scsd_for_mode_full` (SCSD projection).

## [1.0.3] — 2026-05-12

Feature release. Closes the gap between data that the pipeline already
computed and data that the standard outputs actually surfaced. No
behavioural change for downstream analysis built on v1.0.2 outputs —
all v1.0.2 sheets, columns and PNG files remain unchanged.

### Added

- **Ca-UMAP embedding** (`embedding.compute_ca_umap_cluster`).
  Models each mode as a point in an N_Ca-dimensional space of C-alpha
  amplitudes (one dimension per backbone residue), then runs UMAP +
  HDBSCAN. Complements the existing global UMAP (Marcus-Hush
  reorganization features) and SS-UMAP (per-secondary-structure-element
  features). Provides a mode-shape fingerprint that highlights modes
  localized to the same region of the protein.
- **Ca-UMAP Excel export** — two new sheets in `_analysis_Embeddings.xlsx`:
  - `Ca_UMAP_clusters`: per-mode 2D coordinates + cluster ID + mode
    type + frequency.
  - `Ca_UMAP_profile`: Z-score profile of the most discriminative Ca
    residues (Fisher-F ranked, top 30) + representative modes per
    cluster.
- **SS-UMAP PNG** (`*_embedding_SS_UMAP.png`). Three-panel layout
  (frequency / mode type / HDBSCAN cluster). SS-UMAP was already
  computed in v1.0.0+ and exported as Excel data, but had no PNG
  rendering. Standard pipelines now emit it automatically when
  `cfg.export_embedding_plots = True`.
- **Ca-UMAP PNG** (`*_embedding_Ca_UMAP.png`). Same three-panel layout
  as SS-UMAP, distinct colour theme.
- **C-alpha amplitude heatmap PNG** (`*_ca_amplitudes_heatmap.png`).
  Log-scale heatmap of residue × frequency, summarising the entire
  backbone-amplitude matrix that was previously only available as a
  raw Excel sheet.
- **Coordination diagnostic sheet** in `_analysis.xlsx`. Lists every
  ligand group and the Gaussian-atom indices assigned to it, with
  per-atom rows (atom name, element, coordinates). Makes it trivial
  to audit the kind of issue that surfaced in the v1.0.2 follow-up
  fix: e.g. "did a crystal water with overlapping residue number
  sneak into a Cys group?" Previously this information was only
  available as a single line in REPORT.txt giving the atom *count*,
  not the atom *identity*.

### Why

The audit that motivated this release showed that the pipeline had
been computing more useful structure than it was exporting:

1. SS-UMAP coordinates went to Excel only — no PNG.
2. Ca-amplitude data went to Excel only — no UMAP, no PNG.
3. Group atom assignments were only counted (not enumerated), making
   the v1.0.2 HOH-overlap bug invisible until it triggered a knock-on
   symmetry violation (His 255 had 11 atoms, His 259 had 10).

The fix is purely additive: no existing artefact changes. All v1.0.2
callers continue to work unchanged because the new arguments are
optional with sensible defaults.

### Tests

- **`tests/test_ca_umap_and_exports.py`** (10 new tests, all passing):
  - `test_ca_umap_basic` — Ca-UMAP runs on synthetic data and returns
    coherent shape, finite coordinates, sensible feature names.
  - `test_ca_umap_orientation_autodetect` — accepts both
    `(n_calpha, n_modes)` and `(n_modes, n_calpha)` matrices.
  - `test_ca_umap_too_few_modes_returns_none` — graceful abort with
    warning when fewer than 5 modes have Ca data.
  - `test_ca_umap_handles_none_input` — clean return when `ca_data`
    is `None` (no Ca atoms in PDB).
  - `test_export_embedding_excel_writes_ca_umap_sheets` — Excel export
    creates `Ca_UMAP_clusters` (and `Ca_UMAP_profile` when clusters
    are found).
  - `test_export_embedding_excel_skips_ca_umap_when_none` — backward
    compatibility: omitting `ca_umap_data` must not create the new
    sheets.
  - `test_coordination_sheet_lists_each_ligand_with_atom_count` —
    Coordination sheet exists, has one summary row per ligand with
    correct atom count, and per-atom detail rows.
  - `test_export_embedding_plots_back_compat` — old call signature
    (no new args) does not raise.
  - `test_export_embedding_plots_renders_ca_umap_png` — actually
    produces a valid PNG.
  - `test_export_embedding_plots_renders_ca_heatmap_png` — heatmap
    PNG is written.

Total test count: 106 unit + 4 E2E smoke = 110 tests, all green.

### Implementation notes

- The Ca-UMAP is run on the **same valid-mode subset** as Ca-amplitudes;
  modes that have no Ca data (e.g. failed amplitude calculations) get
  a NaN row in the padded output, matching the SS-UMAP convention.
- `Z2d_full` from `compute_ca_umap_cluster` is shape `(n_modes, 2)` —
  full mode list, with NaN rows for invalid modes — so consumers can
  index it parallel to `results` without an extra mapping table.
- Feature normalization is z-score per residue (mean 0, std 1 across
  all valid modes). Without this, the few high-amplitude residues
  would dominate the UMAP distance and the embedding would degenerate
  to a low-rank projection of those few residues. Identical convention
  to SS-UMAP.
- All four new PNGs are emitted only when `cfg.export_embedding_plots
  = True` (same gate as the existing UMAP PNG). The Excel sheets are
  unconditional — they cost negligible additional run time.

## [1.0.2] — 2026-05-12

Bug-fix release. One real analysis bug fixed; two diagnostic warnings
added as defense in depth. The mathematical pipeline is unchanged.
Systems whose PDB files have no hydrogens interleaved between heavy
atoms (e.g. the validation system mitoNEET-H87C used in v1.0.0) produce
numerically identical results to v1.0.1.

### Fixed

- **Index-mismatch bug in `geometry._build_group_map`**.
  `pdb_to_gaus_h` is keyed by indices into `pdb_data["all_h"]` (the full
  atom list, hydrogens included), but `_build_group_map` enumerated a
  pre-filtered `pdb_heavy` list and used those positional indices to
  look up `pdb_to_gaus_h`. For PDBs where hydrogen atoms appear
  interleaved between heavy atoms (Apd1 QM/MM hessians, His-NH inline
  with backbone heavy atoms), the indices silently drift apart starting
  at the first interleaved hydrogen. This corrupts the per-residue
  center list for every residue after the drift, and in combination
  with the eigenvector centers map `c2l` in `core.py` produces residue
  rows that are 100% zero across all modes in the
  `Gruppen_OOP/INP/Winkel/Tors` sheets, with no warning.

  Symptom in Apd1 WT runs (v1.0.0 and v1.0.1):
  - Apd1-prot: `His 259` row is 100% zero in all four `Gruppen_*` sheets
  - Apd1-deprot: `Cys 207` row is 100% zero in all four `Gruppen_*` sheets

  The fix iterates `pdb_data["all_h"]` (the same list `pdb_to_gaus_h`
  was built from) and filters out hydrogens inside the comprehension.

- **Companion bug: crystal waters with overlapping residue numbers**.
  Once the index-mismatch fix above was applied, a second, previously
  hidden bug surfaced: `_build_group_map` filtered only on `rnum`
  (residue number), so a crystal water (HOH) placed in a different
  chain that happens to share a residue number with an amino-acid
  ligand (common in QM/MM PDB files, where waters are numbered
  separately and overlap freely with the protein chain) was silently
  pulled into the ligand's atom group, contaminating its OOP/INP
  values with a non-zero water contribution.

  In v1.0.0 and v1.0.1 the index-mismatch bug accidentally excluded
  these HOH atoms by hitting a stale index that fell outside
  `pdb_to_gaus_h`. With the correct indexing in v1.0.2, the HOH atoms
  were correctly mapped and thus incorrectly pulled into the group.
  The combined fix in v1.0.2 also filters on `rname` (residue type),
  requiring it to match the ligand's residue type (`HIS`, `CYS`, ...).
  Symptom: in the first v1.0.2 Apd1-prot run, the user observed
  `His 255: 11 atoms assigned` (10 amino-acid atoms + 1 stray HOH O)
  and `Cys 216: 7 atoms assigned` (6 + 1 stray HOH O); the second
  amendment brings these to the expected 10 and 6.

### Added (diagnostic safety nets)

- **Warning in `core.analyze_mode`** when a group in
  `coord_info.group_map` resolves to zero eigenvector indices. With
  the index-mismatch fix above this branch should be unreachable for
  canonical Cys/His systems, but the warning is kept as a safety net
  for future ligand types (Asp/Glu/Ser/Thr) or unusual PDBs. The
  warning is deduplicated via a module-level set
  (`_WARNED_EMPTY_GROUP`) so it fires at most once per group per run,
  not once per mode (which would mean thousands of duplicates).
- **Warning in `export._ws_gruppen`** when a residue row in
  `Gruppen_OOP`, `Gruppen_INP` or `Gruppen_Winkel` ends up all zero
  across all modes. This catches the same class of silent zero-row
  bug at the export stage, independent of where the data was lost.
  Torsion is excluded from this check because its raw values are
  legitimately small (~10⁻³) and an empty torsion list returns a
  meaningful 0.
- **New `core.reset_warning_state()` helper**, called automatically by
  `_run_analysis_single` before each cluster analysis, to reset the
  per-run warning deduplication state. Required for long-lived Python
  sessions (notebooks) that call `run_analysis` repeatedly.

### Tests

- **`tests/test_group_map_indexing.py`** (6 new tests, all passing):
  - `test_build_group_map_returns_all_heavy_atoms_with_interleaved_hydrogens` —
    direct regression test that reproduces the Apd1 index-mismatch
    scenario on a synthetic 9-atom system and verifies all heavy atoms
    are recovered after the fix
  - `test_build_group_map_no_interleaved_hydrogens_unchanged` —
    ensures the canonical "heavy first, then H" PDB layout still works
  - `test_build_group_map_excludes_overlapping_water` —
    reproduces the HOH-with-overlapping-rnum scenario on a synthetic
    11-atom system and verifies the HOH oxygen does NOT leak into the
    ligand's group_map
  - `test_warned_empty_group_dedup` — verifies the per-run deduplication
    of the `_WARNED_EMPTY_GROUP` set
  - `test_ws_gruppen_warns_on_all_zero_row` — verifies the export-stage
    warning fires for OOP/INP/Winkel but not Torsion
  - `test_ws_gruppen_no_warning_when_all_good` — verifies no false
    positives when all residues have nonzero values

Total test count: 96 unit + 4 E2E smoke = 100 tests, all green.

### Acknowledgements

The bug was found by side-by-side comparison of Apd1 wildtype runs
(protonated vs deprotonated `[2Fe-2S]` cluster) and a HH255_259CC Cys4
mutant. The independent diagnosis arrived at the same conclusion from
two angles: a code-pattern audit of silent `.get(_, {}).get(_, 0)`
fallbacks in `export.py` (third-party review), and an atom-list
divergence audit between `_build_group_map` and
`find_coordinating_residues` (this project). Both pointed at the same
root cause, which is reassuring: the fix sits where both audits met.

## [1.0.1] — 2026-05-07

Documentation cleanup and consistency release. **No analysis changes** —
the analysis pipeline is mathematically identical to v1.0.0; running
v1.0.1 produces numerically identical results. Two minor user-facing
text changes are noted under "Code" below.

### Added

- **`docs/Manual_EN.pdf`** — 95-page complete English translation of
  the German manual. Covers all theory chapters (NRVS, Marcus-Hush,
  Huang-Rhys, SCSD, UMAP/HDBSCAN), all 48 configuration fields, four
  worked workflow scenarios, validation chapter (model matrix +
  mitoNEET-H87C), troubleshooting, full Excel sheet reference, and
  bibliography (61 entries). Both `Manual_DE.pdf` and `Manual_EN.pdf`
  are now the authoritative full references; the 17-page `Manual.pdf`
  remains as a quick-start English overview.
- **"How to start a run" subsection** in both manuals (`Manual_EN.tex`
  §4.1 and `Manual_DE.tex` §4.1 "Wie startet man einen Lauf?"),
  documenting:
  - TOML as the first positional CLI argument (no `--config` flag).
  - Three valid forms for Windows paths in TOML strings (forward
    slashes, single-quote literals, escaped backslashes).
  - Explicit warning that Python's `r"..."` raw-string syntax is
    **not** valid TOML and produces a parser error when copied from
    Python code into a TOML file.
  - CLI overrides for parameter scans
    (e.g.\ `modenanalyse-2fe2s run.toml --temp-k 40`).

### Fixed (documentation)

- **README.md**: `--config <file>` flag examples replaced with the
  correct positional invocation `modenanalyse-2fe2s file.toml`. Added
  Windows TOML path pitfalls warning. Output structure listing
  corrected (`*_REPORT.txt` instead of `*_BEFUND.txt`; removed
  reference to the no-longer-existing `*_analysis_NIS.xlsx` workbook).
  Citation block reorganized to point primarily at the concept DOI.
- **`docs/QUICKSTART.md`** Step 6: rewritten with positional TOML CLI
  usage, three valid Windows path forms, and CLI override examples.
- **`example_run.py`**: corrected stale `--config full_template.toml`
  example to positional form.
- **`docs/Manual.tex`** (short EN): same `--config` correction; DOI
  block updated to use the concept DOI.
- **`examples/full_template.toml`**: header claim "53 configuration
  fields" corrected to 48 (the five removed fields belonged to the
  retired NIS spectrum subsystem). Path note expanded to show all
  three valid TOML path forms with the Python-`r"..."`-is-not-TOML
  warning.
- **`docs/Manual_DE.tex`**:
  - Same "53 Konfigurationsfelder" → 48 correction.
  - Stale filename `BEFUND.txt` replaced with the actual filename
    written by the code, `*_REPORT.txt`. The German domain term
    *Befund* (diagnostic report) is preserved as a word where used in
    sentences.
  - Removed reference to a "drittes Excel-Workbook (NIS-Excel)" that
    has not been produced by the runner since the NIS subsystem was
    retired before v1.0.0.
  - **Removed 11 obsolete NIS configuration fields** from the
    "Complete parameter reference" table (`analyze_nis`,
    `nis_lineshape`, `nis_fwhm_gauss`, `nis_fwhm_lorentz`,
    `nis_n_points`, `nis_freq_min`, `nis_freq_max`, `nis_phonon_order`,
    `nis_n_theta`, `nis_split_elastic`, `nis_inmemory_max_atom_modes`).
    These are silently dropped with a `UserWarning` when present in a
    legacy TOML, so listing them as active configuration was misleading.
- **`docs/SHEET_MAPPING_DE_EN.md`**: filename mapping table corrected
  (`*_BEFUND.txt` → `*_REPORT.txt` is no longer a per-version
  difference; both versions write `*_REPORT.txt`). Removed the stale
  `*_analysis_NIS.xlsx` row. Column headers updated from
  `(v1.0.0_de)`/`(v1.0.0_en)` to `(DE edition)`/`(EN edition)`, since
  the mapping describes a permanent property of the two language
  editions and not a version-specific feature.

### Fixed (code, user-facing strings only)

These two changes affect log/template text that the user sees, but
not any computed value or data flow.

- **`src/modenanalyse_2fe2s/runner.py`**: replaced the runtime log
  message `"v3.7 reorg modulations active …"` with
  `"Reorg modulations active …"`. The "v3.7" referred to an internal
  pre-1.0 development series and was confusing in v1.0.x logs.
- **`src/modenanalyse_2fe2s/cli.py`** (`_write_template`): the TOML
  template generated by `--write-template` previously had a brief
  Windows-path note ("forward slashes or doubled backslashes") and
  contained two stale comments about the retired NIS subsystem. The
  template now shows the three valid TOML path forms with an explicit
  `r"..."`-is-not-TOML warning, and no longer references NIS.
- **`pyproject.toml`** entry point name: changed from
  `modenanalyse_2fe2s` (underscore) to `modenanalyse-2fe2s` (hyphen),
  matching the form used everywhere in the documentation. Before
  this fix, following the documented invocation
  `modenanalyse-2fe2s run.toml` produced "command not found"; users
  had to know that pip had installed the script under the
  underscored name. The hyphen form has been the documented form
  since v1.0.0; this commit makes the actual entry point match.

### Cleaned up

- **`src/modenanalyse_2fe2s/config.py`**: removed multi-line German
  comment block describing the v3.x development history. The
  `__version__` line is now a clean three-line comment listing the
  v1.0.0 and v1.0.1 release notes only. The block stripped no code,
  only comments. The `_legacy_silent_drop` set itself is kept (it
  ensures legacy v3.6/v3.7 TOMLs still load with a `UserWarning`).
- **`src/modenanalyse_2fe2s/config.py`** `to_toml()` header: German
  comment "modenanalyse_2fe2s -- Konfiguration" / "Creates with
  to_toml() von Version" replaced with English equivalents, for
  consistency with the rest of the EN edition.
- **`install.bat` / `install.ps1`**: removed hardcoded "v1.0.0"
  from the install banner (the actual installed version is reported
  later from `__version__`); removed the `pip uninstall -y nisspec3`
  line (an internal pre-1.0 predecessor package that no end user
  ever installed); expanded the "First run" hint to show both the
  template-generation command and the actual run command.
- **`README.md`**: Python version badge corrected from "3.10+" to
  "3.11+", matching `pyproject.toml`'s `requires-python = ">=3.11"`.
- **`docs/QUICKSTART.md`**: Prerequisites updated from "Python 3.10
  or later" to "Python 3.11 or later".
- **`pyproject.toml`**: `Documentation` URL repointed from the
  17-page `docs/Manual.pdf` to the 95-page `docs/Manual_EN.pdf`,
  which is now the authoritative English reference. Removed the
  "Update once the GitHub repo is public" comment, which is no
  longer applicable.
- **`examples/multi_window_template.toml`**: removed the version
  stamp "(modenanalyse_2fe2s v1.0.0)" from the file header. This
  template is not version-specific.
- **`examples/full_template.toml`**: corrected stale comment "NIS
  and pDOS are still computed globally over all modes" to
  "the cluster/embedding analysis is run independently within each
  window" (the NIS subsystem was retired before v1.0.0).
- **`.zenodo.json`**: version field bumped to `"1.0.1"`.
- **`CITATION.cff`**: version bumped, primary DOI repointed to the
  concept DOI (which always resolves to the latest release); the
  v1.0.0 version-DOI is kept as a secondary identifier for explicit
  historical reference.

### Result

After v1.0.1, all user-facing surfaces (README, QUICKSTART, both full
manuals, the short Manual, tutorial, sheet-mapping reference, example
scripts, TOML templates, runtime log messages, and the
`--write-template` output) are mutually consistent and reflect the
actual behavior of the v1.0.x codebase. All 94 tests pass (90 unit +
4 end-to-end smoke tests).

## [1.0.0] — 2026-05-07 — First public release

Initial release. `modenanalyse_2fe2s` is a Python package for vibrational
spectroscopic analysis of [2Fe-2S] clusters from DFT frequency
calculations (Gaussian and ORCA), with focus on reorganization energies
for electron transfer (ET) and proton-coupled electron transfer (PCET).

### Main features

- **Mode analysis**: Eigenvector-based classification of all normal modes
  into in-plane / out-of-plane / torsional based on the cluster plane.
- **Reorganization energies**: Marcus-Hush reorganization per bonding
  channel (Fe-Fe, Fe-S, Fe-N, N-H, H-acceptor) and per mode, with
  low-temperature asymptotic
  $\lambda_X(i) \to (\hbar\omega_i/4) \cdot \alpha_X^2(i)$ at $T \to 0$.
- **PDOS / NIS-compatible spectra**: Mode-resolved partial vibrational
  densities at 0.5 cm⁻¹ resolution with Gaussian broadening.
- **SCSD**: Symmetry-coordinate structural decomposition of the cluster
  core into $D_{2h}$ irreps (Kingsbury & Senge, 2024).
- **PCET score**: Hydrogen-bond detection via donor-acceptor geometry,
  multi-mode coupling to the reaction coordinate.
- **Embeddings & clustering**: PCA / UMAP over 31 mode features,
  HDBSCAN clustering of the mode landscape.
- **Multi-cluster**: First-class support for dimers and multi-[2Fe-2S]
  systems via `analyze_all_clusters = true` — a single run analyzes all
  detected clusters automatically in separate subfolders.
- **Backend support**: Gaussian log files (`.log`) and ORCA Hessian
  (`.hess`).
- **Output**: Excel with 25 structured sheets (Origin-compatible),
  publication-quality matplotlib plots (300 DPI), `REPORT.txt` with a
  summarizing run report.

### Validation

Validated against mitoNEET-H87C QM/MM Hessian (Cys$_4$ ligands, 422
atoms, 1266 modes after filter, 11 UMAP clusters, characteristic
reorganization values $\Lambda_\text{FeFe} \approx 29$,
$\Lambda_\text{FeS} \approx 338$ cm⁻¹, $\Lambda_\text{HA} = 0$ as
expected for the H87C mutation removing the only histidine ligand).

### Tests

92+ tests (88+ unit tests + 4 end-to-end smoke tests), all passing.

### Documentation

Manual with ~24,700 words on ~100 pages: theory chapter with
Marcus-Hush and Huang-Rhys connections, worked example with complete
mode analysis run, all algorithms documented as pseudocode, full
Excel sheet reference, configuration reference for 48 fields.

### License

GNU General Public License v3.0 or later (GPL-3.0-or-later). When
using this software, please cite the reference given in `CITATION.cff`.
