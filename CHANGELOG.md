# Changelog

All notable changes to `modenanalyse_2fe2s` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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
