# Quickstart Guide

A practical introduction to `modenanalyse_2fe2s` in 15 minutes.

## What you'll learn

By the end of this guide, you will:

1. Have the package installed and tested
2. Run a complete analysis on a small example system
3. Understand the structure of the output files
4. Know how to interpret the key results (reorganization energies, mode classifications)
5. Be ready to analyze your own [2Fe-2S] cluster systems

## Prerequisites

- Python 3.11 or later
- A Gaussian-16 `.log` file (with `freq=hpmodes`) **or** an ORCA-6 `.hess` file
- (Optional) A matching PDB file for ligand detection
- ~500 MB free disk space for dependencies

## Step 1: Installation (2 minutes)

```bash
# From the unpacked release directory:
pip install modenanalyse_2fe2s-1.0.1-py3-none-any.whl
```

This pulls in `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`,
`scikit-learn`, `umap-learn`, `hdbscan`, and `scsdpy`.

Verify:

```bash
modenanalyse-2fe2s --help
```

You should see the help message with all CLI options.

## Step 2: Run on a test system (3 minutes)

The package ships with a small test [2Fe-2S] cluster (Cys₄ ligation,
~50 atoms) under `tests/data/`. Try a quick analysis:

```python
from modenanalyse_2fe2s import Config, run_analysis

cfg = Config(
    log_file   = "tests/data/Cys_2Fe-2S_red3_hpfrq_opt.log.xz",  # xz-compressed test data
    output_dir = "./quickstart_results",
    temp_k     = 5.0,        # low-temperature NRVS standard
    freq_max   = 800.0,      # NRVS measurement window
)
run_analysis(cfg)
```

This takes about 1 minute and produces ~5 output files.

For full systems with a PDB file:

```python
cfg = Config(
    log_file   = "your_freq.log",
    pdb_file   = "your_structure.pdb",
    output_dir = "./results",
    temp_k     = 5.0,
    freq_max   = 800.0,
)
run_analysis(cfg)
```

## Step 3: Understand the output (5 minutes)

After running, look in `./quickstart_results/0-800_cm-1/`:

```
0-800_cm-1/
├── Cys_2Fe-2S_red3_hpfrq_opt_REPORT.txt              # human-readable summary
├── Cys_2Fe-2S_red3_hpfrq_opt_analysis.xlsx           # main results (11 sheets)
├── Cys_2Fe-2S_red3_hpfrq_opt_analysis_Embeddings.xlsx
├── Cys_2Fe-2S_red3_hpfrq_opt_analysis_interp0.05.xlsx
└── Cys_2Fe-2S_red3_hpfrq_opt_embedding_UMAP.png      # 2D mode landscape
```

### The REPORT file

Open `*_REPORT.txt` first — this is the human-readable summary. Sections:

- **CONFIGURATION**: what was run with which parameters
- **GEOMETRY**: detected Fe-Fe, Fe-S distances, Fe-S-Fe angles, cluster
  normal direction, planarity (folding residual)
- **COORDINATING AMINO ACIDS**: detected Cys/His/etc. ligands (needs PDB)
- **MODE ANALYSIS**: how many modes were found, filtered, analyzed
- **MODE DISTRIBUTION**: in-plane / out-of-plane / torsional partitioning
- **CLUSTER ANALYSIS**: UMAP/HDBSCAN clusters in the mode landscape
- **NOTES**: detailed run information including reorganization totals

A typical NOTES line you'll want to look at:

```
[i] Marcus-Hush total reorg (system sums): Lambda_FeFe=21.62 cm-1, Lambda_FeS=206.39 cm-1
```

These are the **system-level reorganization totals** Λ_X for each
bonding channel — the central quantitative result of the analysis.

### The Excel workbook

Open `*_analysis.xlsx` in Excel, LibreOffice, or Origin. The workbook
has 11 sheets:

| Sheet                     | What it contains                                                  |
|---------------------------|-------------------------------------------------------------------|
| `Mode_analysis`           | One row per mode: freq, sym, type, OOP%/INP%, eigenvector amps   |
| `Core_scores`             | Cluster-core mode scores (4 cluster atoms only)                   |
| `Equilibrium_distances`   | All Fe-Fe, Fe-S, Fe-N, etc. distances                             |
| `SCSD`                    | D2h symmetry-coordinate decomposition (Kingsbury & Senge 2024)    |
| `Reorganization_energy`   | Per-mode λ_X(i) and Δr_X(i) for all 5 channels                    |
| `Reorg_total`             | **System totals Λ_X = sum_i λ_X(i)** — the headline numbers      |
| `Reorg_per_bond`          | Per-bond reorganization (each individual Fe-S bond etc.)          |
| `Modulation_spectra`      | Frequency-resolved M_X(ω) — for plotting                          |
| `Lambda_cumulative`       | Cumulative Λ_X(ω) — for plotting                                  |
| `B_factors`               | Debye-Waller B-factors per atom                                   |
| `Info`                    | Run metadata                                                       |

### Reading the Reorg_total sheet

This is the most important sheet. Open it:

```
Channel  Lambda_total_pair_cm1  Lambda_total_mode_cm1  n_modes_contributing
FeFe     54.12                  21.62                  66
FeN       0.00                   0.00                   0
FeS     286.60                 206.39                  66
NH        0.00                   0.00                   0
HA        0.00                   0.00                   0
```

Interpretation for this Cys₄ system:

- **Λ_FeFe ≈ 22 cm⁻¹**: cluster-breathing reorganization (relevant for ET)
- **Λ_FeS ≈ 206 cm⁻¹**: dominant Fe-S stretching reorganization (also ET)
- **Λ_FeN = 0**: no histidine ligands (Cys₄ system)
- **Λ_NH = 0, Λ_HA = 0**: no PCET pathway (no His present)

For a Rieske-type [2Fe-2S] (His₂Cys₂), Λ_FeN and Λ_HA would be non-zero,
indicating PCET capability.

The two columns differ in convention:

- `Lambda_total_pair_cm1`: uses the **pair-reduced mass** μ_pair (pure
  bond-vibration view; same units as a diatomic harmonic oscillator)
- `Lambda_total_mode_cm1`: uses the **mode-reduced mass** μ_mode (the
  full DFT normal mode reduced mass, includes coupling to other atoms)

For **Marcus-Hush comparisons across systems**, use `_mode_cm1`.
For **bond-localized intuition**, use `_pair_cm1`.

## Step 4: Multi-cluster systems (2 minutes)

For dimers and multi-[2Fe-2S] systems (glutaredoxin dimers, multi-domain
constructs):

```python
cfg = Config(
    log_file   = "dimer_freq.hess",
    pdb_file   = "dimer.pdb",
    output_dir = "./results_dimer",
    analyze_all_clusters = True,    # <-- this is the key flag
    temp_k     = 5.0,
    freq_max   = 800.0,
)
run_analysis(cfg)
```

The package detects all [2Fe-2S] clusters automatically and produces:

```
results_dimer/
├── cluster_0/0-800_cm-1/    # full analysis for cluster #0
├── cluster_1/0-800_cm-1/    # full analysis for cluster #1
└── multi_cluster_summary.txt
```

The `multi_cluster_summary.txt` lists each cluster's geometry and run
status. Each cluster gets its own complete workbook + REPORT.

For single-cluster systems, the `analyze_all_clusters=True` flag falls
back automatically to single-cluster mode.

## Step 5: Multi-window analysis (2 minutes)

NRVS spectra are typically interpreted in three frequency bands. Run all
three at once:

```python
cfg = Config(
    log_file   = "system.log",
    output_dir = "./results",
    freq_windows = [
        (0,   100),     # acoustic / bending
        (100, 300),     # Fe-Fe + Fe-S region
        (300, 500),     # ligand modes
        (0,   500),     # full window for context
    ],
    temp_k = 5.0,
)
run_analysis(cfg)
```

You get one subfolder per window plus an overall analysis at the top
level. Useful for separating contributions of different vibrational
classes.

## Step 6: Use the CLI for production (1 minute)

For reproducible runs, use a TOML configuration:

```bash
modenanalyse-2fe2s --write-template my_run.toml
# edit my_run.toml as needed
modenanalyse-2fe2s my_run.toml
```

The TOML file path is the first **positional** argument; no `--config`
flag is needed.

### Pitfalls with paths in TOML

TOML strings with double quotes interpret backslashes as escape
characters, so Windows paths like `"D:\Data\file.log"` will fail
parsing. Three valid alternatives:

```toml
# 1. Forward slashes (works on Windows too)
log_file = "D:/Data/file.log"

# 2. Single quotes = TOML literal string (no escaping)
log_file = 'D:\Data\file.log'

# 3. Escaped backslashes
log_file = "D:\\Data\\file.log"
```

Note: Python's raw-string syntax `r"..."` is **not** valid TOML --
remove the `r` prefix when copying paths from Python code into a
TOML file.

### CLI overrides

CLI arguments can override individual TOML fields without editing
the file:

```bash
# Use my_run.toml but force temp_k = 40 K
modenanalyse-2fe2s my_run.toml --temp-k 40.0

# Use my_run.toml but write to a different output directory
modenanalyse-2fe2s my_run.toml --output-dir ./results_test
```

This is useful for parameter scans where most fields stay the same.

## What's next?

- **For NRVS interpretation**: examine the `Modulation_spectra` sheet to
  identify which modes contribute to which spectral band, and the
  `Reorg_per_bond` sheet to see which individual bonds dominate the
  reorganization.

- **For PCET analysis** (Rieske-type proteins): with His ligands present,
  the `NH` and `HA` channels become active. The `HA` (H-acceptor)
  channel is the proton-transfer reaction coordinate; non-zero Λ_HA
  indicates PCET capability.

- **For comparative studies**: run the analysis on multiple systems and
  compare the `Reorg_total` sheets. The `Lambda_total_mode_cm1` column
  is the right comparison metric.

- **For SCSD analysis**: examine the `SCSD` sheet for the D2h
  symmetry-coordinate decomposition of the cluster-core distortion.
  This is most useful for comparing mutants or oxidation states.

- **For mode landscape exploration**: the `*_embedding_UMAP.png` shows
  all modes in a 2D projection of the 31-feature space. Modes close
  together share similar mechanical character.

## Troubleshooting

- **"No [2Fe-2S] cluster found"**: check that your log file contains
  Standard orientation coordinates. For ORCA, ensure the `.hess` file
  has all required blocks.

- **"Kabsch alignment failed"**: your PDB may have a different chain
  identifier than the default `'A'`. Set `pdb_chain="X"` (or empty
  string `""` to allow all chains).

- **"openpyxl missing"**: `pip install openpyxl` — required for Excel
  output.

- **`Lambda_FeN` and `Lambda_HA` are zero on a Rieske protein**: check
  that the PDB has His ligands present, that the PDB is correctly
  aligned to the QM region (REPORT shows Kabsch RMSD), and that
  `pcet_enabled=True` (the default).

## Where to learn more

- `README.md` — package overview and feature list
- `docs/Manual.pdf` -- English manual (~17 pages) with theory,
  configuration reference, output reference, multi-cluster workflow,
  validation, and troubleshooting (recommended starting point)
- `docs/Manual_DE.pdf` -- complete reference manual in German (96 pages)
- `docs/SHEET_MAPPING_DE_EN.md` — sheet/column mapping if you also use
  the German edition
- `examples/full_template.toml` — annotated configuration template
- `examples/multi_cluster_template.toml` — multi-cluster TOML example
- `tests/` — small examples that can be adapted

## Citation

If you use this package in a publication, please cite as in
`CITATION.cff`. After Zenodo upload, the package will have a permanent
DOI for citation.
