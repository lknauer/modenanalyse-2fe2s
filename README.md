# modenanalyse_2fe2s

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Version](https://img.shields.io/badge/version-1.0.2-green.svg)](CHANGELOG.md)
[![Tests](https://img.shields.io/badge/tests-94%20passing-brightgreen.svg)](tests/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20074194.svg)](https://doi.org/10.5281/zenodo.20074194)

**Mode-resolved vibrational analysis for [2Fe-2S] clusters**, with focus on
reorganization energies for electron transfer (ET) and proton-coupled
electron transfer (PCET).

## What this tool does

`modenanalyse_2fe2s` analyzes DFT frequency calculations of [2Fe-2S] cluster
proteins and produces a publication-ready breakdown of the vibrational mode
landscape. It targets [2Fe-2S] clusters of any ligation pattern (Rieske-type,
ferredoxin-type, mitoNEET, etc.) and is backend-agnostic with support for
both Gaussian (`.log`) and ORCA (`.hess`) outputs.

For each normal mode, the tool computes:

- **In-plane / out-of-plane / torsional classification** based on
  eigenvector projection onto the cluster plane.
- **Marcus-Hush reorganization contributions** $\lambda_X(i)$ per
  bonding channel $X \in \{\text{Fe-Fe}, \text{Fe-S}, \text{Fe-N},
  \text{N-H}, \text{H-acceptor}\}$, with low-temperature asymptotic
  $\lambda_X(i) \to (\hbar\omega_i/4) \cdot \alpha_X^2(i)$ at $T \to 0$.
- **Partial vibrational density of states** ($^{57}$Fe-PVDOS) compatible
  with NRVS / NIS spectroscopic measurements.
- **SCSD symmetry-coordinate decomposition** of the cluster core
  into $D_{2h}$ irreps (after Kingsbury & Senge, 2024).
- **PCET score** with hydrogen-bond detection via donor-acceptor geometry.
- **31-feature mode embeddings** with UMAP and HDBSCAN clustering.

System-level aggregates (total $\Lambda_X = \sum_i \lambda_X(i)$ per channel)
plus full per-mode breakdowns are written to multi-sheet Excel workbooks
(Origin-compatible; a typical run produces 11–16 sheets in the main
analysis workbook depending on whether a PDB is supplied) with
publication-quality matplotlib plots (300 DPI).

## Multi-cluster support (first-class)

For dimers and multi-[2Fe-2S] systems (glutaredoxin dimers, multi-domain
constructs, etc.), set `analyze_all_clusters = true` in the config. The
tool detects all [2Fe-2S] clusters automatically and runs the full
pipeline once per cluster, with separate output subfolders
`cluster_0/`, `cluster_1/`, ... and a `multi_cluster_summary.txt` listing
each cluster's geometry and run status.

## Installation

```bash
pip install modenanalyse_2fe2s-1.0.1-py3-none-any.whl
```

Or from source:

```bash
git clone <repository>
cd modenanalyse_2fe2s
pip install -e .
```

Windows users: PowerShell installer in `install.ps1`, batch file in
`install.bat`.

## Quick start

For a hands-on 15-minute introduction, see [`docs/QUICKSTART.md`](docs/QUICKSTART.md).

### From a config file (recommended)

```bash
modenanalyse-2fe2s examples/full_template.toml
```

The TOML file is the first positional argument — no `--config` flag needed.

> **Windows users:** in TOML, paths must use either forward slashes
> (`"D:/Data/file.log"`) or single quotes for backslashes
> (`'D:\Data\file.log'`). Python-style raw strings (`r"..."`) are
> **not** valid TOML and will produce a parser error.

### Programmatic API

```python
from modenanalyse_2fe2s import Config, run_analysis

cfg = Config(
    log_file   = "system_freq.log",
    pdb_file   = "system.pdb",
    output_dir = "./results",
    temp_k     = 5.0,             # low-temperature NRVS standard
    freq_max   = 800.0,           # cm-1
)
run_analysis(cfg)
```

### Command-line interface (without TOML)

```bash
modenanalyse-2fe2s \
    --log-file system.log \
    --pdb system.pdb \
    --output-dir ./results \
    --temp-k 5.0 \
    --freq-max 800.0
```

CLI arguments can also override individual fields when a TOML is given:

```bash
modenanalyse-2fe2s run.toml --temp-k 40.0   # uses run.toml but at 40 K
```

For multi-cluster systems:

```bash
modenanalyse-2fe2s examples/multi_cluster_template.toml
```

## Output structure

```
results/
└── 0-800_cm-1/
    ├── system_REPORT.txt                 # human-readable run report
    ├── system_analysis.xlsx              # main workbook (16 sheets)
    ├── system_analysis_Embeddings.xlsx   # UMAP + HDBSCAN
    ├── system_analysis_interp0.05.xlsx   # interpolated pDOS
    └── plots/
        ├── lambda_total_per_channel.png
        ├── pdos_combined.png
        └── ...
```

For multi-cluster runs:

```
results/
├── cluster_0/0-800_cm-1/...
├── cluster_1/0-800_cm-1/...
└── multi_cluster_summary.txt
```

## Documentation

The documentation comes in three complementary formats:

- **`docs/Manual.pdf`** (17 pages) -- a purpose-written English overview
  covering theory, configuration, output reference, multi-cluster
  workflow, validation, and troubleshooting. Recommended starting point
  for new users.
- **`docs/Manual_EN.pdf`** (95 pages) -- the complete English reference
  manual, with full theory chapters (NRVS, Marcus-Hush, Huang-Rhys,
  SCSD, UMAP/HDBSCAN), all 48 configuration fields, four worked
  workflow scenarios, validation chapter (model matrix +
  mitoNEET-H87C), troubleshooting, full Excel sheet reference, and
  bibliography (61 entries).
- **`docs/Manual_DE.pdf`** (98 pages) -- the original German reference
  manual; the same content as `Manual_EN.pdf` but with the German
  glossary and the *Befund*-Datei terminology preserved.

For a 15-minute hands-on introduction, see
[`docs/QUICKSTART.md`](docs/QUICKSTART.md).

For a complete worked example with real numbers from a mNT-H87C
QM/MM Hessian (422 atoms, 1266 modes), see
[`docs/tutorial_mnt_h87c.md`](docs/tutorial_mnt_h87c.md). This shows
every input, every output, and how to validate your installation
against known-good headline numbers (Λ_FeS ≈ 338 cm⁻¹, etc.).

The original German manual (`docs/Manual_DE.pdf`, 102 pages) is also
shipped for users who prefer the German source.

**Note on internal language**: Most source-code comments are now in
English; some longer passages may still contain German residuals from
the original development. Public API names (functions, classes,
parameters) are all in English.

A German edition of this package (`modenanalyse_2fe2s_v1.0.1_de.zip`) is
available with German user-facing strings. Both editions produce
**numerically identical** output for the same input. See
`docs/SHEET_MAPPING_DE_EN.md` for a side-by-side mapping of sheet names,
column headers, and classification labels — useful for porting Origin
or Excel post-processing scripts between the two editions.

## Validation

Validated against mitoNEET-H87C QM/MM Hessian (Cys$_4$ ligation, 422 atoms,
1266 modes after filter, 11 UMAP clusters). Characteristic reorganization
totals: $\Lambda_\text{FeFe} \approx 29$, $\Lambda_\text{FeS} \approx 338$
cm⁻¹, $\Lambda_\text{HA} = 0$ as expected for the H87C mutant
(no histidine ligand left).

## Tests

92+ tests covering geometry, reorganization, multi-cluster wrapper, and
end-to-end smoke tests:

```bash
python -m pytest tests/                # unit tests (~5 sec)
python -m pytest tests/ -m slow        # smoke tests (~80 sec)
```

## Citation

If you use `modenanalyse_2fe2s` in your research, please cite the
software using its concept DOI (which always resolves to the latest
released version):

> Knauer, L. (2026). *modenanalyse_2fe2s: Mode-resolved bonding
> reorganization analysis for [2Fe-2S] iron-sulfur clusters from DFT
> frequency calculations* [Software]. Zenodo.
> https://doi.org/10.5281/zenodo.20074194

For citing a specific version, use the version-DOI listed on the
[Zenodo record](https://doi.org/10.5281/zenodo.20074194). The DOI
for v1.0.0 is [10.5281/zenodo.20074193](https://doi.org/10.5281/zenodo.20074193).

BibTeX:

```bibtex
@software{knauer_modenanalyse_2fe2s_2026,
  author       = {Knauer, Lukas},
  title        = {{modenanalyse\_2fe2s: Mode-resolved bonding
                   reorganization analysis for [2Fe-2S] iron-sulfur
                   clusters from DFT frequency calculations}},
  month        = may,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {1.0.1},
  doi          = {10.5281/zenodo.20074194},
  url          = {https://doi.org/10.5281/zenodo.20074194}
}
```

A machine-readable citation entry is provided in `CITATION.cff`.

## License

GNU General Public License v3.0 or later (GPL-3.0-or-later). See
`LICENSE` for details.

## Author

Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau.
Vibrational analysis of [2Fe-2S] clusters in Rieske-type proteins
and mitoNEET.
