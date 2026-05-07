# Changelog

All notable changes to `modenanalyse_2fe2s` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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
