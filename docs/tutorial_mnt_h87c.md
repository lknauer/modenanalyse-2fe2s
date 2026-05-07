# Worked Example: mitoNEET-H87C QM/MM Hessian

This tutorial walks through a complete analysis of a real
QM/MM-derived [2Fe-2S] system, showing every input, every output, and
the interpretation of the headline numbers. Use this to validate your
own installation against known-good results.

## The system

**mitoNEET (CDGSH iron-sulfur protein 1, CISD1) with the H87C mutation.**
mitoNEET in its wild-type form is a [2Fe-2S] cluster protein with an
unusual His-Cys$_3$ ligation pattern. The H87C mutation replaces the
single histidine ligand by cysteine, yielding a Cys$_4$-coordinated
cluster -- a ferredoxin-like ligation environment in an otherwise
mitoNEET-shaped binding pocket.

Why this is a useful test system:

- Cys$_4$ ligation means **no PCET pathway** -- the H-acceptor and
  N-H channels should be exactly zero. This is a sharp consistency
  check.
- The Fe-S$_4$ environment dominates the inner-sphere reorganization,
  giving a clear $\Lambda_\text{FeS}$ signal.
- The cluster geometry is well-resolved, and the system is large
  enough (422 atoms, 1266 modes) to exercise the full pipeline at
  realistic scale.

## Input files

Two files (both shipped with this tutorial as supplementary
references; you can also use your own QM/MM Hessian):

- `mNT_H87C_qmmm-orca_B3LYP.hess` -- ORCA-6 Hessian file from a
  B3LYP QM/MM frequency calculation.
- `mNT_H87C_qmmm-orca_B3LYP.pdb` -- the matching PDB for ligand
  identification and Kabsch alignment.

## Running the analysis

```python
from modenanalyse_2fe2s import Config, run_analysis

cfg = Config(
    log_file   = "mNT_H87C_qmmm-orca_B3LYP.hess",
    pdb_file   = "mNT_H87C_qmmm-orca_B3LYP.pdb",
    output_dir = "./mnt_h87c_results",
    temp_k     = 5.0,         # cryogenic NRVS standard
    freq_max   = 800.0,       # NRVS measurement window
    pdb_chain  = "",          # accept any chain
)
run_analysis(cfg)
```

This takes about 5 minutes on modern hardware. The console output ends
with:

```
  DONE  (4.8 min / 288 s)
======================================================================
  modes analyzed:   736
  HP eigvec:          736/736
  In-plane:           509
  Out-of-plane:       91
  Torsional:          136
  Warnings:          2
```

736 modes were analyzed (the rest filtered above the 800 cm$^{-1}$
window). The ratio In-plane:Out-of-plane:Torsional reflects the
overall mechanical character of the protein -- 509 in-plane modes
dominate, consistent with a soft protein matrix that mostly buckles
parallel to its surface.

## Output structure

```
mnt_h87c_results/
└── 0-800_cm-1/
    ├── mNT_H87C_qmmm-orca_B3LYP_REPORT.txt              # human-readable
    ├── mNT_H87C_qmmm-orca_B3LYP_analysis.xlsx           # 16 sheets
    ├── mNT_H87C_qmmm-orca_B3LYP_analysis_SS.xlsx
    ├── mNT_H87C_qmmm-orca_B3LYP_analysis_Embeddings.xlsx
    ├── mNT_H87C_qmmm-orca_B3LYP_analysis_interp0.05.xlsx
    ├── mNT_H87C_qmmm-orca_B3LYP_analysis_SS_interp0.05.xlsx
    └── mNT_H87C_qmmm-orca_B3LYP_embedding_UMAP.png
```

## Reading the REPORT

The REPORT file is the first thing to look at. Key lines:

### Geometry section

```
Fe-Fe:    2.7795  ± 0.0014 A
Fe1-S1:   2.2094  ± 0.0014 A
Fe2-S1:   2.2594  ± 0.0014 A
Fe1-S2:   2.2521  ± 0.0014 A
Fe2-S2:   2.3215  ± 0.0014 A
Cluster normal n_hat: (-0.8393, -0.1171, -0.5309)
Folding residual:     0.0530 A
```

These are the cluster equilibrium parameters. For a reduced [2Fe-2S]:

- **Fe-Fe ≈ 2.78 Å** is in the canonical range for reduced
  ferredoxin/Cys$_4$ clusters (literature: 2.69-2.80 Å).
- **Fe-S asymmetry**: the Fe-S distances span 2.21-2.32 Å. In an
  ideal $D_{2h}$ cluster they would all be equal; the asymmetry here
  reflects the protein-induced strain. Fe2-S2 is the longest at
  2.32 Å -- this slight elongation correlates with the trans Cys-87
  cysteine, which is the mutated position.
- **Folding residual = 0.053 Å** is small (literature criterion: a
  folded cluster has folding > 0.1 Å). This Cys$_4$ cluster is
  essentially planar.

### Coordination

```
Fe1 ← Cys 72 (S, SG, d=2.347 A)
Fe1 ← Cys 74 (S, SG, d=2.311 A)
Fe2 ← Cys 83 (S, SG, d=2.288 A)
Fe2 ← Cys 87 (S, SG, d=2.327 A)   ← the mutated position
PCET reorg: NOT active (no His)
```

All four cysteines are correctly identified, and PCET is reported as
inactive. This is the expected behavior for a Cys$_4$ cluster.

### Reorganization headline

```
Marcus-Hush total reorg (system sums):
  Lambda_FeFe = 28.68 cm-1
  Lambda_FeS  = 337.89 cm-1
```

These are the system totals $\Lambda_X$ summed over all 736 modes,
using the mode-reduced mass convention.

## The Reorg_total sheet

Open `mNT_H87C_qmmm-orca_B3LYP_analysis.xlsx` and navigate to
`Reorg_total`. You will see:

| Channel | Λ_pair (cm⁻¹) | Λ_mode (cm⁻¹) | n_modes |
|---------|---------------|---------------|---------|
| FeFe    | 800.7         | 28.68         | 736     |
| FeN     | 0             | 0             | 0       |
| FeS     | 6882.5        | 337.89        | 736     |
| NH      | 0             | 0             | 0       |
| HA      | 0             | 0             | 0       |

**Interpretation**:

- **$\Lambda_\text{FeFe}^\text{mode} = 28.7$ cm$^{-1}$** is the
  cluster-breathing reorganization. This is small in absolute terms
  but contributes to the Marcus inner-sphere reorganization for the
  intra-cluster electron transfer between the two Fe centers.

- **$\Lambda_\text{FeS}^\text{mode} = 337.9$ cm$^{-1}$** dominates.
  This reflects the Fe-S stretching reorganization expected when an
  electron localizes on or moves between the two Fe centers --
  essentially the inner-sphere $\lambda$ for ET to and from the
  cluster.

- **$\Lambda_\text{FeN}, \Lambda_\text{NH}, \Lambda_\text{HA}$ are all
  zero**: no histidine ligands → no Fe-N bonds → no PCET pathway.
  This is a strong consistency check; if any of these were nonzero on
  this Cys$_4$ system, it would indicate a bug in the ligand
  identification or the PCET geometry pipeline.

### About the two columns

The pair-reduced and mode-reduced columns differ by ~30×. This factor
is **not a numerical error** -- it reflects the conceptual difference
between:

- **Λ_pair** treats each bond as an isolated diatomic with reduced
  mass $\mu_\text{pair} = m_a m_b/(m_a + m_b)$. Useful for direct
  comparison to a model two-body harmonic oscillator, but ignores the
  collective character of normal modes in extended systems.

- **Λ_mode** uses the actual DFT mode reduced mass $\mu_\text{mode}$,
  which incorporates the inertia of all atoms participating in the
  mode. This is much larger than $\mu_\text{pair}$ for delocalized
  modes -- and since $\lambda \propto 1/\mu$ at fixed $\Delta r$, the
  mode value is correspondingly smaller.

For Marcus-Hush reorganization energies in collective protein modes,
**use Λ_mode**. The pair value is shown alongside for diagnostic and
educational purposes.

## Where the reorganization comes from

Open the `Lambda_cumulative` sheet, plot
`Lambda_FeS_mode_cum_cm1` versus frequency. You will see a roughly
sigmoidal curve with most of the rise between 200-500 cm$^{-1}$. This
spectral region is dominated by Fe-S stretches, exactly where the
Bergner-Pelmenschikov / Wang-bestiary literature places the
characteristic Fe-S stretching bands of [2Fe-2S] clusters.

Open `Modulation_spectra` and plot `M_FeS` versus frequency. The
peaks correspond to specific Fe-S stretching modes; their weighted
integrals give the cumulative $\Lambda$.

## SCSD decomposition

The `SCSD` sheet shows the cluster-core $D_{2h}$ symmetry-coordinate
decomposition for each mode. The header includes the methodology
citation:

```
METHOD: Kingsbury & Senge, Chem. Sci. 15, 13638 (2024).
Canonical D2h reference (Fe-Fe=2.73A, Fe-S=2.20A).
Axes: x=Fe-Fe, y=S-S, z=normal.
SCSD values are ORTHOGONAL and directly comparable between structures.
```

For most low-frequency modes (acoustic, ligand-shell), the cluster
core barely moves -- you will see "No core contribution" in the
`Core mode (SCSD)` column. The strongly cluster-localized modes (in
the 200-450 cm$^{-1}$ region) show clear $A_g$, $B_{1u}$, $B_{2g}$
etc. characters.

## UMAP mode landscape

Open `mNT_H87C_qmmm-orca_B3LYP_embedding_UMAP.png`. You will see a
2D projection of all 736 modes with two color schemes (left: by
frequency; right: by mode type). HDBSCAN finds 6 distinct clusters in
the embedding plus 17 noise points. The clusters group modes with
similar mechanical character; navigating from one cluster to another
in the embedding traces the spectrum from low-frequency soft motion
through ligand bending to high-frequency Fe-S stretching.

## Comparison with NRVS literature

Published NRVS data on mitoNEET (Gee, Pelmenschikov, Mons et al.,
Biochemistry 60, 2419 (2021)) report Fe-S stretching bands in the
280-340 cm$^{-1}$ region for the wild-type protein. The H87C mutant
is computational here and not directly compared in literature, but the
cumulative $\Lambda_\text{FeS}$ rising across this exact region is a
direct visual confirmation that the mode pipeline correctly captures
the dominant inner-sphere reorganization.

## Validation summary

For your own installation, running on this same QM/MM Hessian should
give you:

| Quantity                    | Expected value     | Tolerance |
|-----------------------------|--------------------|-----------|
| Number of modes analyzed    | 736                | exact     |
| Cluster Fe-Fe distance      | 2.7795 Å           | ±0.001 Å  |
| Folding residual            | 0.053 Å            | ±0.005 Å  |
| Cysteine ligands identified | 4                  | exact     |
| Λ_FeFe (mode)               | 28.7 cm⁻¹          | ±0.5 cm⁻¹ |
| Λ_FeS (mode)                | 337.9 cm⁻¹         | ±2 cm⁻¹   |
| Λ_FeN, Λ_NH, Λ_HA           | 0 exactly          | exact     |
| HDBSCAN clusters            | 6                  | ±1        |

If your numbers match within tolerance, the installation is working
correctly. If not, check the WARNINGS section of the REPORT and the
Kabsch RMSD; small differences (under 1 cm$^{-1}$) can result from
numerical-precision differences across BLAS implementations.

## What to do next

- For your own [2Fe-2S] system: copy the configuration above and
  replace the file paths.
- For Rieske-type proteins (His$_2$Cys$_2$): expect non-zero
  $\Lambda_\text{FeN}$ and (with H-bond acceptors nearby) non-zero
  $\Lambda_\text{HA}$. These quantify PCET capability.
- For dimers (e.g. glutaredoxin dimers with two clusters): set
  `analyze_all_clusters = True` in your `Config`.
- For comparative studies: collect Λ_FeS and Λ_FeFe for multiple
  systems and compare. The mode-reduced values are the right metric.

## Files referenced

All paths relative to package root:

- `tests/data/Cys_2Fe-2S_red3_hpfrq_opt.log.xz` -- a smaller
  test system (Cys$_4$ model, ~50 atoms) for fast CI testing
- `examples/full_template.toml` -- annotated configuration template
- `examples/multi_cluster_template.toml` -- multi-cluster example
- `docs/Manual.pdf` -- 17-page English overview
- `docs/Manual_full_EN.pdf` -- 103-page complete reference
