# -*- coding: utf-8 -*-
"""
example_run.py -- example script for modenanalyse_2fe2s.

Typical Spyder/IPython/Jupyter use:

  1. Install package (see README.md)
  2. Copy this file and adjust paths to your data
  3. In Spyder:   %runfile <path>/example_run.py
     In IPython:  exec(open("example_run.py").read())
     From shell:  python example_run.py

Alternatively via TOML configuration (recommended for reproducible
runs -- see examples/full_template.toml):

    modenanalyse-2fe2s --config full_template.toml
"""

from modenanalyse_2fe2s import Config, run_analysis


# =============================================================================
# Adjust your data paths here:
# =============================================================================
cfg = Config(
    # Required
    log_file   = r"PATH/TO/gaussian_freq.log",   # Gaussian .log with freq=hpmodes
    output_dir = r"PATH/TO/output_directory",    # target folder for Excel/PNG/REPORT

    # Optional, but recommended
    pdb_file   = r"PATH/TO/structure.pdb",       # PDB for ligand detection

    # Frequency range (None = unbounded)
    freq_min   = None,
    freq_max   = 800.0,                          # 800 cm-1 = NRVS measurement range

    # Sample temperature (typical values:
    #   5 K = low-temperature NRVS (liquid He),
    #  40 K = standard NRVS,
    #  80 K = cryogenic,
    # 300 K = room temperature)
    temp_k     = 5.0,

    # Multi-window mode: each window gets its own subfolder with
    # separate analysis.xlsx and reorg aggregates, plus an overall
    # analysis at top level.
    freq_windows = [
        (0,   100),
        (100, 300),
        (300, 500),
        (0,   500),
    ],

    # Multi-cluster selection (default 0 = closest Fe-Fe pair)
    cluster_index = 0,

    # Multi-cluster mode (NEW in v1.0.0):
    # Set to True to automatically analyze ALL [2Fe-2S] clusters in
    # the system. One subfolder cluster_0/, cluster_1/, ... is
    # created per cluster + multi_cluster_summary.txt in the base
    # folder with status per cluster.
    # For single-cluster systems the flag has no effect (falls back
    # automatically to single-cluster mode).
    # analyze_all_clusters = True,

    # Modules
    analyze_scsd = True,    # SCSD D2h symmetry decomposition
    pcet_enabled = True,    # PCET reorg pipeline (His-mediated)
)
# =============================================================================


if __name__ == "__main__":
    run_analysis(cfg)
