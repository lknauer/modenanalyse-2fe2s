# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Part of modenanalyse_2fe2s — see LICENSE in repository root.

# -*- coding: utf-8 -*-
"""
modenanalyse-2fe2s -- Command-line wrapper.

Three call patterns:

1. **TOML config** (recommended):
    modenanalyse-2fe2s run.toml

2. **TOML with override** (e.g., only change output-dir):
    modenanalyse-2fe2s run.toml --output-dir results_v2/

3. **Pure CLI args** (without TOML):
    modenanalyse-2fe2s --log-file dimer.log --output-dir results --temp-k 40

Write a template TOML:
    modenanalyse-2fe2s --write-template run_template.toml

Programmatic call from Spyder/Jupyter still works directly:
    >>> from modenanalyse_2fe2s import Config, run_analysis
    >>> cfg = Config.from_toml("run.toml")     # or Config(...)
    >>> run_analysis(cfg)
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Sequence

from .config import Config, __version__
from .runner import run_analysis


def _parse_freq_window(s: str) -> tuple:
    """Parses ein frequency window ``"lo-hi"`` -> ``(lo, hi)``."""
    parts = s.split("-")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"frequency window muss Format 'lo-hi' haben, not {s!r}.")
    try:
        return (float(parts[0]), float(parts[1]))
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"frequency window {s!r}: {e}") from e


def make_parser() -> argparse.ArgumentParser:
    """Returns the argparse parser for ``modenanalyse-2fe2s``."""
    p = argparse.ArgumentParser(
        prog="modenanalyse-2fe2s",
        description=("Quantum-chemical normal-mode analysis of "
                     "[2Fe-2S] iron-sulfur proteins.\n\n"
                     "Usage:\n"
                     "  modenanalyse-2fe2s run.toml          # TOML config\n"
                     "  modenanalyse-2fe2s --log-file ... --output-dir ...  # CLI config\n"
                     "  modenanalyse-2fe2s --write-template run.toml   # Create template"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version",
                   version=f"%(prog)s {__version__}")

    # TOML config (primary)
    p.add_argument("config", nargs="?", default=None,
                   help="Path to a TOML configuration file. "
                        "If given, takes precedence over CLI args; "
                        "CLI args can override individual fields.")

    # Template
    p.add_argument("--write-template", metavar="PATH", default=None,
                   help="Writes an annotated template TOML to PATH "
                        "and exits (no run).")

    # Input (required in pure-CLI mode, optional as override)
    g_in = p.add_argument_group(
        "Input (required without TOML, optional with TOML as override)")
    g_in.add_argument("--log-file", "-l", default=None,
                      help="Path to Gaussian .log file (freq=hpmodes).")
    g_in.add_argument("--output-dir", "-o", default=None,
                      help="Output directory for results.")
    g_in.add_argument("--pdb", "-p", default=None, dest="pdb_file",
                      help="PDB file with the same cluster (ligand resolution).")
    g_in.add_argument("--pdb-chain", default=None,
                      help="Chain identifier in the PDB. Default 'A'.")

    # Frequency (all as optional override)
    g_fr = p.add_argument_group("Frequency range (override)")
    g_fr.add_argument("--freq-min", type=float, default=None,
                      help="Lower frequency limit in cm^-1.")
    g_fr.add_argument("--freq-max", type=float, default=None,
                      help="Upper frequency limit in cm^-1.")
    g_fr.add_argument("--freq-window", action="append", default=None,
                      type=_parse_freq_window, metavar="LO-HI",
                      help=("Frequency windows (multi-window mode). "
                            "Can be specified multiple times. Example: "
                            "--freq-window 0-100 --freq-window 100-300"))

    # Thermal
    g_th = p.add_argument_group("Thermodynamics (override)")
    g_th.add_argument("--temp-k", type=float, default=None,
                      help="Sample temperature in Kelvin. Default 80.")

    # Cluster
    g_cl = p.add_argument_group("Cluster selection (override)")
    g_cl.add_argument("--cluster-index", type=int, default=None,
                      help=("Index of the cluster to analyze "
                            "(0 = closest Fe-Fe pair)."))

    # Toggle
    g_t = p.add_argument_group("Module on/off (override)")
    g_t.add_argument("--no-nis",  action="store_true",
                     help="(Deprecated, ignored) NIS spectra are not "
                          "produced by this tool; use a dedicated NIS "
                          "package on the same log file.")
    g_t.add_argument("--no-scsd", action="store_true",
                     help="Skip SCSD analysis.")
    g_t.add_argument("--no-pcet", action="store_true",
                     help="Skip PCET reorg pipeline (NH/HA channels).")
    g_t.add_argument("--no-cache", action="store_true",
                     help="Ignore file cache (full scan).")

    return p


def _write_template(path: str) -> None:
    """Writes an annotated template TOML."""
    template = f"""# modenanalyse_2fe2s -- configuration template
# Created with --write-template (v{__version__}).
#
# Fields with '#' at the start are commented out (default values active).
#
# Paths on Windows: TOML strings in double quotes interpret backslashes
# as escape characters. Use one of these three forms:
#   "D:/data/system.log"      forward slashes (works on Windows too)
#   'D:\\data\\system.log'      single quotes = TOML literal string (no escaping)
#   "D:\\\\data\\\\system.log"    escaped backslashes
# Note: Python's r"..." raw-string syntax is NOT valid TOML.

[input]
log_file   = "PATH/TO/gaussian.log"
output_dir = "PATH/TO/output_directory"
pdb_file   = "PATH/TO/structure.pdb"     # optional, but recommended
# pdb_chain = "A"
# logname_suffix = ""                    # for multiple runs on the same log file

[freq]
# Lower/upper filter limit for the mode analysis (in cm^-1).
# freq_min = 0.0
freq_max = 800.0

# Optional: multi-window mode. One Excel subfolder per window;
# the cluster/embedding analysis is run independently within each window.
# freq_windows = [[0, 100], [100, 300], [300, 500], [500, 700]]

[thermo]
temp_k = 40.0       # sample temperature in Kelvin

[cluster]
# For multi-[2Fe-2S] systems: 0 = closest cluster, 1 = next-closest, ...
cluster_index = 0
# fe_s_cutoff   = 3.0
# fe_fe_cutoff  = 3.5

[oop]
# Binary OOP/INP threshold in percent
# mode_type_threshold = 60.0
# Three thresholds for 7-level classification
# mode_type_detail_thresholds = [60.0, 75.0, 90.0]

[scsd]
analyze_scsd = true              # requires installed scsdpy
analyze_sse   = true

[pcet]
pcet_enabled               = true
# For purely Cys-coordinated clusters (e.g. HH->CC mutation) the
# HA modulation calculation is skipped (no H, no PT).
# pcet_hbond_cutoff_a      = 4.0    # H-bond acceptor search radius (A)
# pcet_acceptor_r0_a       = 2.8    # optimal distance for Gauss weighting
# pcet_acceptor_sigma_a    = 0.4    # tolerance of Gauss weighting
# reorg_spectrum_sigma_cm1 = 5.0    # broadening of M_X(omega) spectra
# reorg_spectrum_step_cm1  = 0.5    # frequency grid step

[embedding]
# Write embedding PNG (one per UMAP analysis)
export_embedding_plots = true
# umap_n_neighbors = 15        # auto heuristic if omitted

[numerics]
# include_hydrogen   = true
# interp_step        = 0.05
# use_cache          = true
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(template)
    print(f"Template written to: {path}")
    print("Adjust fields, then start with:")
    print(f"  modenanalyse-2fe2s {path}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point for ``modenanalyse-2fe2s`` (see pyproject.toml)."""
    args = make_parser().parse_args(argv)

    # Write template?
    if args.write_template:
        _write_template(args.write_template)
        return 0

    # Choose config source: TOML or pure CLI args
    if args.config:
        try:
            cfg = Config.from_toml(args.config)
        except Exception as exc:
            print(f"Error reading {args.config}: {exc}", file=sys.stderr)
            return 2
    else:
        # Pure-CLI mode: log_file and output_dir required
        if not args.log_file or not args.output_dir:
            print("Error: Without TOML config, --log-file and --output-dir "
                  "are required.\nFor a template TOML:\n"
                  "  modenanalyse-2fe2s --write-template run.toml",
                  file=sys.stderr)
            return 2
        cfg = Config(
            log_file   = args.log_file,
            output_dir = args.output_dir,
        )

    # CLI overrides (applied only if explicitly set)
    if args.log_file is not None:    cfg.log_file = args.log_file
    if args.output_dir is not None:  cfg.output_dir = args.output_dir
    if args.pdb_file is not None:    cfg.pdb_file = args.pdb_file
    if args.pdb_chain is not None:   cfg.pdb_chain = args.pdb_chain
    if args.freq_min is not None:    cfg.freq_min = args.freq_min
    if args.freq_max is not None:    cfg.freq_max = args.freq_max
    if args.freq_window:             cfg.freq_windows = args.freq_window
    if args.temp_k is not None:      cfg.temp_k = args.temp_k
    if args.cluster_index is not None: cfg.cluster_index = args.cluster_index
    if args.no_nis:
        import warnings as _w
        _w.warn("--no-nis is deprecated and ignored. NIS spectra are "
                "not produced by this tool.",
                DeprecationWarning, stacklevel=2)
    if args.no_scsd:   cfg.analyze_scsd = False
    if args.no_pcet:   cfg.pcet_enabled = False
    if args.no_cache:  cfg.use_cache    = False

    return int(run_analysis(cfg) or 0)


if __name__ == "__main__":
    sys.exit(main())
