# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Part of modenanalyse_2fe2s -- see LICENSE in repository root.

# -*- coding: utf-8 -*-
"""

config.py
======================
Configuration for the [2Fe-2S] normal-mode analysis.

Only this file (or rather: an instance of :class:`Config` constructed
in user code or loaded from TOML) needs to be customized. All other
modules consume the configuration unchanged.

Programmatic usage::

    from modenanalyse_2fe2s import Config, run_analysis
    cfg = Config(log_file=..., output_dir=..., temp_k=5.0)
    run_analysis(cfg)

TOML-based usage::

    modenanalyse-2fe2s run.toml
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class Config:
    """Complete configuration for an analysis run.

    Parameters
    ----------
    log_file : str
        Path to the Gaussian ``.log`` or ORCA ``.hess`` file.
    pdb_file : str
        Path to the PDB file for coordination detection and SS analysis.
    output_dir : str
        Output directory. Empty = directory of the .log file.
    freq_min : float or None
        Lower frequency cutoff in cm^-1. ``None`` means no filter.
        In multi-window mode (``freq_windows`` set): ignored;
        all modes are analyzed.
    freq_max : float or None
        Upper frequency cutoff in cm^-1. ``None`` means no filter.
        In multi-window mode: ignored.
    freq_windows : list of (float, float) or None
        Explicit frequency windows for multi-window operation.
        Each window (lo, hi) creates a subfolder
        ``"{lo:.0f}-{hi:.0f}_cm-1/"`` with separate Excel files.
        Per window: separate reorg aggregates, modulation spectra,
        and embedding. An overall run across all modes is also
        placed in the parent ``output_dir``.

        Example::

            freq_windows = [(0, 100), (100, 300), (300, 500), (500, 700)]

        ``None`` (default) = single window, legacy behavior.
    temp_k : float or None
        Temperature in Kelvin for thermal amplitude.
        ``None`` activates classical mode (amplitude = 1.0).
    pdb_chain : str
        Chain ID for ATOM records (e.g., ``"A"``). Empty = all chains.
        HETATM records (FES clusters) are always read.
    fe_coord_cutoff_n : float
        Maximum Fe-N distance in Angstrom for coordination detection.
    fe_coord_cutoff_s : float
        Maximum Fe-S distance in Angstrom for coordination detection.
    fe_coord_cutoff_o : float
        Maximum Fe-O distance in Angstrom for coordination detection.
    include_hn_vibration : bool
        Output H-N bond modulation for protonated histidines.
    interp_step : float
        Grid spacing for interpolated Excel in cm^-1.
    analyze_ss : bool
        Secondary-structure analysis (requires PDB).
    analyze_scsd : bool
        SCSD decomposition (requires ``pip install scsdpy``).
    ss_chain : str
        Chain for SS analysis. Empty = all chains.
    fe_s_cutoff : float
        Maximum Fe-S distance within the cluster in Angstrom.
    fe_fe_cutoff : float
        Maximum Fe-Fe distance within the cluster in Angstrom.
    include_hydrogen : bool
        Include H atoms in the eigenvector analysis.
    sigma_eigvec : float
        Uncertainty per eigenvector element (half-LSB).
    sigma_coord : float
        Uncertainty per coordinate in Angstrom.
    amplitude_threshold : float
        Minimum amplitude for OOP calculation.
    amplitude : float
        Fallback amplitude in classical mode (temp_k=None).
    coord_match_tol : float
        Tolerance for PDB-Gaussian atom matching in Angstrom.
    scan_chunk_mb : int
        Chunk size when scanning the .log file in MB.
    show_errors_in_excel : bool
        Output sigma columns in Excel.
    interp_boundary_mode : str
        Boundary treatment: ``'context'``, ``'zero'``, or ``'nearest'``.
    interp_edge_extend : float
        Grid extension beyond freq_min/freq_max in cm^-1.
    interp_context_cm1 : float
        Width of context-mode region beyond the window in cm^-1.
    use_cache : bool
        Enable byte-offset cache.
    strict_cluster : bool
        At ``True``: abort if cluster geometry outside literature values.
    umap_n_neighbors : int or None
        Override automatic n_neighbors for UMAP.

    Notes
    -----
    Coordinating amino acids are detected automatically from the PDB
    geometry. A ``residues`` parameter is no longer required.
    """

    # -- Input files -----------------------------------------------------
    # log_file: path to DFT frequency calculation. Supported formats:
    #   - Gaussian-16  .log  (with or without ``freq=hpmodes``)
    #   - ORCA          .hess (auto-detected)
    log_file:   str = r""
    pdb_file:   str = r""
    output_dir: str = r""

    # -- Frequency range -------------------------------------------------
    freq_min: Optional[float] = None
    freq_max: Optional[float] = None

    # -- Multi-window ----------------------------------------------------
    # Explicit frequency windows for separate subfolder outputs.
    # None = single window (legacy behavior, freq_min/max apply).
    # Set = analyze all modes, NIS global, per-window export.
    freq_windows: Optional[List[Tuple[float, float]]] = None

    # -- Thermal --------------------------------------------------------
    temp_k: Optional[float] = None

    # -- PDB chain ------------------------------------------------------
    pdb_chain: str = "A"

    # -- Fe-ligand detection --------------------------------------------
    fe_coord_cutoff_n: float = 2.50
    fe_coord_cutoff_s: float = 2.70
    fe_coord_cutoff_o: float = 2.50
    include_hn_vibration: bool = True

    # -- OOP/INP classification (Hardening v3.0, point 6) ---------------
    #
    # OOP% = fraction of squared displacement along the cluster normal;
    # INP% = 100 - OOP%.
    #
    # ``mode_type`` (binary, 3 levels): compatible with older analyses
    # and Excel coloring. Threshold is :data:`mode_type_threshold` (default
    # 60%).
    #
    # ``mode_type_detail`` (7 levels): finer partition for system
    # comparison, with three symmetric thresholds
    # :data:`mode_type_detail_thresholds`. default (60, 75, 90):
    #   OOP% >= 90       -> "Pure OOP"
    #   OOP% in [75, 90) -> "Strong OOP"
    #   OOP% in [60, 75) -> "Majority OOP"
    #   OOP% in [40, 60) -> "Mixed"
    #   OOP% in [25, 40) -> "Majority INP"
    #   OOP% in [10, 25) -> "Strong INP"
    #   OOP% < 10        -> "Pure INP"
    mode_type_threshold: float = 60.0
    mode_type_detail_thresholds: Tuple[float, float, float] = (60.0, 75.0, 90.0)

    # -- Universal significance convention (v3.5) -----------------------
    # Classification of each quantity with sigma:
    #   |X|/sigma <= low      -> "trivial"      (red in compact table)
    #   low < |X|/sigma <= high -> "significant" (yellow)
    #   |X|/sigma > high      -> "high"         (green)
    # Thresholds are also used for difference significance
    # (e.g. bend_oop vs. bend_inp), with sigma_diff = sqrt(s1^2 + s2^2).
    significance_threshold_low:  float = 1.0
    significance_threshold_high: float = 3.0

    # -- Output mode (v3.5) --------------------------------------------
    # analysis_compact: compact table with cell-wise color coding
    #                   (standard for daily work)
    # analysis_full:    additional *_full sheet per main sheet with
    #                   all values + sigmas + explicit significance
    #                   class columns (for publication / SI)
    analysis_compact: bool = True
    analysis_full:    bool = False

    # -- Secondary structure + SCSD ------------------------------------
    analyze_ss:           bool  = True
    analyze_scsd:         bool  = True
    ss_chain:             str   = ""

    # -- Cluster geometry ----------------------------------------------
    fe_s_cutoff:          float = 3.0
    fe_fe_cutoff:         float = 3.5
    # -- Multi-cluster selection (Hardening v3.0 #10) -----------------
    # In multi-iron systems (e.g., glutaredoxin dimer, domain constructs
    # with two [2Fe-2S] clusters), ``find_all_clusters`` finds all
    # clusters and sorts them by Fe-Fe distance (closest cluster first).
    #
    # Zwei Modi:
    # (1) ``cluster_index`` (Single-Cluster-Modus, Default): selects einen
    #     einzelnen Cluster for the Analyse. ``cluster_index = 0`` is das
    #     engste clusters, kompatibel with altem Verhalten.
    # (2) ``analyze_all_clusters = True`` (v1.0.0+): the Tool laeuft
    #     automatically per Cluster durch, jethe cluster bekommt einen
    #     Subordner ``cluster_<N>/`` below ``output_dir``. Ein
    #     Multi-Cluster-Summary-Sheet listet per Cluster the Status- und
    #     Reorg-Total-values. ``cluster_index`` is in doing so ignoriert.
    cluster_index:        int   = 0
    analyze_all_clusters: bool  = False

    # ── PCET / ET-Multi-Feature-Score (Hardening v3.0 #8) ─────────────────
    #
    # Zwei getrennte Scores fuer:
    #   PCET: Proton-coupled electron transfer (histidine-N-H modulierte
    #         H-Bond-Geometrie). Aktiv only for His-koordinierten Clustern.
    #         For rein Cys-koordinierten Systemen NaN (S-ligands koennen
    #         not protoniert/deprotoniert werden).
    #   ET:   Electron transfer (Cluster-Modes, Fe-Fe-Modulation).
    #         Aktiv for alle [2Fe-2S]-Cluster.
    #
    # Beide Scores are geometrische Mittel from drei Komponenten:
    # Lokalisierung x Cluster-Kopplung x tanh(geometrische Modulation/d_0).
    pcet_enabled:                 bool        = True
    #: Maximaler H...acceptor-Abstand for H-Bond-acceptor-Suche [A].
    pcet_hbond_cutoff_a:          float       = 4.0
    #: Gauss-center for acceptor-Gewichtung (typicaler H-Bond-Abstand).
    pcet_acceptor_r0_a:           float       = 2.8
    #: Gauss-Stddev for acceptor-Gewichtung.
    pcet_acceptor_sigma_a:        float       = 0.4
    
    # ── v3.7 Marcus-Hush Reorg-Spektren ───────────────────────────────────
    #: Gauss-Verbreiterung for Modulations-Spektren M_X(omega) in cm-1.
    #: 5 cm-1 corresponds to typicaler DFT-Frequenz-Unschaerfe + NRVS-Aufloesung.
    reorg_spectrum_sigma_cm1:     float       = 5.0
    #: Frequenz-Raster-Schrittweite for Modulations-Spektren in cm-1.
    #: 0.5 cm-1 gives feine Banden for vertretbarer file-Groesse.
    reorg_spectrum_step_cm1:      float       = 0.5

    # ── Eigenvektor-Analyse ───────────────────────────────────────────────
    include_hydrogen:     bool  = True
    sigma_eigvec:         float = 5e-4
    sigma_coord:          float = 1e-3
    amplitude_threshold:  float = 0.001
    amplitude:            float = 1.0
    coord_match_tol:      float = 1.0
    scan_chunk_mb:        int   = 16
    show_errors_in_excel: bool  = True

    # ── Interpolation ─────────────────────────────────────────────────────
    interp_step:          float = 0.05
    interp_boundary_mode: str   = "context"
    interp_edge_extend:   float = 0.5
    interp_context_cm1:   float = 5.0

    # ── Cache + Clustering ────────────────────────────────────────────────
    use_cache:            bool           = True
    strict_cluster:       bool           = False
    umap_n_neighbors:     Optional[int]  = None
    logname_suffix:       str            = ""
    #: Embedding-PNGs (UMAP-Streudiagramme) als Bilddateien
    #: schreiben. default True. If ``umap-learn`` not installiert ist,
    #: fehlen UMAP-Bilder (eine Warnung in the REPORT weist darauf hin).
    #: Auf False setzen, um matplotlib-Setup-Zeit to sparen.
    export_embedding_plots: bool         = True

    # ── Hilfsmethoden ─────────────────────────────────────────────────────

    @property
    def coord_cutoffs(self) -> Dict[str, float]:
        """Fe-Koordinations-Cutoffs after element in Angstrom."""
        return {
            "N": self.fe_coord_cutoff_n,
            "S": self.fe_coord_cutoff_s,
            "O": self.fe_coord_cutoff_o,
        }

    @property
    def fe_coord_cutoffs(self) -> Dict[str, float]:
        """Alias for coord_cutoffs (Abwaertskompatibilitaet with geometry.py)."""
        return self.coord_cutoffs

    def get_windows(self) -> List[Tuple[float, float]]:
        """Returns the effektiven frequency window .

        In multi-Fenster-Modus (``freq_windows`` gesetzt) are die
        expliziten Window zurueckgiven.  Im Einzelfenster-Modus
        (``freq_windows=None``) is ein einzelnes Window aus
        ``freq_min``/``freq_max`` gebildet -- identical with the alten
        Verhalten.

        Returns
        -------
        list of (float, float)
            frequency window als (lo, hi)-Tupel.  ``hi=inf`` bedeutet
            no oberes Limit (corresponds to ``freq_max=None``).
        """
        if self.freq_windows:
            return list(self.freq_windows)
        # Einzelfenster: altes Verhalten
        lo = self.freq_min if self.freq_min is not None else 0.0
        hi = self.freq_max if self.freq_max is not None else float("inf")
        return [(lo, hi)]

    def validate(self) -> List[str]:
        """Checks alle Konfigurationswerte and gibt Fehlermeldungen ."""
        errors = []
        if not self.log_file:
            errors.append("log_file is leer.")
        elif not os.path.isfile(self.log_file):
            errors.append(f"log_file not gefunden: {self.log_file}")
        if self.pdb_file and not os.path.isfile(self.pdb_file):
            errors.append(f"pdb_file not gefunden: {self.pdb_file}")
        if (self.freq_min is not None and self.freq_max is not None
                and self.freq_min >= self.freq_max):
            errors.append(
                f"freq_min ({self.freq_min}) muss kleiner als "
                f"freq_max ({self.freq_max}) sein.")
        # freq_windows validation
        if self.freq_windows is not None:
            if not isinstance(self.freq_windows, (list, tuple)):
                errors.append("freq_windows muss a List of (lo, hi)-Tupeln sein.")
            else:
                for i, w in enumerate(self.freq_windows):
                    if (not isinstance(w, (list, tuple)) or len(w) != 2
                            or not all(isinstance(x, (int, float)) for x in w)):
                        errors.append(
                            f"freq_windows[{i}] muss ein (float, float)-Tuple sein.")
                    elif w[0] >= w[1]:
                        errors.append(
                            f"freq_windows[{i}]: lo ({w[0]}) muss < hi ({w[1]}) sein.")
        if self.temp_k is not None and self.temp_k <= 0:
            errors.append(f"temp_k muss positiv sein ({self.temp_k}).")
        if self.interp_step <= 0:
            errors.append(f"interp_step muss positiv sein ({self.interp_step}).")
        if self.interp_boundary_mode not in ("context", "zero", "nearest"):
            errors.append(
                f"interp_boundary_mode muss 'context', 'zero' or 'nearest' sein "
                f"(ist: '{self.interp_boundary_mode}').")
        if self.scan_chunk_mb <= 0:
            errors.append(f"scan_chunk_mb muss positiv sein ({self.scan_chunk_mb}).")
        if self.sigma_eigvec < 0:
            errors.append(f"sigma_eigvec darf not negativ sein ({self.sigma_eigvec}).")
        if self.sigma_coord < 0:
            errors.append(f"sigma_coord darf not negativ sein ({self.sigma_coord}).")
        if self.coord_match_tol <= 0:
            errors.append(f"coord_match_tol muss positiv sein ({self.coord_match_tol}).")
        if self.fe_fe_cutoff <= 0 or self.fe_s_cutoff <= 0:
            errors.append("fe_fe_cutoff and fe_s_cutoff muessen positiv sein.")
        if self.cluster_index < 0:
            errors.append(
                f"cluster_index muss >= 0 sein ({self.cluster_index}).")
        # PCET/Reorg-Modulations-validation (v3.7)
        if self.pcet_hbond_cutoff_a <= 0:
            errors.append(
                f"pcet_hbond_cutoff_a muss > 0 sein ({self.pcet_hbond_cutoff_a}).")
        for f in ("pcet_acceptor_r0_a", "pcet_acceptor_sigma_a",
                  "reorg_spectrum_sigma_cm1", "reorg_spectrum_step_cm1"):
            v = getattr(self, f, None)
            if v is not None and v <= 0:
                errors.append(f"{f} muss > 0 sein ({v}).")
        if any(getattr(self, f) <= 0 for f in
               ("fe_coord_cutoff_n", "fe_coord_cutoff_s", "fe_coord_cutoff_o")):
            errors.append(
                "fe_coord_cutoff_n/s/o muessen positiv sein (Abstandsgrenzen in A).")
        if self.interp_edge_extend < 0:
            errors.append(
                f"interp_edge_extend darf not negativ sein ({self.interp_edge_extend}).")
        if self.interp_context_cm1 < 0:
            errors.append(
                f"interp_context_cm1 darf not negativ sein ({self.interp_context_cm1}).")
        if not isinstance(self.strict_cluster, bool):
            errors.append("strict_cluster muss True or False sein.")
        if not isinstance(self.show_errors_in_excel, bool):
            errors.append("show_errors_in_excel muss True or False sein.")
        if not isinstance(self.analyze_ss, bool):
            errors.append("analyze_ss muss True or False sein.")
        if not isinstance(self.analyze_scsd, bool):
            errors.append("analyze_scsd muss True or False sein.")
        if not isinstance(self.include_hn_vibration, bool):
            errors.append("include_hn_vibration muss True or False sein.")
        if self.amplitude_threshold < 0:
            errors.append(f"amplitude_threshold darf not negativ sein.")
        if self.umap_n_neighbors is not None and self.umap_n_neighbors < 2:
            errors.append(
                f"umap_n_neighbors muss None or >= 2 sein ({self.umap_n_neighbors}).")
        if not isinstance(self.use_cache, bool):
            errors.append("use_cache muss True or False sein.")
        if not isinstance(self.include_hydrogen, bool):
            errors.append("include_hydrogen muss True or False sein.")
        # OOP/INP classifications-Schwellen
        if not 50.0 < self.mode_type_threshold < 100.0:
            errors.append(
                f"mode_type_threshold muss in (50, 100) sein "
                f"({self.mode_type_threshold}).")
        thr = self.mode_type_detail_thresholds
        if (len(thr) != 3
                or not all(50.0 < t < 100.0 for t in thr)
                or not (thr[0] < thr[1] < thr[2])):
            errors.append(
                f"mode_type_detail_thresholds muss a streng aufsteigende "
                f"Sequenz (low, mid, high) with valuesn in (50, 100) sein "
                f"({thr}).")
        return errors

    def freq_label(self) -> str:
        """Returns a kurzes Label for the frequency range .

        Wird als Unterordner-Name used, e.g. ``"0-100_cm-1"``.
        Sind weder ``freq_min`` still ``freq_max`` gesetzt, is ein
        leerer String zurueckgiven (kein Unterordner).

        Returns
        -------
        str
            Label in the Format ``"{lo:.0f}-{hi:.0f}_cm-1"`` or ``""``.
        """
        if self.freq_min is None and self.freq_max is None:
            return ""
        lo = f"{self.freq_min:.0f}" if self.freq_min is not None else "0"
        hi = f"{self.freq_max:.0f}" if self.freq_max is not None else "max"
        return f"{lo}-{hi}_cm-1"

    def outdir(self) -> str:
        """Returns the Ausgabeverzeichnis .

        If ``freq_min`` or ``freq_max`` gesetzt sind, is ein
        Unterordner with the frequency range-Label angeplaces.

        Returns
        -------
        str
            Ausgabeverzeichnis inklusive Frequenz-Unterordner.
        """
        base = self.output_dir if self.output_dir \
               else os.path.dirname(os.path.abspath(self.log_file))
        label = self.freq_label()
        if label:
            path = os.path.join(base, label)
            os.makedirs(path, exist_ok=True)
            return path
        return base

    def outname(self, suffix: str) -> str:
        """Creates a completeen Ausgabepfad.

        Parameters
        ----------
        suffix : str
            fileendung or -suffix, e.g. ``"_analysis.xlsx"``.

        Returns
        -------
        str
            Vollstaendiger path ``<outdir>/<basename><suffix>``.
        """
        base = os.path.splitext(os.path.basename(self.log_file))[0]
        return os.path.join(self.outdir(), base + suffix)

    # ── TOML-Support (v3.2) ───────────────────────────────────────────────

    @classmethod
    def from_toml(cls, path: str) -> "Config":
        """Creates a ``Config`` from a TOML-file.

        Die TOML-file kann beliebig in Sektionen gegliedert sein
        (e.g. ``[input]``, ``[freq]``, ``[nis]``, ``[pcet]``); alle
        key are flach in the ``Config`` uebernommen.
        Unbekannte key fuehren to a ``ValueError``.

        Parameters
        ----------
        path : str
            Path to the TOML-file.

        Returns
        -------
        Config

        Beispiel
        --------
        A minimale ``run.toml``::

            [input]
            log_file   = "D:/Daten/dimer.log"
            pdb_file   = "D:/Daten/dimer.pdb"
            output_dir = "D:/Daten/results"

            [freq]
            freq_max = 800.0
            freq_windows = [[0, 100], [100, 300], [300, 500], [500, 700]]

            [thermo]
            temp_k = 40.0
        """
        try:
            import tomllib                  # Python 3.11+
        except ImportError as exc:          # pragma: no cover
            raise RuntimeError(
                "Config.from_toml required Python >= 3.11 (tomllib). "
                "Aktuell: " + str(exc)) from exc

        with open(path, "rb") as f:
            raw = tomllib.load(f)

        # Sektionen flach machen: Top-Level Keys + alle Sub-Tables
        flat: Dict[str, object] = {}
        for k, v in raw.items():
            if isinstance(v, dict):
                # Sub-Tabelle (e.g. [freq])
                for kk, vv in v.items():
                    if kk in flat:
                        raise ValueError(
                            f"TOML: key '{kk}' doppelt defined "
                            f"(Top-Level and in Sektion '{k}').")
                    flat[kk] = vv
            else:
                flat[k] = v

        # freq_windows in TOML als List[List[float]] -> List[Tuple[float,float]]
        if "freq_windows" in flat and flat["freq_windows"] is not None:
            try:
                flat["freq_windows"] = [
                    (float(lo), float(hi)) for (lo, hi) in flat["freq_windows"]
                ]
            except Exception as exc:
                raise ValueError(
                    f"TOML: 'freq_windows' muss a List of [lo, hi]-Paaren "
                    f"be. Error: {exc}") from exc

        # mode_type_detail_thresholds: TOML list -> Tuple
        if "mode_type_detail_thresholds" in flat and \
                isinstance(flat["mode_type_detail_thresholds"], list):
            flat["mode_type_detail_thresholds"] = tuple(
                float(x) for x in flat["mode_type_detail_thresholds"])

        # Legacy-Parameter from frueheren Versionen ignorieren with Warnung.
        # These key haben in the aktuellen Version keinen Effekt mehr;
        # sie are akzeptiert, so that alte TOMLs not crashen.
        _legacy_silent_drop = {
            # v3.6 PCET-Score-Parameter (entfernt in v3.7.0)
            "pcet_dr_normalization_a", "pcet_strong_threshold",
            "pcet_moderate_threshold", "pcet_band_edges_cm1",
            "cpet_lambda0_cm1", "pt_lambda0_cm1", "et_lambda0_cm1",
            # legacy embedding parameters (removed when t-SNE was dropped)
            "tsne_perplexity",
            # legacy NIS fields (NIS calculations are not part of this tool;
            # use a dedicated NIS package on the same log file).
            "analyze_nis", "nis_fwhm_gauss", "nis_fwhm_lorentz",
            "nis_lineshape", "nis_n_points", "nis_freq_min", "nis_freq_max",
            "nis_phonon_order", "nis_n_theta", "nis_split_elastic",
            "nis_inmemory_max_atom_modes", "nis_resolution_cm1",
        }
        import warnings as _w
        for _k in list(flat.keys()):
            if _k in _legacy_silent_drop:
                _w.warn(
                    f"TOML: key '{_k}' is obsolete and ignored. "
                    f"NIS spectra are not produced by this tool.",
                    UserWarning, stacklevel=2)
                flat.pop(_k)

        # validation against the Felder the Datenklasse
        import dataclasses as _dc
        valid_fields = {f.name for f in _dc.fields(cls)}
        unknown = set(flat) - valid_fields
        if unknown:
            raise ValueError(
                f"TOML: Unbekannte key: {sorted(unknown)}. "
                f"Allowed are: {sorted(valid_fields)}.")

        return cls(**flat)

    def to_toml(self, path: str, *, group: bool = True) -> None:
        """Writes the ``Config`` als TOML-file.

        By default grouped into sections
        (``[input]``, ``[freq]``, ``[thermo]``, ``[cluster]``,
        ``[ligand]``, ``[oop]``, ``[scsd]``, ``[pcet]``,
        ``[embedding]``, ``[numerics]``). Mit ``group=False`` als
        flache file.

        Parameters
        ----------
        path : str
            Zielpfad for the TOML-file.
        group : bool
            If True (Default), Felder in thematische Sektionen
            gruppieren.
        """
        import dataclasses as _dc
        flat = {f.name: getattr(self, f.name) for f in _dc.fields(self)}

        # Gruppen-Mapping (musst only Top-Level-Felder enthalten)
        groups = {
            "input":     ["log_file", "output_dir", "pdb_file", "pdb_chain",
                          "logname_suffix"],
            "freq":      ["freq_min", "freq_max", "freq_windows"],
            "thermo":    ["temp_k", "amplitude"],
            "cluster":   ["fe_s_cutoff", "fe_fe_cutoff", "cluster_index",
                          "analyze_all_clusters", "strict_cluster"],
            "ligand":    ["fe_coord_cutoff_n", "fe_coord_cutoff_s",
                          "fe_coord_cutoff_o", "include_hn_vibration"],
            "oop":       ["mode_type_threshold",
                          "mode_type_detail_thresholds"],
            "output":    ["analysis_compact", "analysis_full",
                          "significance_threshold_low",
                          "significance_threshold_high"],
            "scsd":      ["analyze_scsd", "analyze_ss", "ss_chain"],
            "pcet":      ["pcet_enabled", "pcet_hbond_cutoff_a",
                          "pcet_acceptor_r0_a",
                          "pcet_acceptor_sigma_a",
                          "reorg_spectrum_sigma_cm1",
                          "reorg_spectrum_step_cm1"],
            "embedding": ["umap_n_neighbors", "export_embedding_plots"],
            "numerics":  ["include_hydrogen", "sigma_eigvec", "sigma_coord",
                          "amplitude_threshold", "coord_match_tol",
                          "interp_step", "interp_boundary_mode",
                          "interp_edge_extend", "interp_context_cm1",
                          "scan_chunk_mb", "use_cache",
                          "show_errors_in_excel"],
        }

        def _emit_value(v) -> str:
            """Serialisiert einen Python-value after TOML."""
            if v is None:
                # TOML hat no None; wir schreiben the Feld nicht
                return ""
            if isinstance(v, bool):
                return "true" if v else "false"
            if isinstance(v, (int, float)):
                return repr(v)
            if isinstance(v, str):
                # Backslashes escapen + in Anfuehrungszeichen
                escaped = v.replace("\\", "\\\\").replace('"', '\\"')
                return f'"{escaped}"'
            if isinstance(v, (list, tuple)):
                # Listen rekursiv
                inner = ", ".join(_emit_value(x) for x in v if _emit_value(x) != "")
                return f"[{inner}]"
            raise TypeError(f"No TOML serializer for type {type(v)}: {v!r}")

        lines: List[str] = [
            "# modenanalyse_2fe2s -- configuration",
            f"# Created with to_toml() from version {__version__}",
            "",
        ]
        used = set()
        if group:
            for section, keys in groups.items():
                section_lines = []
                for k in keys:
                    if k not in flat:
                        continue
                    used.add(k)
                    s = _emit_value(flat[k])
                    if s == "":
                        # None weglassen
                        continue
                    section_lines.append(f"{k} = {s}")
                if section_lines:
                    lines.append(f"[{section}]")
                    lines.extend(section_lines)
                    lines.append("")
            # Restliche Felder (falls neue dazukommen, the hier nicht
            # gemappt sind) ans Ende
            rest = [k for k in flat if k not in used]
            if rest:
                lines.append("[other]")
                for k in rest:
                    s = _emit_value(flat[k])
                    if s != "":
                        lines.append(f"{k} = {s}")
        else:
            for k, v in flat.items():
                s = _emit_value(v)
                if s != "":
                    lines.append(f"{k} = {s}")

        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")


__version__ = "1.0.2"
# v1.0.1: Documentation cleanup release. No code changes — analysis
# pipeline is identical to v1.0.0 and produces numerically identical
# results. See CHANGELOG.md for the full list of documentation fixes.
# v1.0.0 (2026-05-07): First public release.
