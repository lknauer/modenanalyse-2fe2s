# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Part of modenanalyse_2fe2s -- see LICENSE in repository root.

# -*- coding: utf-8 -*-
"""
runner.py
=====================
Entry point of the [2Fe-2S] normal-mode analysis.

The program is designed exclusively for [2Fe-2S] clusters (Rieske,
ferredoxin, mitoNEET, etc.). Other cluster types ([4Fe-4S],
[3Fe-4S], P-cluster) would require a different cluster-detection
logic and reference geometries.

Multi-cluster systems (dimers, multi-domain constructs) are
supported as a first-class feature via ``Config.analyze_all_clusters``.

Running
-------
::

    modenanalyse-2fe2s run.toml

or programmatically::

    from modenanalyse_2fe2s import Config, run_analysis
    cfg = Config(log_file=..., output_dir=..., temp_k=5.0)
    run_analysis(cfg)

Output files
------------
``_analysis.xlsx``
    Main analysis: mode_analysis, group amplitudes, Fe-ligand amplitudes,
    His H-N, equilibrium geometry, SCSD, reorganization energies,
    B-factors, info.
``_analysis_SS.xlsx``
    Secondary-structure amplitudes (if PDB available).
``_analysis_Embeddings.xlsx``
    UMAP coordinates + HDBSCAN clusters + C-alpha amplitudes.
``_analysis_interp{step}.xlsx``
    Core analysis on uniform frequency grid (step = interp_step).
    Symmetric boundary treatment: context modes left and right of the
    window.
``_analysis_SS_interp{step}.xlsx``
    Secondary-structure amplitudes interpolated (if PDB and analyze_ss=True).
``_embedding_*.png``
    Embedding figure (UMAP).
``_REPORT.txt``
    Run report: configuration, coordination, mode distribution, warnings.

Frequency subfolders
--------------------
With ``freq_min`` / ``freq_max`` set, outputs are placed in a subfolder
(e.g. ``100-300_cm-1/``). Multiple windows can thus be collected under
a common ``output_dir``.

Cache
-----
With ``use_cache=True`` (default), the log scan is cached after the
first run. Subsequent runs (e.g. different frequency windows on the
same file) skip the full scan completely and are significantly faster.
"""
from __future__ import annotations
import os, sys, time, warnings
from typing import Dict, List
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")
warnings.filterwarnings("ignore", category=UserWarning,    message=".*n_jobs.*")

import numpy as np

from .config import Config
from .logio import (
    RunLog, scan_log, read_std_orient, read_all_meta, get_eigvec, parse_pdb,
    load_scan_cache, save_scan_cache,
    check_hp_std_frequency_consistency,
    _ELEM,
)
from .geometry import (
    find_cluster, find_all_clusters, cluster_normal, compute_dist_ref,
    kabsch_align, find_coordinating_residues,
    build_ss_center_map, get_calpha_centers,
    detect_ss_dssp, detect_ss_phipsi,
    CoordInfo,
)
from .core import (
    analyze_mode_with_fallback, analyze_all_ss,
    compute_scsd_for_mode_full, _get_scsd_model,
    reset_warning_state,
)
from .embedding import (
    build_feature_matrix, compute_embeddings, compute_ss_umap_cluster,
    compute_ca_umap_cluster,
)
from .export import (

    ExportPayload, export_all,
    # Einzelfunktionen bleiben importiert for Abwaertskompatibilitaet
    export_main_excel, export_ss_excel, export_embedding_excel,
    export_interpolated_excel, export_ss_interp_excel,
    export_embedding_plots,
)


# ===========================================================================
def _build_ca_data(ca_pre_list, amps_by_mode, res_list):
    """Assembliert Calpha-Amplituden-Matrix for a Ergebnisliste.

    Parameters
    ----------
    ca_pre_list : list of (center, resname)
    amps_by_mode : dict of {mode_num: list of float}
    res_list : list of dict

    Returns
    -------
    tuple (centers, resnames, matrix) or None
    """
    if not ca_pre_list or not amps_by_mode:
        return None
    n_ca = len(ca_pre_list)
    ca_matrix = np.zeros((n_ca, len(res_list)))
    for mi, r in enumerate(res_list):
        amps = amps_by_mode.get(r["number"])
        if amps is not None and len(amps) == n_ca:
            ca_matrix[:, mi] = amps
    return ([c for c, _ in ca_pre_list],
            [rn for _, rn in ca_pre_list],
            ca_matrix)


def _make_synthetic_zero_mode(freq_cm1: float,
                                coord_info,
                                fe_c: List[int],
                                s_c:  List[int]) -> Dict:
    """Build a synthetic null-mode result dict (v1.0.4 FUND 13).

    Produces a Result-dict that matches the structure of a real mode
    analysed by ``analyze_mode_with_fallback`` but has all observable
    fields set to zero. Used as an explicit decay anchor for the
    interpolated pDOS when there is no real mode above ``freq_max``
    (i.e., the DFT spectrum ends inside the analysis range).

    Why this exists. Without an upper anchor, ``np.interp`` would simply
    return its ``right=0.0`` boundary value, producing a step at
    ``freq_max``. With a synthetic zero anchor at
    ``freq_max + interp_context_cm1``, the interpolation *decays
    smoothly* toward zero across the buffer zone, which is the
    physically correct behaviour (no modes there = no contribution).

    Parameters
    ----------
    freq_cm1 : float
        Frequency at which to place the synthetic zero. Conventionally
        ``cfg.freq_max + cfg.interp_context_cm1``.
    coord_info : CoordInfo
        Used only to mirror the group/ligand dict structure of a real
        mode (all entries get zero values).
    fe_c, s_c : list of int
        Cluster atom centres; mirrored into ``cl_com``/``cl_exp``/etc.

    Returns
    -------
    dict
        A result dict shaped exactly like a real mode but with all
        observables zero. The ``_evg``/``_centers``/``_c2l`` keys are
        omitted so the B-factor accumulator and SS-analysis loop in the
        runner skip this entry (both check for ``r.get("_evg")``).
    """
    # Zero per-group/per-ligand/per-his_hn dicts (same keys as a real
    # mode, all values zero).
    zero_group = {k: 0. for k in
        ("oop", "inp", "angle", "torsion", "total",
         "s_oop", "s_inp", "s_angle", "s_tors")}
    groups = {gn: dict(zero_group) for gn in coord_info.group_map.keys()}

    zero_fe_lig = {
        "stretch": 0., "bend": 0., "bend_inp": 0., "bend_oop": 0.,
        "s_stretch": 0., "s_bend": 0., "s_bend_inp": 0., "s_bend_oop": 0.,
        "bend_significance": "trivial", "lig_element": "?"}
    fe_lig = {l.res_label: dict(zero_fe_lig, lig_element=l.lig_element)
                for l in coord_info.ligands}

    his_hn = {}  # no protonated-His H-N stretching for synthetic mode

    return {
        "number":    -1,       # sentinel: not a real Gaussian mode index
        "freq":      float(freq_cm1),
        "red_mass":  0.,
        "frc_const": 0.,
        "sym":       "A",
        "precision": "synthetic",
        "mode_type": "synthetic_zero",
        "mode_type_detail": "synthetic_zero",
        # Ring 2: ligand sphere (all zero)
        "lig_oop_pct":   0., "lig_inp_pct":   0., "lig_d": 0.,
        "s_lig_oop":     0., "s_lig_d":       0.,
        # Ring 3: secondary sphere (all zero)
        "second_oop_pct": 0., "second_inp_pct": 0., "second_d": 0.,
        "s_second_oop":   0., "s_second_d":     0.,
        # Ring 1: cluster core (all zero)
        "kern_oop":  0., "kern_inp":  0., "kern_d":    0.,
        "s_kern_oop": 0., "s_kern_d":  0.,
        "cl_com":    np.zeros(3), "cl_exp":  0., "cl_rot": 0.,
        "groups":    groups,
        "fe_lig":    fe_lig,
        "his_hn":    his_hn,
        "u_rms":     0.,
        "kern_primary":   "synthetic", "kern_secondary": "synthetic",
        "kern_scores":    {},
        "kern_loc":       0.,
        # Reorganization channels: empty dict means
        # compute_total_reorganization / compute_modulation_spectra skip
        # this mode automatically (their loops check truthiness).
        "reorg_per_mode": {},
        "reorg_subchannels": [],
        "pts_ref":   None,
        "pts_dist":  None,
        # NOTE: NO "_evg" / "_centers" / "_c2l" keys. The B-factor and
        # SS-analysis loops in the runner explicitly check
        # `r.get("_evg") is not None` and skip otherwise -- so this
        # synthetic mode contributes 0 to B-factors and 0 to SS by
        # construction, without us having to fabricate eigenvectors.
    }


# ===========================================================================

def _run_analysis_single(cfg: "Config") -> int:
    """Single-Cluster-Pipeline: analysiert genau einen Cluster (cfg.cluster_index).

    Interner Helfer; the oeffentliche API is :func:`run_analysis`, die
    for ``cfg.analyze_all_clusters=True`` einen Multi-Cluster-Loop um
    diese Funktion herumplaces.

    Parameters
    ----------
    cfg : Config
        Vollstaendige Konfiguration. Typical creation via
        ``from modenanalyse_2fe2s import Config`` and Anpassen der
        Felder ``log_file``, ``pdb_file``, ``output_dir``, etc.

    Returns
    -------
    int
        Exit-Code: 0 = erfolgreich, 1 = Error (REPORT enthaelt Details).

    Ablauf
    ------
    1. Config validieren and Gaussian-Log scannen.
    2. PDB laden, Cluster erkennen, Kabsch-Alignment.
    3. Jede Mode analysieren (HP-Eigenvectors beforezugt, Fallback on Standard).
    4. Marcus-Hush-Reorg-Aggregate berechnen (Lambda_X total + Spektren).
    5. Embeddings and Cluster-Analyse.
    6. Excel fileen + REPORT.txt write.

    Beispiel
    --------
    >>> from modenanalyse_2fe2s import Config, run_analysis
    >>> cfg = Config(
    ...     log_file   = r"D:\\Daten\\dimer.log",
    ...     pdb_file   = r"D:\\Daten\\dimer.pdb",
    ...     output_dir = r"D:\\Daten\\results",
    ...     temp_k     = 40.0,
    ... )
    >>> rc = run_analysis(cfg)
    """
    runlog = RunLog(cfg)
    t_start = time.time()

    def _abort(msg: str = "") -> None:
        """Writes BEFUND and beendet the Programm with Exit-Code 1.

        Stellt sicher dass also for fruehen Errorabbruechen
        a BEFUND-file for the Errordiagnose present ist.
        """
        if msg:
            runlog.error(msg)
        try:
            # Fallback: if cfg.outname() fehlschlägt (e.g. log_file leer),
            # schreibe BEFUND ins aktuelle directory
            try:
                befund_err = cfg.outname("_REPORT.txt")
            except Exception:
                import tempfile as _tmp, os as _os
                befund_err = _os.path.join(
                    _os.getcwd(), "modenanalyse_FEHLER_REPORT.txt")
            runlog.write_befund(befund_err)
            print(f"  REPORT written: {befund_err}")
        except Exception as _be:
            print(f"  [WARNING] REPORT could not be written: {_be}")
        sys.exit(1)

    # ── Konfiguration validieren ──────────────────────────────────────────
    errs = cfg.validate()
    if errs:
        print("\n[ERROR] Configuration invalid:")
        for e in errs: print(f"  - {e}")
        _abort()

    # Reset per-run warning state (so warnings.warn() in analyze_mode is
    # emitted at least once per run even when run_analysis() is called
    # repeatedly from the same Python session, e.g. in a notebook).
    reset_warning_state()

    os.makedirs(cfg.outdir(), exist_ok=True)
    base = os.path.splitext(os.path.basename(cfg.log_file))[0]

    print("\n" + "="*70)
    print("  runner.py  -  [2Fe-2S] normal-mode analysis")
    print(f"  Temperature: {cfg.temp_k} K")
    print("="*70)
    print(f"  File:   {os.path.basename(cfg.log_file)}")
    print(f"  PDB:    {os.path.basename(cfg.pdb_file) if cfg.pdb_file else '(none)'}")
    filt = (f"{'-' if cfg.freq_min is None else cfg.freq_min} - "
           f"{'-' if cfg.freq_max is None else cfg.freq_max} cm-1")
    print(f"  Filter: {filt}  |  sigma_ev={cfg.sigma_eigvec:.1e}  sigma_coord={cfg.sigma_coord:.1e}")

    # Module status
    for mod in ["openpyxl","scipy","sklearn","umap","hdbscan"]:
        try:
            __import__(mod if mod != "sklearn" else "sklearn")
            runlog.module_status[mod] = True
        except ImportError:
            runlog.module_status[mod] = False
    try:
        from scsd.scsd import scsd_model as _
        runlog.module_status["scsdpy"] = True
    except ImportError:
        runlog.module_status["scsdpy"] = False

    # ── Format-Erkennung (v3.3) ──────────────────────────────────────────
    # Gaussian .log or ORCA .hess? For ORCA are Phase 1 (Scan) und
    # Phase 3 (Block metadata) through direkte Extraktion from dem
    # ParseResult ersetzt -- ein .hess hat no streaming-relevanten
    # blocks and no HP/Std-Trennung.
    from .orca_io import (is_orca_input, load_orca_hess,
                          parseresult_to_atoms, parseresult_to_blocks,
                          get_eigvec_orca)
    _is_orca = is_orca_input(cfg.log_file)
    _orca_pr = None  # OrcaHessResult; only beplaces if ORCA

    if _is_orca:
        print("\n  Phase 1: reading ORCA .hess...")
        runlog.info(f"Format detection: ORCA .hess ({cfg.log_file})")
        try:
            _orca_pr = load_orca_hess(cfg.log_file)
        except Exception as exc:
            runlog.error(f"Error reading ORCA .hess: {exc}")
            _abort()
        runlog.info(f"ORCA: {_orca_pr.n_atoms} atoms, {_orca_pr.n_modes} modes")
        print(f"    {_orca_pr.n_atoms} atoms, {_orca_pr.n_modes} modes loaded")

        # Dummy-Offsets, so that the Konsistenzpruefungen unten gluecklich sind
        so_off = [0]
        nc_off = []
        fr_off = [0]
    else:
        # ── Phase 1: Scan (oder Cache laden) ─────────────────────────────────
        print("\n  Phase 1: scan...")
        _cached_early = load_scan_cache(cfg.log_file) if cfg.use_cache else None
        if _cached_early is not None:
            _, _, _, so_off, nc_off, fr_off = _cached_early
            print(f"    Cache: scan offsets loaded (full scan skipped)")
            runlog.info("Cache used: full scan skipped (byte offsets from cache)")
        else:
            so_off, nc_off, fr_off = scan_log(cfg.log_file, cfg)
            print(f"    SO: {len(so_off)}  NC: {len(nc_off)}  "
                  f"Freq groups: {len(fr_off)}")

        if not so_off: runlog.error("Keine 'Standard orientation' found."); _abort()
        if not fr_off: runlog.error("No frequency blocks found.");        _abort()

    # ── Phase 2: Geometrie ────────────────────────────────────────────────
    print("\n  Phase 2: geometry...")

    if _is_orca:
        atoms_h, idx_map_h = parseresult_to_atoms(_orca_pr, include_hydrogen=True)
        atoms,   idx_map   = parseresult_to_atoms(_orca_pr, include_hydrogen=False)
    else:
        # Atome MIT H (fuer Koordinationserkennung)
        atoms_h, idx_map_h = read_std_orient(cfg.log_file, so_off[-1],
                                              include_hydrogen=True)
        # Atome OHNE H (fuer alle weiteren Berechnungen)
        atoms, idx_map = read_std_orient(cfg.log_file, so_off[-1],
                                          include_hydrogen=False)
    _atom_msg = f"{len(atoms)} heavy atoms ({len(atoms_h)} total with H)"
    runlog.info(_atom_msg)
    print(f"    {_atom_msg}")

    # Cluster finden
    runlog.info("Phase 2: cluster detection")
    print("  Cluster detection...")

    # Pre-Scan: alle Cluster auflisten (Hardening v3.0 #10), so that the User
    # for Multi-Cluster-Systemen the Auswahl nachvollziehen kann.
    try:
        _all_clusters = find_all_clusters(atoms, cfg)
    except Exception as e:
        _all_clusters = []
        runlog.warn(f"Cluster-pre-scan failed: {e}")

    if len(_all_clusters) > 1:
        runlog.info(f"Multi-cluster system detected: {len(_all_clusters)} "
                    f"[2Fe-2S] clusters found.")
        for _i, (_fe, _s, _g) in enumerate(_all_clusters):
            _marker = "  -> ANALYSIERT" if _i == cfg.cluster_index else ""
            runlog.info(
                f"  Cluster #{_i}: Fe={_fe}, S={_s}, "
                f"Fe-Fe={_g['fe_fe']:.3f} A, "
                f"Fe-S={_g['fe_s_min']:.3f}-{_g['fe_s_max']:.3f} A"
                f"{_marker}")
        print(f"    {len(_all_clusters)} clusters found, "
              f"analyzing #{cfg.cluster_index}")
        # v1.0.0: note about multi-cluster mode
        runlog.info(
            f"IMPORTANT: only cluster #{cfg.cluster_index} is analyzed. "
            f"For all {len(_all_clusters)} clusters automatically: "
            f"set analyze_all_clusters = true in the config (TOML) "
            f"or Config(analyze_all_clusters=True, ...) in Python.")
        print(f"    Note: only cluster #{cfg.cluster_index} is "
              f"analysiert. For all Cluster: analyze_all_clusters=true")
    elif len(_all_clusters) == 1:
        _g = _all_clusters[0][2]
        runlog.info(
            f"One cluster found: Fe-Fe={_g['fe_fe']:.3f} A, "
            f"Fe-S={_g['fe_s_min']:.3f}-{_g['fe_s_max']:.3f} A")

    # Warnungen from find_cluster/_check_cluster_geometry abfangen und
    # in RunLog umleiten (Hardening #11: einheitliches Warning-Konzept)
    import warnings as _warn_mod
    with _warn_mod.catch_warnings(record=True) as _caught_warnings:
        _warn_mod.simplefilter("always")
        try:
            fe_c, s_c = find_cluster(atoms, cfg)
        except ValueError as e:
            runlog.error(f"Cluster detection: {e}"); _abort()
    for _w in _caught_warnings:
        runlog.warn(f"Cluster-Geometry: {_w.message}")

    # cluster normal + folding residual (B0, v3.5) -- the Normale is die
    # Bezugsachse for alle nachfolgenden OOP/INP-Berechnungen. Das Residual
    # quantifiziert the Abweichung von the Best-Fit-Ebene (sigma_3 the SVD).
    normal   = cluster_normal(atoms, idx_map, fe_c, s_c)
    dist_ref = compute_dist_ref(atoms, idx_map, fe_c, s_c, cfg)
    # Abstände in the Fortschrittslog ausgeben (aus compute_dist_ref ausgelagert)
    for _dk, (_dv, _ds) in dist_ref.items():
        print(f"    {_dk}: {_dv:.6f} \u00b1 {_ds:.1e} A")

    # cluster normal berichten (Komponenten + SVD-folding residual)
    try:
        _pts_cl = np.array([
            [atoms[idx_map[c]]["x"], atoms[idx_map[c]]["y"], atoms[idx_map[c]]["z"]]
            for c in (fe_c + s_c) if c in idx_map
        ])
        if _pts_cl.shape[0] >= 3:
            _ctr = _pts_cl.mean(0)
            _, _sv, _ = np.linalg.svd(_pts_cl - _ctr)
            _residual = float(_sv[-1]) / np.sqrt(max(_pts_cl.shape[0] - 1, 1))
        else:
            _residual = float("nan")
    except Exception:
        _residual = float("nan")
    print(f"    Cluster normal n_hat: ({normal[0]:+.4f}, {normal[1]:+.4f}, {normal[2]:+.4f})")
    print(f"    Folding residual:     {_residual:.4f} A  (RMS distance of 4 cluster atoms from best-fit plane)")

    # Cluster geometry strukturiert in the RunLog store (Hardening #3 + B0 v3.5)
    # Important: hier *erweitern*, not ueberschreiben — sonst gehen die
    # cluster normal and the folding residual von oben verloren.
    runlog.geometry["cluster_normal"]            = [float(normal[0]),
                                                    float(normal[1]),
                                                    float(normal[2])]
    runlog.geometry["cluster_plane_residual_a"]  = _residual
    runlog.geometry["n_heavy"] = len(atoms)
    runlog.geometry["n_total"] = len(atoms_h)
    runlog.geometry["fe_fe"]   = dist_ref.get("Fe-Fe", (0., 0.))[0]
    for k,(v,_) in dist_ref.items():
        runlog.geometry[k.lower().replace("-","_")] = v

    # document geometric plausibility in REPORT
    _fe_fe = dist_ref.get("Fe-Fe", (0.,0.))[0]
    if _fe_fe > 0:
        _fe_s_vals = [v for k,(v,_) in dist_ref.items() if "S" in k]
        _fe_s_avg  = float(np.mean(_fe_s_vals)) if _fe_s_vals else 0.
        # Fe-S Einzelwerte strukturiert speichern (Hardening #3a)
        runlog.geometry["fe_s_distances"] = {
            k: v for k,(v,_) in dist_ref.items() if "S" in k}
        runlog.geometry["fe_s_mean"] = _fe_s_avg
        # Fe-S-Fe Winkel berechnen if alle vier Atome present (Hardening #3b)
        if len(fe_c) >= 2 and len(s_c) >= 2:
            _ang_list = []
            for _sc in s_c:
                if _sc in idx_map and fe_c[0] in idx_map and fe_c[1] in idx_map:
                    _ps = np.array([atoms[idx_map[_sc]]["x"],
                                    atoms[idx_map[_sc]]["y"],
                                    atoms[idx_map[_sc]]["z"]])
                    _pf1 = np.array([atoms[idx_map[fe_c[0]]]["x"],
                                     atoms[idx_map[fe_c[0]]]["y"],
                                     atoms[idx_map[fe_c[0]]]["z"]])
                    _pf2 = np.array([atoms[idx_map[fe_c[1]]]["x"],
                                     atoms[idx_map[fe_c[1]]]["y"],
                                     atoms[idx_map[fe_c[1]]]["z"]])
                    _v1 = _pf1 - _ps; _v2 = _pf2 - _ps
                    _n1 = float(np.linalg.norm(_v1)); _n2 = float(np.linalg.norm(_v2))
                    if _n1 > 0 and _n2 > 0:
                        _ang = float(np.degrees(np.arccos(
                            np.clip(np.dot(_v1,_v2)/(_n1*_n2), -1., 1.))))
                        _ang_list.append(_ang)
            if _ang_list:
                runlog.geometry["fe_s_fe_angles_deg"] = _ang_list
                runlog.geometry["fe_s_fe_mean_deg"]   = float(np.mean(_ang_list))
        runlog.info(
            f"Cluster geometry: Fe-Fe={_fe_fe:.3f} A, "
            f"Fe-S(mean)={_fe_s_avg:.3f} A"
            + (f", Fe-S-Fe={runlog.geometry.get('fe_s_fe_mean_deg',0.):.1f} deg"
               if "fe_s_fe_mean_deg" in runlog.geometry else "")
            + f" - {'OK' if 2.40 < _fe_fe < 3.10 else 'UNUSUAL (outside 2.40-3.10 A)'}")
        if not (2.40 < _fe_fe < 3.10):
            msg = (f"Fe-Fe = {_fe_fe:.3f} A outside typicalem [2Fe-2S]-Bereich "
                   f"(2.40-3.10 A). Cluster detection pruefen.")
            if cfg.strict_cluster:
                runlog.error(msg + " [strict_cluster=True: Abbruch]")
                _abort()
            else:
                runlog.warn(msg)

    # ── PDB laden + Koordination erkennen ────────────────────────────────
    print("  PDB + coordination...")
    pdb_data  = None
    coord_info = CoordInfo(ligands=[], group_map={}, pdb_to_center={},
                            his_ligand_labels=[])
    ss_elements: list = []
    ss_center_map: dict = {}

    if cfg.pdb_file and os.path.isfile(cfg.pdb_file):
        pdb_data = parse_pdb(cfg.pdb_file, chain_filter=cfg.pdb_chain)

        R, t, kabsch_rmsd, kabsch_ok = kabsch_align(
            pdb_data, atoms, idx_map, fe_c, s_c, runlog)
        runlog.group_match["kabsch_rmsd"] = kabsch_rmsd

        if kabsch_ok:
            print(f"    Kabsch alignment OK  (RMSD≈{kabsch_rmsd:.3f} A)")
            coord_info = find_coordinating_residues(
                pdb_data, atoms_h, idx_map_h,
                atoms, idx_map,
                fe_c, s_c, R, t, cfg, runlog)
        else:
            runlog.warn("Kabsch failed -- automatic coordination detection skipped.")

        # REPORT: Koordinations-Zusammenfassung
        runlog.coord_summary = {
            "ligands": [
                {
                    "res_label":      lig.res_label,
                    "element":        lig.lig_element,
                    "aname":          lig.lig_aname,
                    "fe_idx":         lig.fe_idx,
                    "bond_len":       lig.bond_len,
                    "his_protonated": lig.his_protonated,
                    "hn_len":         lig.hn_len,
                    "hn_n_type":      (
                        "ND1" if lig.his_hn_center == lig.lig_center
                        else "NE2"
                        if lig.his_hn_center is not None else None),
                    "hn_via":         (
                        "PDB" if lig.his_protonated
                        and lig.h_center is not None
                        and lig.his_hn_center is not None
                        else "Gaussian-Fallback"),
                }
                for lig in coord_info.ligands
            ],
            "cluster_distances": {k: v for k, (v, _) in dist_ref.items()},
        }

        # Secondary structure
        if cfg.analyze_ss:
            ss_elements = pdb_data["ss_elements"]
            ch = cfg.ss_chain or ""
            ss_elements = [e for e in ss_elements
                           if not ch or e.get("chain", "") == ch]

            # Fallback 1: DSSP (H-Brücken, Kabsch & Sander 1983 [11])
            # Fallback 2: phi/psi-Winkel (nur if DSSP leer)
            _ss_auto   = False
            _ss_method = ""
            if not ss_elements:
                ss_elements = detect_ss_dssp(
                    pdb_data, chain_filter=ch or "A")
                if ss_elements:
                    _ss_auto   = True
                    _ss_method = "DSSP (H-bonds)"
                else:
                    import warnings as _w_phipsi
                    with _w_phipsi.catch_warnings(record=True) as _phipsi_warns:
                        _w_phipsi.simplefilter("always")
                        ss_elements = detect_ss_phipsi(
                            pdb_data, chain_filter=ch or "A")
                    for _pw in _phipsi_warns:
                        runlog.warn(f"SS phi/psi: {_pw.message}")
                    if ss_elements:
                        _ss_auto   = True
                        _ss_method = "phi/psi-Winkel"

                if _ss_auto:
                    runlog.warn(
                        f"PDB contains no HELIX/SHEET records. "
                        f"{len(ss_elements)} SS elements automatically "
                        f"per {_ss_method} erkannt "
                        f"(Kabsch & Sander 1983 [11]). "
                        f"Ergebnisse können von HELIX/SHEET-basierten "
                        f"valuesn leicht abweichen.")
                    print(f"    [INFO] SS auto-detected ({_ss_method}): "
                          f"{len(ss_elements)} elements")
                else:
                    runlog.warn(
                        "PDB contains no HELIX/SHEET records. "
                        "Weder DSSP still phi/psi-Erkennung fanden "
                        "SS elements. SS-Analyse skipped.")
                    print("    [WARNING] No SS elements detected "
                          "(no HELIX/SHEET, no DSSP, no phi/psi).")

            if ss_elements:
                ss_center_map = build_ss_center_map(
                    ss_elements, pdb_data, coord_info.pdb_to_center, ch)
                n_ss = sum(len(v) for v in ss_center_map.values())
                src_lbl = (f"auto, {_ss_method} [11]"
                           if _ss_auto else "HELIX/SHEET")
                print(f"    {len(ss_elements)} SS elements "
                      f"[{src_lbl}] ({n_ss} atoms assigned)")
                if n_ss == 0:
                    runlog.warn(
                        "SS analysis: SS elements detected, but 0 atoms "
                        "assigned (pdb_to_center empty). "
                        "SS-Export skipped.")
                    ss_center_map = {}
    else:
        if cfg.pdb_file:
            runlog.warn(f"PDB not found: {cfg.pdb_file}")
        print("    No groups -- core analysis only.")

    # PCET/ET-Multi-Feature-Score-Vorberechnung (v3.7: only noch
    # Geometrie-Erkennung, no Score-Berechnung mehr) ──────────────
    # Wir bauen the PCET- and ET-Atom listn einmal before the modenschleife,
    # so that jede Mode the Modulationen without erneute Geometriesuche
    # berechnen kann. Die eigentliche Marcus-Hush-Reorganisations-
    # Berechnung occurs in core.py via reorganization.py.
    if cfg.pcet_enabled:
        from .pcet_et import build_pcet_info, build_et_info
        coord_info.atoms_h    = atoms_h
        coord_info.idx_map_h  = idx_map_h
        coord_info.pcet_info  = build_pcet_info(
            coord_info, atoms_h, idx_map_h, fe_c, s_c, cfg)
        coord_info.et_info    = build_et_info(coord_info, fe_c, s_c)
        if coord_info.pcet_info.enabled:
            runlog.info(
                f"Reorg modulations active "
                f"({coord_info.pcet_info.diagnose}). "
                f"Computes per mode: dr_X and lambda_X = (1/2) mu omega^2 dr^2 "
                f"for all bonding channels X = FeFe, FeN, FeS, NH, HA. "
                f"system aggregation: Total-Reorg, Modulations-Spektren, "
                f"kumulative Lambda(omega).")
            print(f"    Reorg modulations active: "
                  f"{coord_info.pcet_info.n_his} His ligands, "
                  f"{sum(len(a) for a in coord_info.pcet_info.acceptor_centers_per_h)} "
                  f"H-Bond-Paare")
        else:
            # Bugfix v1.0.4 (post-release Apd1 audit): pcet_info.n_his
            # counts only PROTONATED His ligands. To distinguish "no
            # His at all" from "His present but all deprotonated" we
            # additionally inspect coord_info.ligands directly. The
            # previous message "no His ligands at cluster" was actively
            # misleading for deprot systems with His ligands.
            n_his_total = sum(
                1 for l in coord_info.ligands
                if l.res_name.upper() == "HIS" and l.lig_element == "N")
            if n_his_total == 0:
                runlog.info(
                    "PCET reorg: NOT active (no His ligands at cluster). "
                    "PCET is physically impossible; S ligands (Cys) "
                    "cannot be protonated/deprotonated under physiological "
                    "conditions. NH and HA channels are skipped.")
                print("    PCET reorg: NOT active (no His)")
            elif coord_info.pcet_info.n_his == 0:
                # His present but none protonated (deprot system)
                runlog.info(
                    f"PCET reorg: NOT active ({n_his_total} His ligand(s) "
                    f"found but none protonated). NH and HA channels are "
                    f"skipped. To enable PCET analysis, supply a Gaussian "
                    f"log with at least one protonated His or check that "
                    f"the H atoms are present in PDB/Gaussian inputs.")
                print(f"    PCET reorg: NOT active ({n_his_total} His "
                      f"present but all deprotonated)")
            else:
                runlog.warn(
                    f"PCET reorg: His ligands found, but no "
                    f"H-bond acceptors within {cfg.pcet_hbond_cutoff_a} A "
                    f"range. HA channel omitted.")
        runlog.info(
            f"ET reorg ligands: {len(coord_info.et_info.ligand_centers)} "
            f"atoms detected (cluster + ligand sphere for FeFe/FeN/FeS "
            f"reorg channels).")

    # SCSD-Modell  (Symmetry-Coordinate Structural Decomposition,
    # Kingsbury & Senge, Chem. Sci. 2024, 15, 13638; see Modulkopf
    # in core.py).
    scsd_model = None
    if cfg.analyze_scsd:
        print("  SCSD model (Kingsbury method, canonical D2h reference, "
              "Fe-Fe=2.73A, Fe-S=2.20A)...")
        # SCSD-UserWarnings in RunLog umleiten
        import warnings as _scsd_warn
        with _scsd_warn.catch_warnings(record=True) as _scsd_caught:
            _scsd_warn.simplefilter("always")
            scsd_model = _get_scsd_model()
        for _sw in _scsd_caught:
            runlog.warn(f"SCSD: {_sw.message}")
        if scsd_model:
            print("    scsdpy OK  (reference: 2Fe2S_canonical, D2h)")
            runlog.info(
                "SCSD method: Kingsbury & Senge, Chem. Sci. 15, 13638 (2024). "
                "Canonical D2h reference: Fe-Fe = 2.73 A, Fe-S = 2.20 A "
                "(means from Rieske/ferredoxin crystal structures). "
                "Axis convention: x = Fe-Fe, y = S-S, z = cluster normal. "
                "The displacement is orthogonally decomposed into D2h irreps; "
                "values are directly comparable between different structures."
            )
        else:
            runlog.warn("scsdpy not installed or SCSD model failed.")

    # ── Phase 3: Block metadata (mit Cache) ────────────────────────────────
    print("\n  Phase 3: Block metadata...")
    if _is_orca:
        # ORCA: ein einziger synthetischer Pseudo-Block with allen Modes
        all_blocks, best_block, cand_map = parseresult_to_blocks(_orca_pr)
        _cache_msg = (f"ORCA: {len(best_block)} modes in synthetic "
                      f"pseudo-block")
        runlog.info(_cache_msg)
        print(f"    {_cache_msg}")
    else:
        _cached = _cached_early if _cached_early is not None else (
            load_scan_cache(cfg.log_file) if cfg.use_cache else None)
        if _cached is not None:
            all_blocks, best_block, cand_map, so_off, nc_off, fr_off = _cached
            _cache_msg = f"Cache geladen: {len(best_block)} Moden"
            runlog.info(_cache_msg)
            print(f"    {_cache_msg}")
        else:
            all_blocks, best_block, cand_map = read_all_meta(
                cfg.log_file, nc_off, fr_off)
            if cfg.use_cache:
                save_scan_cache(cfg.log_file, all_blocks, best_block,
                                cand_map, so_off, nc_off, fr_off,
                                runlog=runlog)
                print(f"    Cache gespeichert")
                runlog.info("Cache written (full scan gecacht)")
    if not best_block:
        runlog.error("No modes found."); _abort()

    # HP/Standard-Frequenz-Konsistenzcheck (Hardening v3.0, Punkt 11):
    # If beide Block-Typen present sind, muessen the Frequenzen exakt
    # uebereinstimmen — bis on the Druckpraezision the Standard-Ausgabe
    # (~ 0.005 cm-1). Abweichungen darvia deuten on Parser-Versatz,
    # korrupte Logdatei or verkettete freq-Jobs hin and are im
    # REPORT festgehalten.
    # For ORCA gibt es no HP/Std-Trennung -- Check skip.
    if _is_orca:
        runlog.hp_std_check = {
            "any_hp": True, "any_std": False, "checked_modes": 0,
            "n_outliers": 0, "max_dev_cm1": 0.0, "mean_dev_cm1": 0.0,
            "outliers": [], "skipped_reason": "ORCA hat no HP/Std-Trennung",
        }
        runlog.info("HP/Std consistency check: not applicable for ORCA input.")
    else:
        _hpstd = check_hp_std_frequency_consistency(cand_map, tolerance_cm1=0.01)
        runlog.hp_std_check = _hpstd
        if not _hpstd["any_hp"] or not _hpstd["any_std"]:
            # Nur ein Block-Typ present — Vergleich not possible
            if _hpstd["any_hp"]:
                runlog.info("HP/Std consistency check: only HP blocks available, "
                            "no comparison with default blocks possible.")
            elif _hpstd["any_std"]:
                runlog.info("HP/Std-Konsistenzcheck: only Standard-blocks "
                            "available, no comparison with HP possible. "
                            "Empfehlung: 'freq=hpmodes' verwenden for hoehere "
                            "Eigenvektor-Praezision.")
        elif _hpstd["checked_modes"] == 0:
            runlog.warn("HP/Std-Konsistenzcheck: no Mode in beiden Block-"
                        "types available -- check block indexing.")
        elif _hpstd["n_outliers"] == 0:
            runlog.info(
                f"HP/Std consistency check OK: {_hpstd['checked_modes']} modes "
                f"compared, max |HP-std| = {_hpstd['max_dev_cm1']:.4f} cm-1, "
                f"mean deviation = {_hpstd['mean_dev_cm1']:.4f} cm-1 "
                f"(tolerance 0.01 cm-1).")
        else:
            # Diskrepanzen: warnen + Top-5 listen
            runlog.warn(
                f"HP/Std-Konsistenzcheck: {_hpstd['n_outliers']}/"
                f"{_hpstd['checked_modes']} modes with |HP-std| > 0.01 cm-1. "
                f"Maximale Abweichung: {_hpstd['max_dev_cm1']:.4f} cm-1. "
                f"Mittlere Abweichung: {_hpstd['mean_dev_cm1']:.4f} cm-1. "
                f"Possible causes: parser offset, corrupted log file, "
                f"verkettete freq-Jobs.")
            for (mn, f_hp, f_std, dev) in _hpstd["outliers"]:
                runlog.warn(
                    f"  Mode {mn}: HP = {f_hp:.4f} cm-1, std = {f_std:.4f} cm-1, "
                    f"|dev| = {dev:.4f} cm-1")

    all_freqs = [f for bi in all_blocks for f in bi.freqs]
    if all_freqs:
        print(f"    All modes: {len(best_block)}, "
              f"f = {min(all_freqs):.2f}-{max(all_freqs):.2f} cm-1")

    # Filter + Modenauswahl
    # In multi-Fenster-Modus (freq_windows gesetzt) are alle positiven
    # modes analyzed; the Fenster-Filterung is done later per Fenster.
    _multi_window = bool(cfg.freq_windows)
    selected = []
    for bi in all_blocks:
        for col,(mn,freq) in enumerate(zip(bi.mode_nums, bi.freqs)):
            if mn not in best_block or best_block[mn] is not bi: continue
            if freq <= 0:                                   # B18: imaginaere Moden
                runlog.add_parse_failure(mn, freq,
                    f"imaginary/zero ({freq:.3f} cm-1), skipped")
                continue
            # Hardening v3.1: Frequenzfilter gilt jetzt immer (auch im
            # Multi-window-Modus). For NRVS-Messbereich bis ~800 cm-1 und
            # verschwindend kleinen Fe-projizierten Beicontribute daruber sind
            # Modes outside of the analysierten Bereichs without Mehrwert; sie
            # auszurechnen kostet only Laufzeit (typ. -70% for freq_max=500).
            if cfg.freq_min is not None and freq < cfg.freq_min: continue
            if cfg.freq_max is not None and freq > cfg.freq_max: continue
            if _multi_window:
                if not any(lo <= freq <= hi for (lo, hi) in cfg.freq_windows):
                    continue
            selected.append((bi, col, mn, freq))
    selected.sort(key=lambda x: x[3])
    print(f"    After filter: {len(selected)} Moden")
    if not selected:
        runlog.error("No modes after filter."); _abort()

    # Modenanzahl-Pruefung
    n_heavy = len(atoms); n_found = len(best_block)
    n_exp   = max(0, 3*n_heavy - 6)
    runlog.mode_stats = {"n_found": n_found, "n_filtered": len(selected), "n_hp": 0}
    if n_found == n_exp:
        print(f"    Modenanzahl: {n_found} OK (3×{n_heavy}−6)")
    elif n_found > 1.05*n_exp:
        print(f"    [INFO] {n_found} modes (HP+standard, deduplication active)")
    else:
        msg = f"Nur {n_found}/{n_exp} modes found."
        print(f"    [WARNING] {msg}"); runlog.warn(msg)

    # Calpha-centers vorab
    ca_pre = []
    if pdb_data:
        ca_pre = get_calpha_centers(pdb_data, atoms, coord_info.pdb_to_center)

    # ── Phase 4: mode analysis ─────────────────────────────────────────────
    print(f"\n  Phase 4: analysis of {len(selected)} modes...")

    # ORCA: ParseResult in cfg injizieren, so that _get_eigvec_smart in core.py
    # the Eigenvectors from the RAM lesen kann (statt from the Gaussian-Log).
    if _is_orca:
        cfg._orca_parse_result = _orca_pr

    results = []; n_fail = 0
    b_accum = np.zeros(len(atoms))  # Debye-Waller-Akkumulator
    for i,(bi,col,mn,freq) in enumerate(selected):
        if (i+1)%100==0 or (i+1)==len(selected):
            _prog = f"Analysis: {i+1}/{len(selected)} modes ({freq:.1f} cm-1)"
            if (i+1)==len(selected): runlog.info(_prog)
            print(f"    [{i+1:4d}/{len(selected)}] {freq:.2f} cm-1 "
                  f"({'HP' if bi.is_hp else 'std'})    ", end="\r")
        try:
            cands = cand_map.get(mn, [bi])
            # Hardening v3.0: Warnungen from analyze_mode (z. B. fehlgeschlagene
            # H-N-Analyse) direkt in the RunLog umleiten, statt sie per
            # default-Warning-Filter aufs stderr to kippen. Konsistent mit
            # Cluster/NIS/SCSD-Phasen.
            import warnings as _w_am
            with _w_am.catch_warnings(record=True) as _am_caught:
                _w_am.simplefilter("always")
                r, fail_reason = analyze_mode_with_fallback(
                    cands, col, cfg.log_file, atoms, idx_map,
                    normal, coord_info, fe_c, s_c, cfg, mode_num=mn)
            for _wam in _am_caught:
                runlog.warn(f"{_wam.message}")

            if r is None:
                runlog.add_parse_failure(mn, freq, fail_reason or "unbekannt")
                n_fail += 1; continue

            # Calpha-Amplituden (inline, solange _evg verfuegbar)
            if ca_pre and r.get("_evg") is not None and r.get("_c2l") is not None:
                from .core import evg_sub_extern as _ext
                evg_r = r["_evg"]; c2l_r = r["_c2l"]
                ca_amps = []
                for c_ca, _ in ca_pre:
                    row = c2l_r.get(c_ca)
                    ca_amps.append(float(np.linalg.norm(evg_r[row]))
                                   if row is not None and row < evg_r.shape[0] else 0.)
                r["_ca_amps"] = ca_amps

            # SCSD
            if (scsd_model and r.get("pts_ref") is not None
                    and r.get("pts_dist") is not None):
                try:
                    r["scsd"] = compute_scsd_for_mode_full(
                        r["pts_ref"], r["pts_dist"], scsd_model, r["u_rms"],
                        sigma_coord=cfg.sigma_coord,
                        sigma_eigvec=cfg.sigma_eigvec)
                except Exception as _e:
                    r["scsd"] = {}
                    runlog.warn(f"SCSD Mode {mn} @ {freq:.2f} cm-1: "
                                f"{type(_e).__name__}: {_e}")

            r.pop("pts_ref", None); r.pop("pts_dist", None)

            # Secondary structure (Bugfix B2: korrekte Indizierung via c2l)
            if ss_elements and ss_center_map and r.get("_evg") is not None:
                evg_r = r["_evg"]; c2l_r = r.get("_c2l",{})
                try:
                    r["ss"] = analyze_all_ss(
                        evg_r, c2l_r, ss_center_map,
                        atoms, idx_map, ss_elements, r["u_rms"],
                        sigma_eigvec=cfg.sigma_eigvec)
                except Exception as _e_ss:
                    r["ss"] = {}
                    runlog.warn(f"SS Mode {mn} @ {freq:.2f} cm-1: "
                                f"{type(_e_ss).__name__}: {_e_ss}")

            # Debye-Waller-Faktor akkumulieren (vor _evg-Pop)
            # Saves also pro-Mode-contributions for the spaetere
            # fensterweise B-Faktor-Berechnung in the multi-window-Modus.
            if r.get("_evg") is not None:
                _evg_r = r["_evg"]; _c2l_r = r.get("_c2l", {})
                _b_contrib_r: dict = {}
                for _bc, _bi in _c2l_r.items():
                    _ba = idx_map.get(_bc)
                    if _ba is not None and _bi < _evg_r.shape[0]:
                        _contrib = float(np.sum(_evg_r[_bi] ** 2))
                        b_accum[_ba] += _contrib
                        _b_contrib_r[_ba] = _b_contrib_r.get(_ba, 0.) + _contrib
                r["_b_contribs"] = _b_contrib_r

            # Interne Arrays entfernen
            r.pop("_evg", None); r.pop("_centers", None); r.pop("_c2l", None)
            results.append(r)

        except Exception as exc:
            runlog.add_parse_failure(mn, freq, str(exc)); n_fail += 1

    n_hp = sum(1 for r in results if r["precision"]=="high")
    n_std = sum(1 for r in results if r["precision"]!="high")
    runlog.mode_stats["n_hp"] = n_hp
    runlog.mode_stats["n_std_fallback"] = n_std
    print(f"\n    {len(results)} modes analyzed "
          f"(HP: {n_hp}, standard fallback: {n_std}, "
          f"failed: {n_fail}).")
    if not results:
        runlog.error("All modes failed."); _abort()

    # Debye-Waller-Faktoren finalisieren
    # B_i = 8π^2 * Σ_l u_rms(l)^2 * |e_{i,l}|^2 / 3  [A^2]
    # Nur from analysierten modes in the frequency window.
    b_factors = 8. * np.pi**2 * b_accum / 3.

    # REPORT: Moden-Verteilung
    _mtypes = {}
    for r in results:
        t = r.get("mode_type", "?")
        _mtypes[t] = _mtypes.get(t, 0) + 1
    _oops = [r["lig_oop_pct"] for r in results if "lig_oop_pct" in r]
    runlog.results_summary = {
        "n_modes":    len(results),
        "mode_types": _mtypes,
        "freq_range": (min(r["freq"] for r in results),
                       max(r["freq"] for r in results)),
        "mean_lig_oop": float(np.mean(_oops)) if _oops else 0.,
        "n_scsd":     sum(1 for r in results if r.get("scsd")),
    }

    # ── Kontext-Moden: symmetrisch an beiden Raendern ────────────────────
    # Rechts (jenseits freq_max): verhindert abrupten Abfall at the oberen Rand
    # Links  (unterhalb freq_min): verhindert Nullsetzung at the unteren Rand
    #   → only if freq_min > erster Gesamtmode (echte modes liegen darunter)
    #
    # v1.0.4 (FUND 12): Wir laden Modes IM Fenster
    # [freq_max, freq_max+interp_context_cm1] als bisher PLUS — falls dieses
    # Fenster keine einzige Mode enthaelt — zusaetzlich die EINE
    # naechstliegende Mode oberhalb freq_max (analog links). Damit hat
    # np.interp() immer mindestens einen Anker jenseits the Grenze, sodass
    # the lineare Interpolation auch dann sauber bis freq_max laeuft, wenn
    # the naechste Mode weiter als interp_context_cm1 entfernt liegt.
    #
    # v1.0.4 (FUND 13, user-requested follow-up): Falls auch das nicht hilft,
    # weil ueberhaupt keine echte Mode oberhalb freq_max existiert (echtes
    # Spektrumsende), fuegen wir eine SYNTHETISCHE Null-Mode bei
    # freq_max + interp_context_cm1 ein. Das ist physikalisch der korrekte
    # Decay-Anker -- oberhalb des DFT-Spektrums existieren wirklich keine
    # Beitraege -- und macht das Verhalten explizit (frueher implizit ueber
    # np.interp's right=0.0).
    context_results = []
    if cfg.freq_max is not None and cfg.interp_context_cm1 > 0:
        ctx_hi = cfg.freq_max + cfg.interp_context_cm1
        # Alle Modes-im-Kontextfenster (nahe Kontextmodes)
        candidates_in_window = [(bi, col, mn, freq)
            for bi in all_blocks
            for col, (mn, freq) in enumerate(zip(bi.mode_nums, bi.freqs))
            if (mn in best_block and best_block[mn] is bi
                and freq > cfg.freq_max and freq <= ctx_hi)]
        # v1.0.4 FUND 12: wenn keine Modes im Kontextfenster, nimm die
        # naechste Mode oberhalb als "minimaler Anker"
        _use_synthetic_upper_zero = False
        if not candidates_in_window:
            all_above = sorted(
                [(bi, col, mn, freq)
                 for bi in all_blocks
                 for col, (mn, freq) in enumerate(zip(bi.mode_nums, bi.freqs))
                 if (mn in best_block and best_block[mn] is bi
                     and freq > cfg.freq_max)],
                key=lambda x: x[3])
            if all_above:
                candidates_in_window = [all_above[0]]
                runlog.info(
                    f"No context modes in [{cfg.freq_max:.1f}, "
                    f"{ctx_hi:.1f}] cm-1; using single anchor mode "
                    f"#{all_above[0][2]} @ {all_above[0][3]:.2f} cm-1 "
                    f"for upper-edge interpolation (v1.0.4 FUND 12).")
            else:
                # v1.0.4 FUND 13: keine echte Mode oberhalb freq_max.
                # Synthetische Null-Mode als Decay-Anker einfuegen.
                _use_synthetic_upper_zero = True
        selected_ctx = sorted(candidates_in_window, key=lambda x: x[3])
        ctx_fail_r = 0
        for bi, col, mn, freq in selected_ctx:
            try:
                cands = cand_map.get(mn, [bi])
                import warnings as _w_ctx
                with _w_ctx.catch_warnings(record=True) as _ctx_caught:
                    _w_ctx.simplefilter("always")
                    r, _ = analyze_mode_with_fallback(
                        cands, col, cfg.log_file, atoms, idx_map,
                        normal, coord_info, fe_c, s_c, cfg, mode_num=mn)
                for _wctx in _ctx_caught:
                    runlog.warn(f"{_wctx.message}")
                if r is not None:
                    context_results.append(r)
            except Exception as _e_ctx:
                ctx_fail_r += 1
                runlog.warn(f"Kontext-Mode rechts {mn} @ {freq:.2f} cm\u207b\xb9: "
                            f"{type(_e_ctx).__name__}: {_e_ctx}")
        # v1.0.4 FUND 13: Append synthetische Null-Mode falls noetig.
        # Wichtig: nur einfuegen, wenn weder echte Kontextmodes noch
        # FUND-12-Fallback-Modes geladen wurden (sonst doppelt anchored).
        if _use_synthetic_upper_zero and not context_results:
            _synth_freq = float(ctx_hi)
            _synth_result = _make_synthetic_zero_mode(
                freq_cm1 = _synth_freq,
                coord_info = coord_info,
                fe_c = fe_c, s_c = s_c)
            context_results.append(_synth_result)
            runlog.info(
                f"Synthetic zero anchor at {_synth_freq:.2f} cm-1 "
                f"(no real modes above freq_max={cfg.freq_max:.1f}). "
                f"Interpolated quantities decay to 0 above this point "
                f"(v1.0.4 FUND 13).")
        if context_results:
            runlog.info(f"{len(context_results)} context modes right "
                        f"({cfg.freq_max:.1f}\u2013{ctx_hi:.1f} cm\u207b\xb9)")
        if ctx_fail_r:
            runlog.warn(f"{ctx_fail_r} context modes right failed "
                        f"(Interpolation at the oberen Rand if applicable ungenau).")

    context_results_left = []
    if cfg.freq_min is not None and cfg.interp_context_cm1 > 0:
        # Returns es ueberhaupt modes unterhalb freq_min?
        all_freqs_below = [freq
                           for bi in all_blocks
                           for (mn, freq) in zip(bi.mode_nums, bi.freqs)
                           if mn in best_block and best_block[mn] is bi
                           and 0 < freq < cfg.freq_min]
        if all_freqs_below:
            ctx_lo = cfg.freq_min - cfg.interp_context_cm1
            candidates_in_window_l = [(bi, col, mn, freq)
                for bi in all_blocks
                for col, (mn, freq) in enumerate(zip(bi.mode_nums, bi.freqs))
                if (mn in best_block and best_block[mn] is bi
                    and freq < cfg.freq_min and freq >= ctx_lo)]
            # v1.0.4 FUND 12: wenn keine Modes im Kontextfenster, nimm die
            # naechste Mode unterhalb als "minimaler Anker"
            if not candidates_in_window_l:
                all_below = sorted(
                    [(bi, col, mn, freq)
                     for bi in all_blocks
                     for col, (mn, freq) in enumerate(zip(bi.mode_nums, bi.freqs))
                     if (mn in best_block and best_block[mn] is bi
                         and 0 < freq < cfg.freq_min)],
                    key=lambda x: -x[3])  # absteigend → naechste UNTER freq_min zuerst
                if all_below:
                    candidates_in_window_l = [all_below[0]]
                    runlog.info(
                        f"No context modes in [{ctx_lo:.1f}, "
                        f"{cfg.freq_min:.1f}] cm-1; using single anchor "
                        f"mode #{all_below[0][2]} @ "
                        f"{all_below[0][3]:.2f} cm-1 for lower-edge "
                        f"interpolation (v1.0.4 FUND 12).")
            selected_ctx_l = sorted(candidates_in_window_l, key=lambda x: x[3])
            ctx_fail_l = 0
            for bi, col, mn, freq in selected_ctx_l:
                try:
                    cands = cand_map.get(mn, [bi])
                    import warnings as _w_ctxl
                    with _w_ctxl.catch_warnings(record=True) as _ctxl_caught:
                        _w_ctxl.simplefilter("always")
                        r, _ = analyze_mode_with_fallback(
                            cands, col, cfg.log_file, atoms, idx_map,
                            normal, coord_info, fe_c, s_c, cfg, mode_num=mn)
                    for _wctxl in _ctxl_caught:
                        runlog.warn(f"{_wctxl.message}")
                    if r is not None:
                        context_results_left.append(r)
                except Exception as _e_ctx:
                    ctx_fail_l += 1
                    runlog.warn(f"Kontext-Mode links {mn} @ {freq:.2f} cm\u207b\xb9: "
                                f"{type(_e_ctx).__name__}: {_e_ctx}")
            if context_results_left:
                runlog.info(f"{len(context_results_left)} context modes left "
                            f"({ctx_lo:.1f}\u2013{cfg.freq_min:.1f} cm\u207b\xb9)")
            if ctx_fail_l:
                runlog.warn(f"{ctx_fail_l} context modes left failed "
                            f"(Interpolation at the unteren Rand if applicable ungenau).")

    # ── Calpha-Daten als Mode-Nummer-Dict speichern ──────────────────────────
    # In multi-Fenster-Modus is the Ca-Matrix spaeter per Window assembliert.
    # In the Einzelfenster-Modus is sie direkt hier gebaut.
    ca_amps_by_mode: dict = {}
    ca_data = None
    if ca_pre:
        runlog.info("Phase 4: C-alpha amplitudes")
        print("  Calpha-Amplituden...")
        for r in results:
            _amps = r.pop("_ca_amps", None)
            if _amps is not None:
                ca_amps_by_mode[r["number"]] = _amps
        if ca_amps_by_mode:
            print(f"    {len(ca_pre)} C-alpha atoms ({len(ca_amps_by_mode)} modes)")
        if not _multi_window:
            # Einzelfenster: Ca-Matrix direkt bauen (altes Verhalten)
            ca_data = _build_ca_data(ca_pre, ca_amps_by_mode, results)

    # ── Cluster-Info aufbauen (wird in allen Ausgabemodi required) ──────
    _cluster_info = []
    for lbl, c_list in [("Fe1", fe_c[:1]), ("Fe2", fe_c[1:2]),
                         ("S1",  s_c[:1]),  ("S2",  s_c[1:2])]:
        if not c_list: continue
        c = c_list[0]
        if c in idx_map:
            a = atoms[idx_map[c]]
            coords = (a["x"], a["y"], a["z"])
            elem = a.get("symbol", "?"); note = ""
        else:
            coords, elem, note = None, "?", "not in idx_map"
        _cluster_info.append((lbl, c, elem, coords, note))

    # ── v3.7: Marcus-Hush Total-Reorg + Modulations-Spektren ──────────────
    # Pro Lauf: aggregiere the Pro-Mode-Modulationen to System-Quantitaeten.
    # Das is einmal global gerechnet and in the Results als Diagnose entry
    # mitgiven (damit the Export-Mode sie without Re-Computing schreiben kann).
    def _compute_v37_aggregates(results_subset):
        """Computes the v3.7 System-Aggregate for a Mode list.
        
        Returns a Dict with keyn:
          'totals'             : compute_total_reorganization-Output
          'spectra_grid_cm1'   : Frequenz-Raster
          'modulation_spectra' : M_X(omega) per Kanal
          'co_modulation_spectra' : C_PCET, C_PT_FeN, C_ET_FeS
          'cumulative_reorg'   : Lambda_X(omega)
          'spectra_sigma_cm1'  : Verbreiterung
          'cumulative_uses_mode_mass' : True/False
        """
        from .reorganization import (
            compute_total_reorganization, compute_modulation_spectra,
            compute_co_modulation_spectra, compute_cumulative_reorganization,
        )
        per_mode_aggs = [r.get("reorg_per_mode") for r in results_subset
                         if r.get("reorg_per_mode")]
        if not per_mode_aggs:
            return None
        freqs = np.array([r["freq"] for r in results_subset
                          if r.get("reorg_per_mode")], dtype=float)
        totals = compute_total_reorganization(per_mode_aggs)
        
        # Sub-Channel-Totals: per einzelner bond (FeS_Cys207, FeS_Cluster_X_Y,
        # FeN_His255, FeN_His259, HA_His255_..., etc.) seine eigene Total-Reorg.
        # So that man sieht, welche bond at the meisten beicontributes.
        sub_totals: dict = {}
        for r in results_subset:
            sub_list = r.get("reorg_subchannels") or []
            for (name, parent, weight,
                 lam_pair, lam_mode, dr_sig) in sub_list:
                if name not in sub_totals:
                    sub_totals[name] = {
                        "parent_channel":         parent,
                        "weight":                 weight,
                        "lambda_total_pair_cm1":  0.0,
                        "lambda_total_mode_cm1":  0.0,
                        "n_modes_contributing":   0,
                    }
                d = sub_totals[name]
                # v3.7.4: konsistente Count -- jede Mode, the zu
                # IRGENDEINEM the beiden Lambdas (pair or mode) beicontributes,
                # is gezaehlt.
                contributed = False
                if np.isfinite(lam_pair) and lam_pair > 0:
                    d["lambda_total_pair_cm1"] += weight * lam_pair
                    contributed = True
                if np.isfinite(lam_mode) and lam_mode > 0:
                    d["lambda_total_mode_cm1"] += weight * lam_mode
                    contributed = True
                if contributed:
                    d["n_modes_contributing"] += 1
        
        # Spektren-Raster: 0 bis max(freq)+10, in 0.5-cm-1-Schritten
        # (passt to the NRVS-Aufloesung; per TOML konfigurierbar)
        sig = float(getattr(cfg, "reorg_spectrum_sigma_cm1", 5.0))
        step = float(getattr(cfg, "reorg_spectrum_step_cm1", 0.5))
        max_f = float(np.max(freqs)) + 10.0 if len(freqs) > 0 else 800.0
        grid = np.arange(0.0, max_f + step, step)
        modspec = compute_modulation_spectra(freqs, per_mode_aggs, grid,
                                              sigma_cm1=sig)
        cospec = compute_co_modulation_spectra(modspec)
        cum = compute_cumulative_reorganization(freqs, per_mode_aggs, grid,
                                                 use_mode_mass=True)
        return {
            "totals":               totals,
            "sub_totals":           sub_totals,
            "spectra_grid_cm1":     grid,
            "modulation_spectra":   modspec,
            "co_modulation_spectra": cospec,
            "cumulative_reorg":     cum,
            "spectra_sigma_cm1":    sig,
            "cumulative_uses_mode_mass": True,
        }
    
    # Globale Aggregate
    _v37_aggregates_global = _compute_v37_aggregates(results)
    if _v37_aggregates_global is not None:
        # In the Results als Diagnose entry in the ERSTEN Result speichern.
        # Der Export-Code holt sie von dort for the System-Sheets.
        results[0]["_v37_aggregates"] = _v37_aggregates_global
        # REPORT entry with the Total-Reorg-valuesn
        from .reorganization import CHANNELS
        _t = _v37_aggregates_global["totals"]
        _lambda_msg = "Marcus-Hush total reorg (system sums): " + ", ".join(
            f"Lambda_{ch}={_t[ch]['lambda_total_mode_cm1']:.2f} cm-1"
            for ch in CHANNELS if _t[ch]["lambda_total_mode_cm1"] > 0
        )
        runlog.info(_lambda_msg)
        print(f"  {_lambda_msg}")

    # ── MULTI-FENSTER-MODUS ───────────────────────────────────────────────
    if _multi_window:
        import copy as _copy
        windows = cfg.get_windows()
        _base_outdir = (cfg.output_dir if cfg.output_dir
                        else os.path.dirname(os.path.abspath(cfg.log_file)))
        print(f"\n  multi-window-Export: {len(windows)} Window + Gesamt")

        # ── Top-Level: Gesamtauswertung 0–max (alle Modes) ───────────────
        # v3.4: Lukas hat berechtigt kritisiert, dass the Top-Level-
        # directory bisher leer blieb. Jetzt exportieren wir hier
        # the volle Auswertung over alle Modes (im output_dir),
        # _additionally_ to the per-Fenster-Subdirs.
        print(f"\n  Overall analysis (all {len(results)} modes)...")
        os.makedirs(_base_outdir, exist_ok=True)

        ges_cfg = _copy.copy(cfg)
        ges_cfg.freq_windows = None     # einmaliger Single-Window-Export

        # Embeddings over alle Modes
        print(f"    Embeddings ({len(results)} Moden)...")
        X_b_g, X_e_g, feat_b_g, feat_e_g = build_feature_matrix(
            results, coord_info)
        emb_coords_g, cluster_data_g = compute_embeddings(
            X_b_g, X_e_g, feat_e_g, results, runlog)
        if cluster_data_g:
            _cl_sum_g = {}
            for _m, (_lb, _ch, _ci) in cluster_data_g.items():
                _la = np.array(_lb)
                _cl_sum_g[_m] = {
                    "n_clusters": len([k for k in _ci if k >= 0]),
                    "n_noise":    int((_la == -1).sum()),
                    "n_total":    len(_la),
                }
            runlog.cluster_summary = _cl_sum_g

        # SS-UMAP over alle
        ss_umap_g = None
        if ss_elements and any(r.get("ss") for r in results):
            try:
                ss_umap_g = compute_ss_umap_cluster(
                    results, ss_elements, runlog=runlog)
            except Exception as _eg:
                runlog.warn(f"SS-UMAP Gesamt: {_eg}")

        # Ca-Daten over alle
        ca_data_g = _build_ca_data(ca_pre, ca_amps_by_mode, results) \
                    if ca_pre is not None else None

        # Ca-UMAP over alle (new in v1.0.3)
        ca_umap_g = None
        if ca_data_g is not None:
            try:
                ca_umap_g = compute_ca_umap_cluster(
                    results, ca_data_g, runlog=runlog)
            except Exception as _eg:
                runlog.warn(f"Ca-UMAP Gesamt: {_eg}")

        ges_payload = ExportPayload(
            results              = results,
            coord_info           = coord_info,
            dist_ref             = dist_ref,
            logname              = base,
            cfg                  = ges_cfg,
            runlog               = runlog,
            cluster_info         = _cluster_info,
            b_factors            = b_factors,
            atoms                = atoms,
            context_results      = [],
            context_results_left = [],
            embedding_coords     = emb_coords_g,
            embed_feat_matrix    = X_e_g,
            embed_feat_names     = feat_e_g,
            cluster_data         = cluster_data_g,
            ca_data              = ca_data_g,
            ss_umap_data         = ss_umap_g,
            ca_umap_data         = ca_umap_g,
        )
        print(f"    Exporting overall analysis...")
        export_all(ges_payload)
        print(f"    \u2713 Gesamt \u2192 {_base_outdir}")

        # ── Pro Fenster: Sub-Auswertung in Unterverzeichnis ──────────────
        n_windows = len(windows)
        for _win_idx, (win_lo, win_hi) in enumerate(windows):
            _hi_fin = win_hi != float("inf")
            # v1.0.4 bugfix: half-open intervals [lo, hi) so a mode with
            # frequency exactly equal to a window boundary lands in
            # exactly ONE window, not two. The very last window stays
            # closed [lo, hi] so its upper edge mode is not lost. With
            # the default 100 cm-1-quantized boundaries this is mostly
            # theoretical (real frequencies are floats), but the
            # half-open convention is the standard for binning and
            # prevents silent double-counting of B-factor and
            # Lambda_cumulative contributions if a quantum chemistry
            # code emits a frequency that lands exactly on a boundary.
            _is_last = (_win_idx == n_windows - 1)
            win_label = (f"{win_lo:.0f}-{win_hi:.0f} cm-1" if _hi_fin
                         else f"{win_lo:.0f}-max cm-1")
            print(f"\n  Window {win_label}...")

            # modes for dieses Window filtern (half-open except last)
            if _is_last or not _hi_fin:
                win_results = [r for r in results
                                if r["freq"] >= win_lo and
                                (not _hi_fin or r["freq"] <= win_hi)]
            else:
                win_results = [r for r in results
                                if win_lo <= r["freq"] < win_hi]
            if not win_results:
                print(f"    No modes -- skipped.")
                runlog.warn(f"Window {win_label}: no modes, skipped.")
                continue

            # Fenster-Config: freq_min/max setzen → outdir() erzeugt Unterordner
            win_cfg = _copy.copy(cfg)
            win_cfg.freq_min     = win_lo
            win_cfg.freq_max     = win_hi if _hi_fin else None
            win_cfg.freq_windows = None   # no rekursiver multi-window
            os.makedirs(win_cfg.outdir(), exist_ok=True)

            # Kontext-modes from Gesamtergebnissen (kein additionallyer Analyselauf)
            ctx_hi   = cfg.interp_context_cm1
            ctx_right = ([r for r in results
                           if _hi_fin and win_hi < r["freq"] <= win_hi + ctx_hi])
            ctx_left  = [r for r in results
                         if win_lo - ctx_hi <= r["freq"] < win_lo]

            # B-Faktoren for dieses Window from gespeicherten Beicontribute
            b_accum_w = np.zeros(len(atoms))
            for r in win_results:
                for _ba, _contrib in r.get("_b_contribs", {}).items():
                    b_accum_w[_ba] += _contrib
            b_factors_w = 8. * np.pi**2 * b_accum_w / 3.

            # v3.7-Aggregate for dieses Fenster
            _v37_w = _compute_v37_aggregates(win_results)
            if _v37_w is not None and win_results:
                win_results[0]["_v37_aggregates"] = _v37_w

            # Ca-Daten for dieses Fenster
            ca_data_w = _build_ca_data(ca_pre, ca_amps_by_mode, win_results)

            # Embeddings for dieses Fenster
            print(f"    Embeddings ({len(win_results)} Moden)...")
            X_b_w, X_e_w, feat_b_w, feat_e_w = build_feature_matrix(
                win_results, coord_info)
            emb_coords_w, cluster_data_w = compute_embeddings(
                X_b_w, X_e_w, feat_e_w, win_results, runlog)

            if cluster_data_w:
                _cl_sum_w = {}
                for _m, (_lb, _ch, _ci) in cluster_data_w.items():
                    _la = np.array(_lb)
                    _cl_sum_w[_m] = {
                        "n_clusters": len([k for k in _ci if k >= 0]),
                        "n_noise":    int((_la == -1).sum()),
                        "n_total":    len(_la),
                    }
                runlog.cluster_summary = _cl_sum_w

            # SS-UMAP for dieses Fenster
            ss_umap_w = None
            if ss_elements and any(r.get("ss") for r in win_results):
                try:
                    ss_umap_w = compute_ss_umap_cluster(
                        win_results, ss_elements, runlog=runlog)
                except Exception as _ew:
                    runlog.warn(f"SS-UMAP Window {win_label}: {_ew}")

            # Ca-UMAP for dieses Fenster (new in v1.0.3)
            ca_umap_w = None
            if ca_data_w is not None:
                try:
                    ca_umap_w = compute_ca_umap_cluster(
                        win_results, ca_data_w, runlog=runlog)
                except Exception as _ew:
                    runlog.warn(f"Ca-UMAP Window {win_label}: {_ew}")

            # Payload: per Fenster
            win_payload = ExportPayload(
                results              = win_results,
                coord_info           = coord_info,
                dist_ref             = dist_ref,
                logname              = base,
                cfg                  = win_cfg,
                runlog               = runlog,
                cluster_info         = _cluster_info,
                b_factors            = b_factors_w,
                atoms                = atoms,
                context_results      = ctx_right,
                context_results_left = ctx_left,
                embedding_coords     = emb_coords_w,
                embed_feat_matrix    = X_e_w,
                embed_feat_names     = feat_e_w,
                cluster_data         = cluster_data_w,
                ca_data              = ca_data_w,
                ss_umap_data         = ss_umap_w,
                ca_umap_data         = ca_umap_w,
            )
            print(f"    Export...")
            export_all(win_payload)
            print(f"    \u2713 {len(win_results)} modes \u2192 {win_cfg.outdir()}")

        # Hilfsfelder aufraeumen
        for r in results: r.pop("_b_contribs", None)

        # BEFUND einmalig in the Wurzel-output_dir
        befund = os.path.join(_base_outdir, base + "_REPORT.txt")
        runlog.write_befund(befund)

    else:
        # ── EINZELFENSTER-MODUS (unveraendertes altes Verhalten) ─────────────
        runlog.info("Phase 4: embeddings/clustering")
        print("  Embeddings...")
        X_b, X_e, feat_b, feat_e = build_feature_matrix(results, coord_info)
        print(f"    Feature matrix: base={X_b.shape[1]}D, extended={X_e.shape[1]}D")
        embedding_coords, cluster_data = compute_embeddings(
            X_b, X_e, feat_e, results, runlog)

        if cluster_data:
            _cl_sum = {}
            for method, (labels, chars, cids) in cluster_data.items():
                import numpy as _np
                lbl_arr = _np.array(labels)
                _cl_sum[method] = {
                    "n_clusters": len([k for k in cids if k >= 0]),
                    "n_noise":    int((lbl_arr == -1).sum()),
                    "n_total":    len(lbl_arr),
                }
            runlog.cluster_summary = _cl_sum

        ss_umap_data = None
        if ss_elements and any(r.get("ss") for r in results):
            print("  SS-UMAP-Clustering...")
            try:
                ss_umap_data = compute_ss_umap_cluster(
                    results, ss_elements, runlog=runlog)
                if ss_umap_data and ss_umap_data[0] is not None:
                    _ss_labels = ss_umap_data[1]
                    import numpy as _np2
                    _ss_arr = _np2.array([l for l in _ss_labels if l != -99])
                    runlog.cluster_summary["SS_UMAP"] = {
                        "n_clusters": len(set(_ss_arr) - {-1}) if len(_ss_arr) else 0,
                        "n_noise":    int((_ss_arr == -1).sum()) if len(_ss_arr) else 0,
                        "n_total":    len(_ss_arr),
                    }
            except Exception as e:
                runlog.warn(f"SS-UMAP failed: {e}")

        # Ca-UMAP-Clustering (new in v1.0.3)
        ca_umap_data = None
        if ca_data is not None:
            print("  Ca-UMAP-Clustering...")
            try:
                ca_umap_data = compute_ca_umap_cluster(
                    results, ca_data, runlog=runlog)
                if ca_umap_data and ca_umap_data[0] is not None:
                    _ca_labels = ca_umap_data[1]
                    import numpy as _np3
                    _ca_arr = _np3.array([l for l in _ca_labels if l != -99])
                    runlog.cluster_summary["Ca_UMAP"] = {
                        "n_clusters": len(set(_ca_arr) - {-1}) if len(_ca_arr) else 0,
                        "n_noise":    int((_ca_arr == -1).sum()) if len(_ca_arr) else 0,
                        "n_total":    len(_ca_arr),
                    }
            except Exception as e:
                runlog.warn(f"Ca-UMAP failed: {e}")

        print("  Export...")
        _payload = ExportPayload(
            results              = results,
            coord_info           = coord_info,
            dist_ref             = dist_ref,
            logname              = base,
            cfg                  = cfg,
            runlog               = runlog,
            cluster_info         = _cluster_info,
            b_factors            = b_factors,
            atoms                = atoms,
            context_results      = context_results,
            context_results_left = context_results_left,
            embedding_coords     = embedding_coords,
            embed_feat_matrix    = X_e,
            embed_feat_names     = feat_e,
            cluster_data         = cluster_data,
            ca_data              = ca_data,
            ss_umap_data         = ss_umap_data,
            ca_umap_data         = ca_umap_data,
        )
        export_all(_payload)
        print("  Embedding PNGs: done")

        for r in results: r.pop("_b_contribs", None)

        befund = cfg.outname("_REPORT.txt")
        runlog.write_befund(befund)

    # ── Abschluss (beide Modi) ─────────────────────────────────────────────
    elapsed = time.time() - t_start
    n_ip  = sum(1 for r in results if r["mode_type"] == "In-plane")
    n_oop = sum(1 for r in results if r["mode_type"] == "Out-of-plane")
    n_mx  = sum(1 for r in results if r["mode_type"] == "Torsional/Mixed")

    print("\n" + "="*70)
    print(f"  DONE  ({elapsed/60:.1f} min / {elapsed:.0f} s)")
    print("="*70)
    print(f"  modes analyzed:   {len(results)}")
    print(f"  HP eigvec:          {n_hp}/{len(results)}")
    if n_std: print(f"  standard fallback:  {n_std}/{len(results)}")
    runlog.info(f"Result: In-plane={n_ip}, Out-of-plane={n_oop}, "
                f"Torsional={n_mx}, Failed={n_fail}")
    print(f"  In-plane:           {n_ip}")
    print(f"  Out-of-plane:       {n_oop}")
    print(f"  Torsional:          {n_mx}")
    if n_fail:           print(f"  Failed:     {n_fail}")
    if runlog.warnings:  print(f"  Warnings:          {len(runlog.warnings)}")
    if runlog.errors:    print(f"  ERRORS:             {len(runlog.errors)}")
    if _multi_window:
        print(f"  Windows:            {len(cfg.get_windows())}")
        print(f"  Output:            {_base_outdir}")
    else:
        print(f"  Output:            {cfg.outdir()}")
    print(f"  Report file:        {os.path.basename(befund)}")

    return 0


# ===========================================================================
#  Multi-Cluster-Wrapper
# ===========================================================================

def run_analysis(cfg: "Config") -> int:
    """Programmatic entry point: analyzes the normal modes of a
    [2Fe-2S] system and writes Excel outputs.

    This function is the official API of ``modenanalyse_2fe2s``.
    It can be used interactively (Spyder, Jupyter), from scripts,
    or via the CLI.

    **Single-cluster mode** (default, ``cfg.analyze_all_clusters=False``)
    delegates directly to :func:`_run_analysis_single` and analyzes
    exactly the cluster selected by ``cfg.cluster_index``.

    **Multi-cluster mode** (``cfg.analyze_all_clusters=True``) is
    intended for dimers and multi-[2Fe-2S] systems (glutaredoxin
    dimers, domain constructs with two clusters, etc.). Before the
    actual run, the atom list is extracted and all [2Fe-2S] clusters
    are identified via :func:`find_all_clusters`. Then
    :func:`_run_analysis_single` is called once per cluster, each
    with its own subfolder ``output_dir/cluster_<N>/``. At the end,
    the wrapper writes ``output_dir/multi_cluster_summary.txt`` with
    the status per cluster.

    Parameters
    ----------
    cfg : Config
        Complete configuration. Typical creation via
        ``from modenanalyse_2fe2s import Config`` and setting the
        fields ``log_file``, ``pdb_file``, ``output_dir``, etc.

    Returns
    -------
    int
        Exit code: 0 = all clusters successful, 1 = at least one
        cluster failed (details in the multi-cluster summary).
    """
    if not cfg.analyze_all_clusters:
        return _run_analysis_single(cfg)

    # ── Multi-Cluster-Modus ────────────────────────────────────────────
    import dataclasses as _dc
    import os as _os
    from pathlib import Path as _Path
    from .geometry import find_all_clusters as _find_all_clusters

    print("\n" + "=" * 72)
    print("  modenanalyse_2fe2s -- multi-cluster analysis")
    print("=" * 72)

    # 1. Pre-Scan: atom list extrahieren (ohne vollen Lauf)
    is_orca = cfg.log_file.lower().endswith(".hess")
    if is_orca:
        from .orca_io import (load_orca_hess as _load_orca_hess,
                              parseresult_to_atoms as _pra)
        _pr = _load_orca_hess(cfg.log_file)
        atoms, _ = _pra(_pr, include_hydrogen=False)
    else:
        from .logio import (scan_log as _scan_log,
                            read_std_orient as _read_so)
        so_off, _, _ = _scan_log(cfg.log_file, cfg)
        if not so_off:
            print("[ERROR] No Standard-orientation found in log file.")
            return 1
        atoms, _ = _read_so(cfg.log_file, so_off[-1], include_hydrogen=False)

    all_clusters = _find_all_clusters(atoms, cfg)
    n_clusters = len(all_clusters)

    if n_clusters == 0:
        print("[ERROR] No [2Fe-2S] cluster found in log file.")
        return 1

    if n_clusters == 1:
        print("\n[i] Only 1 cluster in system -- analyze_all_clusters=True has\n"
              "    no effect (runs as single cluster).\n")
        return _run_analysis_single(cfg)

    print(f"\n[i] {n_clusters} [2Fe-2S]-clusters found in system:")
    for _i, (_fe, _s, _g) in enumerate(all_clusters):
        print(f"    Cluster #{_i}: Fe-Fe = {_g['fe_fe']:.3f} A, "
              f"Fe-S = {_g['fe_s_min']:.3f}--{_g['fe_s_max']:.3f} A")

    # 2. Pro Cluster ein Lauf with eigenem Subordner
    base_outdir = cfg.output_dir
    _Path(base_outdir).mkdir(parents=True, exist_ok=True)
    summaries = []
    overall_status = 0

    for cl_idx in range(n_clusters):
        sub_outdir = _os.path.join(base_outdir, f"cluster_{cl_idx}")
        print("\n" + "-" * 72)
        print(f"  Cluster #{cl_idx} -> {sub_outdir}")
        print("-" * 72)

        # Config-Kopie with angepasstem cluster_index and output_dir.
        # IMPORTANT: analyze_all_clusters on False setzen, sonst infinite 
        # Rekursion.
        cfg_sub = _dc.replace(
            cfg,
            cluster_index=cl_idx,
            output_dir=sub_outdir,
            analyze_all_clusters=False,
        )

        try:
            rc = _run_analysis_single(cfg_sub)
            status = "OK" if rc == 0 else f"FAILED (rc={rc})"
            if rc != 0:
                overall_status = 1
        except Exception as e:
            status = f"EXCEPTION: {type(e).__name__}: {e}"
            overall_status = 1
        summaries.append((cl_idx, sub_outdir, status,
                          all_clusters[cl_idx][2]))

    # 3. Multi-Cluster-Summary schreiben
    summary_file = _os.path.join(base_outdir, "multi_cluster_summary.txt")
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write(f"modenanalyse_2fe2s -- multi-cluster analysis\n")
        f.write(f"=" * 60 + "\n\n")
        f.write(f"Log file: {cfg.log_file}\n")
        f.write(f"PDB:     {cfg.pdb_file or '(none)'}\n")
        f.write(f"Number of clusters: {n_clusters}\n\n")
        for cl_idx, sub_outdir, status, geom in summaries:
            f.write(f"Cluster #{cl_idx}\n")
            f.write(f"  Geometry: Fe-Fe={geom['fe_fe']:.3f} A, "
                    f"Fe-S={geom['fe_s_min']:.3f}--{geom['fe_s_max']:.3f} A\n")
            f.write(f"  Status:    {status}\n")
            f.write(f"  Output:    {sub_outdir}\n\n")

    print("\n" + "=" * 72)
    print(f"  Multi-cluster analysis completed")
    print(f"  Summary: {summary_file}")
    print("=" * 72 + "\n")

    return overall_status
