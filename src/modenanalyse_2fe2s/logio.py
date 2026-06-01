# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

# -*- coding: utf-8 -*-
"""
logio.py
==================
file-I/O for Gaussian ``.log``-fileen and PDB-fileen.

Oeffentliche functions
-----------------------
scan_log
    Scannt a Gaussian-Logdatei and gibt Byte-Offsets .
read_std_orient
    Reads the atom list from a Standard-orientation block.
read_all_meta
    Reads Metadaten aller Frequenz-Gruppen.
get_eigvec
    Reads the eigenvector of a single mode.
parse_pdb
    Reads a PDB-file complete.

classes
-------
RunLog
    Sammelt Meldungen and Statistiken a Analyselaufs.
BlockInfo
    Metadaten a Frequenz-Blocks.

Bugfixes (gegenvia Vorversion)
---------------------------------
B3  scan_log:        Tail-Overlap-Duplikate dedupliziert.
B4  data_offset:     ``f.tell()`` is VOR ``readline()`` gespeichert.
B5  data_offset:     ``-1`` als error sentinel validates.
B8  parse_pdb:       HETATM-Records are ungefiltert geladen
                     (FES kann in anderer Kette liegen als Protein).
B9  read_std_eigvec: Flexible Regex; no ``break`` for kurzen Zeilen.
B10 read_std_orient: Optionale Type-Spalte in Regex.
B11 read_std_eigvec: Doppelte Center are dedupliziert.
B12 read_hp_eigvec:  ``ci_idx < 1`` or ``> 3`` → ``break``.
B13 _el_from_name:   Zwei-Zeichen-elements (FE, ZN, …) korrekt erkannt.
B14 parse_pdb:       Nur erstes MODEL-Record is geladen.
"""
from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .config import Config



# ---------------------------------------------------------------------------
# element-Tabelle
# ---------------------------------------------------------------------------

_ELEM: Dict[int, str] = {
    1:"H",  6:"C",  7:"N",  8:"O",  9:"F",  11:"Na", 12:"Mg",
    14:"Si",15:"P", 16:"S", 17:"Cl",19:"K",  20:"Ca",
    25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu", 30:"Zn",
    34:"Se",35:"Br",53:"I", 79:"Au",80:"Hg",
}

_TWO_CHAR_ELEM: frozenset = frozenset([
    "FE","ZN","MG","CA","MN","CU","NI","CO","MO","SE","CL","BR",
    "NA","AL","SI","CR","GA","AS","RB","SR","AG","CD","SN","CS",
    "BA","IR","PT","AU","HG","PB",
])

_HIS: frozenset = frozenset(["HIS","HIE","HID","HIP","HSD","HSE","HSP"])
_CYS: frozenset = frozenset(["CYS","CYX","CYM","CYD"])


# ===========================================================================
# RunLog
# ===========================================================================

class RunLog:
    """Sammelt alle Meldungen and Statistiken a Analyselaufs.

    Parameters
    ----------
    cfg : Config
        Konfigurationsobjekt of the laufenden Analyselaufs.

    Attributes
    ----------
    warnings : list of str
        Alle Warnmeldungen (``warn()``-Aufrufe).
    errors : list of str
        Alle Errormeldungen (``error()``-Aufrufe).
    parse_failures : list of dict
        Fehlgeschlagene modes with ``mode_num``, ``freq``, ``reason``.
    module_status : dict of {str: bool}
        Verfuegbarkeit optionaler Python-Module.
    output_files : list of (str, float)
        Erzeugte output fileen als ``(path, Groesse_MB)``.
    geometry : dict
        Geometrische Kennzahlen of the cluster.
    group_match : dict
        Ergebnis of the Kabsch-Alignments and the Koordinationserkennung.
    mode_stats : dict
        Statistiken the modenanalyse.
    """

    def __init__(self, cfg: Config) -> None:
        """Initializes the RunLog with empty log.

        Parameters
        ----------
        cfg : Config
            Analysekonfiguration (wird for the BEFUND referenziert).
        """
        self.cfg              = cfg
        self._start           = time.time()
        self._entries: List[Tuple[str, str]] = []
        self.warnings:  List[str]  = []
        self.errors:    List[str]  = []
        self.parse_failures: List[Dict] = []
        self.module_status:  Dict[str, bool] = {}
        self.output_files:   List[Tuple[str, float]] = []
        self.geometry:       Dict = {}
        self.group_match:    Dict = {}
        self.mode_stats:     Dict = {}
        # Erweiterte REPORT-Felder
        self.coord_summary:    Dict = {}  # ligands, H-N, Cluster-Geometrie
        self.results_summary:  Dict = {}  # Modentypen, OOP-Verteilung
        self.cluster_summary:  Dict = {}  # Cluster per Methode
        self.match_stats:      Dict = {}  # PDB matching: n_matched, n_total, mean_d, max_d, n_ambiguous

    # ------------------------------------------------------------------
    def info(self, msg: str) -> None:
        """Adds an INFO entry (appears in REPORT under NOTES).

        Parameters
        ----------
        msg : str
            Meldungstext.
        """
        self._entries.append(("INFO", msg))

    def warn(self, msg: str) -> None:
        """Adds a WARN entry (appears in REPORT under WARNINGS).

        Parameters
        ----------
        msg : str
            Meldungstext.
        """
        self._entries.append(("WARN", msg))
        self.warnings.append(msg)

    def error(self, msg: str) -> None:
        """Adds an ERROR entry (appears in REPORT under ERRORS).

        Parameters
        ----------
        msg : str
            Meldungstext.
        """
        self._entries.append(("ERROR", msg))
        self.errors.append(msg)

    def add_parse_failure(self, mode_num: int, freq: float,
                          reason: str) -> None:
        """Dokumentiert a fehlgeschlagene Mode.

        Parameters
        ----------
        mode_num : int
            Modenummer from the Gaussian-Ausgabe.
        freq : float
            Frequency the mode in cm⁻^1.
        reason : str
            Kurze Errorbeschreibung.
        """
        self.parse_failures.append(
            {"mode_num": mode_num, "freq": freq, "reason": reason})

    def add_output(self, path: str) -> None:
        """Registriert a erzeugte output file.

        Parameters
        ----------
        path : str
            Vollstaendiger Path to the output file.
        """
        try:
            mb = os.path.getsize(path) / 1e6
        except OSError:
            mb = 0.0
        self.output_files.append((path, mb))

    def write_befund(self, path: str) -> None:
        """Writes the run report as a text file.

        Parameters
        ----------
        path : str
            Target file for the report (will be overwritten).
        """
        elapsed = time.time() - self._start
        lines: List[str] = []
        sep = "─" * 62

        def h(title: str) -> None:
            """Adds a formatted section header to the REPORT lines."""
            lines.append(f"\n{sep}\n  {title}\n{sep}\n")

        lines.append("RUN REPORT - runner\n")
        lines.append(f"Created:  {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        lines.append(f"Runtime:  {elapsed/60:.1f} min ({elapsed:.0f} s)\n")

        cfg = self.cfg
        h("CONFIGURATION")
        lines.append(f"  log_file:    {cfg.log_file}\n")
        lines.append(f"  pdb_file:    {cfg.pdb_file or '(none)'}\n")
        lines.append(f"  output_dir:  {cfg.outdir()}\n")
        lines.append(
            f"  freq_filter: "  # freq_min/max: None-safe
            f"{'-' if cfg.freq_min is None else cfg.freq_min} - "
            f"{'-' if cfg.freq_max is None else cfg.freq_max} cm-1\n")
        # Bugfix v1.0.4 (post-release Apd1 audit): also report freq_windows
        # if set; previously only freq_min/freq_max were shown, which made
        # the actual filter invisible when only freq_windows was used.
        _fw = getattr(cfg, "freq_windows", None)
        if _fw:
            try:
                _fw_str = ", ".join(f"[{lo:.0f}-{hi:.0f}]" for lo, hi in _fw)
            except Exception:
                _fw_str = str(_fw)
            lines.append(f"  freq_windows: {_fw_str} cm-1\n")
        lines.append(f"  temp_k:      {cfg.temp_k} K\n")
        lines.append(f"  interp_step: {cfg.interp_step} cm-1\n")
        lines.append(f"  interp_boundary_mode: {getattr(cfg, 'interp_boundary_mode', 'context')}\n")
        lines.append(f"  pdb_chain:   '{cfg.pdb_chain}'\n")
        lines.append(
            f"  H atoms:     "
            f"{'yes' if cfg.include_hydrogen else 'no'}\n")

        h("INSTALLED MODULES")
        # Some modules have different names in pip than in import
        _PIP_NAMES = {
            "umap":     "umap-learn",
            "sklearn":  "scikit-learn",
            "scsdpy":   "scsdpy",
            "hdbscan":  "hdbscan",
        }
        for mod, ok in self.module_status.items():
            lines.append(f"  {'OK     ' if ok else 'MISSING'} - {mod}\n")
            if not ok:
                pip_name = _PIP_NAMES.get(mod.split()[0], mod.split()[0])
                lines.append(f"          → pip install {pip_name}\n")

        h("GEOMETRY")
        g = self.geometry
        if g:
            lines.append(f"  Atoms (heavy): {g.get('n_heavy','?')}\n")
            lines.append(f"  Atoms (total): {g.get('n_total','?')}\n")
            if "fe_fe" in g:
                lines.append(f"  Fe-Fe:          {g['fe_fe']:.6f} A\n")
            # Fe-S individual values structured
            for k, v in g.get("fe_s_distances", {}).items():
                lines.append(f"  {k:12s}:     {v:.6f} A\n")
            if "fe_s_mean" in g:
                lines.append(f"  Fe-S (mean):    {g['fe_s_mean']:.6f} A\n")
            # Fe-S-Fe angles
            for ai, ang in enumerate(g.get("fe_s_fe_angles_deg", [])):
                lines.append(f"  Fe-S{ai+1}-Fe:      {ang:.2f} deg\n")
            if "fe_s_fe_mean_deg" in g:
                lines.append(f"  Fe-S-Fe (mean):  {g['fe_s_fe_mean_deg']:.2f} deg\n")
            # Cluster normal + folding residual
            cn = g.get("cluster_normal")
            if cn is not None and len(cn) == 3:
                lines.append(
                    f"  Cluster normal n_hat: "
                    f"({cn[0]:+.4f}, {cn[1]:+.4f}, {cn[2]:+.4f})\n")
            res = g.get("cluster_plane_residual_a")
            if res is not None:
                lines.append(
                    f"  Folding residual:    {res:.4f} A "
                    f"(RMS distance of 4 cluster atoms from best-fit plane)\n")

        h("COORDINATING AMINO ACIDS (auto-detected)")
        gm  = self.group_match
        cs  = self.coord_summary
        if gm:
            rmsd = gm.get("kabsch_rmsd")
            if rmsd is not None:
                lines.append(f"  Kabsch RMSD:    {rmsd:.4f} A\n")
        # Ligands from coord_summary (richer)
        for lig in cs.get("ligands", []):
            prot_info = ""
            if lig.get("his_protonated"):
                hn_n   = lig.get("hn_n_type", "?")
                hn_len = lig.get("hn_len", 0.)
                via    = lig.get("hn_via", "?")
                prot_info = (f"  → H-N detected: {hn_n}-H = {hn_len:.3f} A"
                             f" (via {via})\n")
            elif lig.get("element") == "N":
                prot_info = "  → His not protonated (no H on ring N)\n"
            lines.append(
                f"  Fe{lig.get('fe_idx',0)+1} ← {lig.get('res_label','?')}"  
                f" ({lig.get('element','?')},"  
                f" {lig.get('aname','?')},"  
                f" d={lig.get('bond_len',0.):.3f} A)\n")
            if prot_info:
                lines.append(prot_info)
        if not cs.get("ligands") and not gm:
            lines.append("  (no PDB / no Kabsch)\n")
        # Cluster distances
        cd = cs.get("cluster_distances", {})
        if cd:
            lines.append("\n  Cluster equilibrium geometry:\n")
            for k, v in cd.items():
                lines.append(f"    {k:<12s}: {v:.6f} A\n")

        # Matching quality
        ms = self.match_stats
        if ms:
            lines.append(
                f"\n  PDB matching: {ms.get('n_matched',0)}/{ms.get('n_total',0)} atoms "
                f"({ms.get('n_matched',0)/max(ms.get('n_total',1),1)*100:.1f}%), "
                f"mean distance {ms.get('mean_d',0.):.4f} A, "
                f"max {ms.get('max_d',0.):.4f} A, "
                f"ambiguous {ms.get('n_ambiguous',0)}\n")

        # Cluster geometry info and SCSD reference
        _cl_geo = [msg for lvl, msg in self._entries
                   if lvl == "INFO" and ("Cluster-Geometrie" in msg or "Cluster geometry" in msg)]
        if _cl_geo:
            lines.append(f"\n  {_cl_geo[0]}\n")
        _scsd_ref = [msg for lvl, msg in self._entries
                     if lvl == "INFO" and ("SCSD-Referenz" in msg or "SCSD reference" in msg)]
        if _scsd_ref:
            lines.append(f"  {_scsd_ref[0]}\n")

        h("MODE ANALYSIS")
        ms = self.mode_stats
        if ms:
            lines.append(f"  Found (total):     {ms.get('n_found','?')}\n")
            lines.append(f"  After filter:      {ms.get('n_filtered','?')}\n")
            lines.append(f"  HP eigvec:         {ms.get('n_hp','?')}/{ms.get('n_filtered','?')}\n")
            if ms.get('n_std_fallback', 0):
                lines.append(f"  default fallback: {ms.get('n_std_fallback','?')}/{ms.get('n_filtered','?')}\n")
            lines.append(
                f"  Failed:            {len(self.parse_failures)}\n")

        # Result distribution
        rs = self.results_summary
        if rs:
            h("MODE DISTRIBUTION")
            total = rs.get("n_modes", 0)
            for mtype, count in rs.get("mode_types", {}).items():
                pct = 100*count/total if total else 0
                lines.append(f"  {mtype:<22s}: {count:5d}  ({pct:.1f}%)\n")
            if "freq_range" in rs:
                lo, hi = rs["freq_range"]
                lines.append(f"  Frequency range analyzed: {lo:.2f} - {hi:.2f} cm-1\n")
            if "mean_oop" in rs:
                lines.append(f"  Mean OOP fraction: {rs['mean_oop']:.1f}%\n")
            if "n_scsd" in rs:
                lines.append(f"  SCSD computed for: {rs['n_scsd']} modes\n")

        # Cluster analysis
        cl = self.cluster_summary
        if cl:
            h("CLUSTER ANALYSIS (embeddings)")
            for method, info in cl.items():
                nc  = info.get("n_clusters", 0)
                ns_ = info.get("n_noise", 0)
                nt  = info.get("n_total", 0)
                lines.append(
                    f"  {method:<20s}: {nc:3d} clusters, "
                    f"{ns_:4d} noise"
                    f"{f' ({ns_*100//nt if nt else 0}%)'  if nt else ''}\n")

        if self.parse_failures:
            h("FAILED MODES")
            for pf in self.parse_failures[:50]:
                lines.append(
                    f"  Mode {pf['mode_num']:6d} @ "
                    f"{pf['freq']:8.3f} cm-1: {pf['reason']}\n")
            if len(self.parse_failures) > 50:
                lines.append(
                    f"  ... and {len(self.parse_failures)-50} more\n")

        if self.warnings:
            h("WARNINGS")
            _max_warn = 100
            for w in self.warnings[:_max_warn]:
                lines.append(f"  ! {w}\n")
            if len(self.warnings) > _max_warn:
                lines.append(
                    f"  ... and {len(self.warnings)-_max_warn} more warnings\n")

        if self.errors:
            h("ERRORS")
            for e in self.errors:
                lines.append(f"  !! {e}\n")

        # INFO entries (assignment decisions, fallbacks, notes)
        info_msgs = [msg for lvl, msg in self._entries if lvl == "INFO"]
        if info_msgs:
            h("NOTES")
            for msg in info_msgs:
                lines.append(f"  [i] {msg}\n")

        h("OUTPUT FILES")
        for pf, mb in self.output_files:
            lines.append(
                f"  {os.path.basename(pf):<55s} {mb:6.1f} MB\n")

        with open(path, "w", encoding="utf-8") as f:
            f.writelines(lines)
        self.info(f"REPORT written: {os.path.basename(path)}")
        print(f"  → {os.path.basename(path)}")


# ===========================================================================
# Phase 1 - Scan
# ===========================================================================

def scan_log(
        filepath: str,
        cfg: Config,
        runlog: Optional["RunLog"] = None,
) -> Tuple[List[int], List[int], List[int]]:
    """Scannt a Gaussian-Logdatei and gibt Byte-Offsets .

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    cfg : Config
        Konfiguration; required ``scan_chunk_mb``.

    Returns
    -------
    so_off : list of int
        Byte-Offsets aller ``Standard orientation:``-blocks.
    nc_off : list of int
        Byte-Offsets aller ``and normal coordinates:``-Zeilen.
    fr_off : list of int
        Byte-Offsets aller ``Frequencies ---``-Zeilen.

    Notes
    -----
    Bugfix B3: Tail-Overlap-Duplikate are through ``seen``-Sets
    dedupliziert.
    """
    so_off: List[int] = []
    nc_off: List[int] = []
    fr_off: List[int] = []
    seen_so: Set[int] = set()
    seen_nc: Set[int] = set()
    seen_fr: Set[int] = set()

    chunk      = cfg.scan_chunk_mb * 1024 * 1024
    tail       = b""
    offset     = 0
    total_size = os.path.getsize(filepath)
    t0         = time.time()
    _size_msg = f"{total_size/1e9:.2f} GB"
    if runlog is not None:
        runlog.info(f"scan_log: filegröße {_size_msg}")
    print(f"    {_size_msg}")

    pat_so = re.compile(rb"Standard orientation:")
    pat_nc = re.compile(rb"and normal coordinates:")
    pat_fr = re.compile(rb"^\s*Frequencies\s*-{2,}", re.MULTILINE)

    with open(filepath, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            data = tail + block
            base = offset - len(tail)

            for m in pat_so.finditer(data):
                off = base + m.start()
                if off not in seen_so:
                    seen_so.add(off); so_off.append(off)
            for m in pat_nc.finditer(data):
                off = base + m.start()
                if off not in seen_nc:
                    seen_nc.add(off); nc_off.append(off)
            for m in pat_fr.finditer(data):
                off = base + m.start()
                if off not in seen_fr:
                    seen_fr.add(off); fr_off.append(off)

            tail    = data[-200:]
            offset += len(block)
            elapsed = max(time.time() - t0, 1e-3)
            speed   = offset / elapsed / 1e6
            eta     = max(0.0, (total_size-offset)/(offset/elapsed)) if offset else 0.0
            print(
                f"    Scan: {offset/total_size*100:.1f}%  "
                f"{speed:.0f} MB/s  ETA {eta:.0f}s    ",
                end="\r")

    print()
    _scan_msg = (f"scan_log: SO={len(so_off)}, NC={len(nc_off)}, "
                 f"FR={len(fr_off)} in {time.time()-t0:.1f}s")
    if runlog is not None:
        runlog.info(_scan_msg)
    return so_off, nc_off, fr_off


# ===========================================================================
# Phase 2 - default orientation
# ===========================================================================


# ===========================================================================
# Cache-Hilfsfunktionen
# ===========================================================================

_CACHE_VERSION = "modenanalyse_cache_v1"   # For Formatänderungen erhöhen


def _cache_key(filepath: str) -> str:
    """Creates a uniqueen Cache-key from path, Groesse, mtime and Python-Version.

    Die Python-Version is Teil of the Keys, da Pickle-Formate zwischen
    Hauptversionen (e.g. 3.11 -> 3.12) inkompatibel sein koennen.
    Nach a Python-Upgrade is the Cache automatically neu gebaut.
    """
    import hashlib, sys
    try:
        st  = os.stat(filepath)
        pyv = f"{sys.version_info.major}.{sys.version_info.minor}"
        raw = f"{os.path.abspath(filepath)}|{st.st_size}|{st.st_mtime}|py{pyv}"
        return hashlib.md5(raw.encode()).hexdigest()[:16]
    except OSError:
        return ""


def _cache_path(filepath: str) -> str:
    """Returns the path the Cache-file zurück."""
    base = os.path.splitext(os.path.basename(filepath))[0]
    key  = _cache_key(filepath)
    if not key:
        return ""
    cache_dir = os.path.join(os.path.dirname(os.path.abspath(filepath)),
                              ".modenanalyse_cache")
    return os.path.join(cache_dir, f"{base}_{key}.pkl")


def load_scan_cache(filepath: str,
                    runlog: Optional["RunLog"] = None) -> Optional[Tuple]:
    """Laedt gecachte Scan-Ergebnisse if present and gueltig.

    Returns
    -------
    tuple or None
        ``(all_blocks, best_block, cand_map, so_off, nc_off, fr_off)``
        or ``None`` if no gültiger Cache existiert.
    

    Sicherheitshinweis: Pickle-fileen (`.pkl`) sollten only aus
    vertrauenswuerdigen Quellen geladen werden. Fuer den
    typicalen Offline-Einsatz with eigenen Log-fileen ist
    the Risiko vernachlaessigbar.
    """
    import pickle
    cp = _cache_path(filepath)
    if not cp or not os.path.isfile(cp):
        return None
    try:
        with open(cp, "rb") as f:
            data = pickle.load(f)
        # Versions- and Integritätsprüfung
        if not isinstance(data, dict):
            return None
        if data.get("_version") != _CACHE_VERSION:
            return None
        if data.get("_key") != _cache_key(filepath):
            return None
        return (data["all_blocks"], data["best_block"], data["cand_map"],
                data["so_off"],    data["nc_off"],    data["fr_off"])
    except Exception as _ce:
        # For jedem Error: Cache loeschen and neu berechnen
        _warn_msg = f"Cache beschaedigt or inkompatibel – is neu erstellt ({_ce})"
        if runlog is not None:
            runlog.warn(_warn_msg)
        else:
            import warnings as _cw
            _cw.warn(_warn_msg, UserWarning)
        try:
            os.unlink(cp)
        except Exception as _ue:
            _del_msg = f"Cache file could not be deleted: {_ue}"
            if runlog is not None:
                runlog.warn(_del_msg)
            else:
                import warnings as _cw2
                _cw2.warn(_del_msg, UserWarning)
        return None


def save_scan_cache(filepath: str,
                    all_blocks: List, best_block: Dict, cand_map: Dict,
                    so_off: List[int], nc_off: List[int],
                    fr_off: List[int],
                    runlog: Optional["RunLog"] = None) -> None:
    """Writes Scan-Ergebnisse in Cache-file.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian-Logdatei (bestimmt Cache-path).
    all_blocks, best_block, cand_map : list/dict
        Scan-Ergebnisse from ``read_all_meta``.
    so_off, nc_off, fr_off : list of int
        Byte-Offsets from ``scan_log``.
    runlog : RunLog, optional
        If angiven, are Cache-Error als ``warn`` entry ins
        BEFUND/RunLog geschrieben (nicht only on stdout).

    Notes
    -----
    Schlaegt stillschweigend fehl if directory not beschreibbar.
    Error are over ``runlog.warn`` (beforezugt) and stdout gemeldet.
    """
    import pickle
    cp = _cache_path(filepath)
    if not cp:
        return
    try:
        cache_dir = os.path.dirname(cp)
        os.makedirs(cache_dir, exist_ok=True)
        data = {
            "_version":  _CACHE_VERSION,
            "_key":      _cache_key(filepath),
            "all_blocks": all_blocks,
            "best_block": best_block,
            "cand_map":   cand_map,
            "so_off":     so_off,
            "nc_off":     nc_off,
            "fr_off":     fr_off,
        }
        # Atomar schreiben: erst tmp, dann rename
        tmp = cp + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, cp)
    except Exception as _e:
        _msg = (f"Cache not schreibbar "
                f"({os.path.basename(filepath)}): {_e}")
        # Einheitlicher Meldeweg: RunLog (BEFUND) hat Vorrang
        if runlog is not None:
            runlog.warn(_msg)
        else:
            import warnings as _sw
            _sw.warn(_msg, UserWarning, stacklevel=2)


def read_std_orient(
        filepath: str,
        so_offset: int,
        include_hydrogen: bool = False,
) -> Tuple[List[Dict], Dict[int, int]]:
    """Reads the atom list from a Standard-orientation block.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    so_offset : int
        Byte-Offset the zugehoerigen ``Standard orientation:`` line.
    include_hydrogen : bool, optional
        If ``True``, are hydrogen-Atome read in.
        default is ``False``.

    Returns
    -------
    atoms : list of dict
        atom list; jedes Dict enthaelt ``center``, ``atomic_num``,
        ``symbol``, ``x``, ``y``, ``z``.
    idx_map : dict of {int: int}
        Mapping Gaussian-Center-Nummer → Index in ``atoms``.

    Raises
    ------
    ValueError
        If no Atome in the Block found wurden.

    Notes
    -----
    Bugfix B10: Regex akzeptiert jetzt also fehlendes ``Atomic Type``-Feld
    (``(?:\\s+\\d+)?`` macht the Type-Spalte optional).
    """
    atoms: List[Dict] = []
    idx_map: Dict[int, int] = {}
    pat = re.compile(
        rb"\s*(\d+)\s+(\d+)(?:\s+\d+)?\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")

    with open(filepath, "rb") as fh:
        fh.seek(so_offset)
        for _ in range(5):
            fh.readline()
        for line in fh:
            m = pat.match(line)
            if not m:
                break
            ctr  = int(m.group(1))
            anum = int(m.group(2))
            x, y, z = float(m.group(3)), float(m.group(4)), float(m.group(5))
            if anum == 1 and not include_hydrogen:
                continue
            entry = {
                "center":     ctr,
                "atomic_num": anum,
                "symbol":     _ELEM.get(anum, f"Z{anum}"),
                "x": x, "y": y, "z": z,
            }
            idx_map[ctr] = len(atoms)
            atoms.append(entry)

    if not atoms:
        raise ValueError("No atoms found in default orientation.")
    return atoms, idx_map


# ===========================================================================
# Phase 3 - Block metadata
# ===========================================================================

@dataclass
class BlockInfo:
    """Metadaten a Frequenz-Blocks in the Gaussian-Logdatei.

    Parameters
    ----------
    offset : int
        Byte-Offset the ``Frequencies ---``-Zeile.

    Attributes
    ----------
    mode_nums : list of int
        Globale Modenummern (je Spalte in the Block).
    freqs : list of float
        Frequenzen in cm⁻^1.
    red_masses : list of float
        Reduzierte Massen in AMU.
    frc_consts : list of float
        Kraftkonstanten in mDyn/A.
    syms : list of str
        Symmetrierassen.
    is_hp : bool
        ``True`` for HP-Format (hohe Praezision).
    data_offset : int
        Byte-Offset the ersten Eigenvektor-data row.
        ``-1`` bedeutet: not bestimmt / invalid.
    """

    offset:      int
    mode_nums:   List[int]   = field(default_factory=list)
    freqs:       List[float] = field(default_factory=list)
    red_masses:  List[float] = field(default_factory=list)
    frc_consts:  List[float] = field(default_factory=list)
    syms:        List[str]   = field(default_factory=list)
    is_hp:       bool        = False
    data_offset: int         = -1


def _read_block_at_freq_line(filepath: str,
                              fr_offset: int) -> Optional[BlockInfo]:
    """Reads Metadaten a einzelnen Frequenz-Blocks.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    fr_offset : int
        Byte-Offset the ``Frequencies ---``-Zeile.

    Returns
    -------
    BlockInfo or None
        Metadaten of the Blocks, or ``None`` if the Block leer ist.

    Notes
    -----
    Bugfix B4: ``data_offset`` is over ``f.tell()`` VOR dem
    ``readline()``-Aufruf bestimmt (statt ``f.tell() - len(line.encode())``).

    Bugfix B5: ``data_offset == -1`` verbleibt als error sentinel und
    is in ``get_eigvec`` explizit geprueft.
    """
    bi = BlockInfo(offset=fr_offset)

    def _parse_dashes(line: str) -> List[float]:
        """Reads Zahlenwerte after the ersten Doppelbindestrich-Marker in a Gaussian-Zeile."""
        after = re.sub(r"^.*?-{2,}\s*", "", line)
        try:
            return [float(x) for x in after.split()]
        except ValueError:
            return []

    with open(filepath, "rb") as fh:
        fh.seek(fr_offset)
        raw = fh.readline()
        if not raw:
            return None
        freq_line = raw.decode("utf-8", errors="replace")
        bi.is_hp  = bool(re.match(r"\s*Frequencies\s*---", freq_line))
        bi.freqs  = _parse_dashes(freq_line)
        n = len(bi.freqs)
        if n == 0:
            return None

        # Modenummern from the Block-Header bestimmen
        try:
            back = min(fr_offset, 300)
            fh.seek(fr_offset - back)
            prev_lines = []
            while True:
                line = fh.readline()
                if not line or fh.tell() > fr_offset:
                    break
                prev_lines.append(line.decode("utf-8", errors="replace"))
            mode_nums: List[int] = []
            for pline in reversed(prev_lines):
                parts = pline.strip().split()
                if parts and all(p.isdigit() for p in parts):
                    try:
                        mode_nums = [int(p) for p in parts]
                        break
                    except ValueError:
                        pass
            bi.mode_nums = (mode_nums if mode_nums and len(mode_nums) == n
                            else list(range(1, n + 1)))
        except Exception:
            bi.mode_nums = list(range(1, n + 1))

        # Header-Zeilen bis to ersten Datenseile lesen
        # Window 60 Zeilen (statt 25) for Berechnungen with vielen
        # optionalen Output-Feldern (Raman, Depolar, NMR, etc.)
        fh.seek(fr_offset)
        fh.readline()   # Frequencies-Zeile skip
        for _ in range(60):
            prev_pos = fh.tell()   # B4-Fix: Position VOR readline
            raw_line = fh.readline()
            if not raw_line:
                break
            line = raw_line.decode("utf-8", errors="replace")
            ll   = line.lower()

            if "red" in ll and "mass" in ll:
                bi.red_masses = _parse_dashes(line)
            elif "frc" in ll or ("force" in ll and "const" in ll):
                bi.frc_consts = _parse_dashes(line)
            elif "coord atom" in ll:
                bi.is_hp        = True
                bi.data_offset  = fh.tell()
                break
            # Standard: Headerzeile "Atom  AN      X      Y      Z"
            elif "atom" in ll and (" an" in ll or "an " in ll):
                bi.is_hp        = False
                bi.data_offset  = fh.tell()
                break
            # Standard: erste data row direkt (kein Header)
            elif re.match(r"\s+\d+\s+\d+\s+[-\d.]", line):
                bi.is_hp        = False
                bi.data_offset  = prev_pos   # B4-Fix
                break

    for lst in (bi.red_masses, bi.frc_consts):
        while len(lst) < n:
            lst.append(0.0)
    while len(bi.syms) < n:
        bi.syms.append("A")

    return bi


def read_all_meta(
        filepath: str,
        nc_offsets: List[int],
        fr_offsets: List[int],
        runlog: Optional["RunLog"] = None,
) -> Tuple[List[BlockInfo], Dict[int, BlockInfo], Dict[int, List]]:
    """Reads Metadaten aller Frequenz-Gruppen in the Logdatei.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    nc_offsets : list of int
        Byte-Offsets aller ``and normal coordinates:``-Zeilen.
    fr_offsets : list of int
        Byte-Offsets aller ``Frequencies ---``-Zeilen.

    Returns
    -------
    blocks : list of BlockInfo
        Alle read inen Frequenz-blocks.
    best_block : dict of {int: BlockInfo}
        Mapping Modenummer → bester Block (HP beforezugt).
    cand_map : dict of {int: list of BlockInfo}
        Mapping Modenummer → alle Kandidaten-blocks.
    """
    if not fr_offsets:
        return [], {}, {}

    nc_sorted = sorted(nc_offsets)
    first_nc  = nc_sorted[0] if nc_sorted else 0
    fr_in     = sorted(f for f in fr_offsets if f > first_nc)

    t0 = time.time()
    n  = len(fr_in)
    raw_blocks: List[BlockInfo] = []
    for i, fr_off in enumerate(fr_in):
        bi = _read_block_at_freq_line(filepath, fr_off)
        if bi and bi.freqs:
            raw_blocks.append(bi)
        if (i + 1) % 500 == 0 or (i + 1) == n:
            if (i + 1) == n or (i + 1) % 500 == 0:
                _prog = f"Metadata: {i+1}/{n} ({(i+1)/n*100:.0f}%)"
                if runlog is not None:
                    runlog.info(_prog)
            print(f"    Metadata: {i+1}/{n} "
                  f"({(i+1)/n*100:.0f}%)    ", end="\r")

    _blk_msg = f"Blocks read: {len(raw_blocks)} in {time.time()-t0:.1f} s"
    if runlog is not None:
        runlog.info(_blk_msg)
    print(f"\n    {_blk_msg}")

    def _which_section(off: int) -> int:
        """Determines the naechstliegenden Normcoord-Offset for einen file-Offset."""
        sec = nc_sorted[0]
        for nc in nc_sorted:
            if nc <= off:
                sec = nc
        return sec

    section_count: Dict[int, int] = {}
    blocks: List[BlockInfo] = []
    for bi in raw_blocks:
        sec   = _which_section(bi.offset)
        section_count.setdefault(sec, 0)
        nm    = len(bi.mode_nums)
        fallback = (bi.mode_nums == list(range(1, nm + 1))
                    and section_count[sec] > 0)
        if fallback or not bi.mode_nums or bi.mode_nums[0] <= 0:
            start = section_count[sec] * nm + 1
            bi.mode_nums = list(range(start, start + nm))
        section_count[sec] += 1
        blocks.append(bi)

    best:     Dict[int, BlockInfo] = {}
    cand_map: Dict[int, List]      = {}
    for bi in blocks:
        for mn in bi.mode_nums:
            if mn not in best or (bi.is_hp and not best[mn].is_hp):
                best[mn] = bi
            cand_map.setdefault(mn, [])
            if bi not in cand_map[mn]:
                cand_map[mn].append(bi)

    n_hp  = sum(1 for bi in blocks if bi.is_hp)
    n_std = len(blocks) - n_hp
    _meta_msg = f"HP groups: {n_hp}  |  default groups: {n_std}"
    if runlog is not None:
        runlog.info(_meta_msg)
    print(f"    {_meta_msg}")
    return blocks, best, cand_map


def check_hp_std_frequency_consistency(
        cand_map: Dict[int, List["BlockInfo"]],
        tolerance_cm1: float = 0.01,
) -> Dict[str, object]:
    r"""Vergleicht the Frequenzen between HP- and Standard-blocksn per Mode.

    For a Gaussian-Job with ``freq=hpmodes`` are for each Mode zwei
    blocks ausgiven: the Standard-Block (displacementen with 2 Dezimal-
    stellen) and the HP-Block (5 Dezimalstellen). Die *Frequenzen* selbst
    are in beiden blocksn Ergebnis derselben Hessian-Diagonalisierung
    and muessen exakt uebereinstimmen — abgesehen von the Druckpraezision
    the Standard-Ausgabe (~ 0.005 cm^-1).

    Diese Funktion vergleicht the Frequenzen for each Mode, the in beiden
    Block-Typen vorkommt, and meldet:

      * the maximale absolute Abweichung
      * the Number of modes, for denen the Abweichung > ``tolerance_cm1`` ist
      * the TOP-5-Ausreisser to Diagnose

    A relevante Diskrepanz (groesser als the Druckpraezision) deutet auf
    Parser-Versatz, korrupte Logdatei or zwei verkettete freq-Jobs hin
    and sollte before weiterer Analyse ueberprueft werden.

    Parameters
    ----------
    cand_map : dict of {int: list of BlockInfo}
        Aus :func:`read_all_meta`. Modennummer -> alle Kandidaten-blocks.
    tolerance_cm1 : float, optional
        Threshold for "diskrepant" in cm^-1. default 0.01 (= zweifache
        Druckpraezision the Standard-Ausgabe).

    Returns
    -------
    dict
        key:
          ``"checked_modes"`` : int — Number of modes with beiden Block-Typen
          ``"max_dev_cm1"``   : float — groesste |HP - std| Differenz
          ``"mean_dev_cm1"``  : float — mittlere |HP - std| Differenz
          ``"n_outliers"``    : int — modes with |dev| > tolerance
          ``"outliers"``      : list of (mode, freq_hp, freq_std, dev) — Top 5
          ``"any_hp"``        : bool — ueberhaupt HP-blocks present?
          ``"any_std"``       : bool — ueberhaupt Standard-blocks present?
    """
    devs: List[Tuple[int, float, float, float]] = []  # (mode, f_hp, f_std, dev)
    any_hp = False
    any_std = False

    for mn, cands in cand_map.items():
        hp_freq: Optional[float] = None
        std_freq: Optional[float] = None
        for bi in cands:
            try:
                col = bi.mode_nums.index(mn)
            except ValueError:
                continue
            if col >= len(bi.freqs):
                continue
            f = float(bi.freqs[col])
            if bi.is_hp:
                if hp_freq is None:    # erster Treffer reicht
                    hp_freq = f
                any_hp = True
            else:
                if std_freq is None:
                    std_freq = f
                any_std = True
        if hp_freq is not None and std_freq is not None:
            dev = abs(hp_freq - std_freq)
            devs.append((mn, hp_freq, std_freq, dev))

    n_checked = len(devs)
    if n_checked == 0:
        return {
            "checked_modes": 0,
            "max_dev_cm1":   0.0,
            "mean_dev_cm1":  0.0,
            "n_outliers":    0,
            "outliers":      [],
            "any_hp":        any_hp,
            "any_std":       any_std,
        }

    abs_devs = [d[3] for d in devs]
    max_dev = max(abs_devs)
    mean_dev = sum(abs_devs) / len(abs_devs)
    outliers = [d for d in devs if d[3] > tolerance_cm1]
    outliers.sort(key=lambda x: -x[3])
    top5 = outliers[:5]
    return {
        "checked_modes": n_checked,
        "max_dev_cm1":   float(max_dev),
        "mean_dev_cm1":  float(mean_dev),
        "n_outliers":    len(outliers),
        "outliers":      top5,
        "any_hp":        any_hp,
        "any_std":       any_std,
    }


# ===========================================================================
# Phase 4 - Eigenvektoren
# ===========================================================================

def read_hp_eigvec(
        filepath: str,
        bi: BlockInfo,
        col: int,
        include_hydrogen: bool,
) -> Tuple[List[int], np.ndarray]:
    """Reads einen Eigenvector from a HP-Format-Block.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    bi : BlockInfo
        Metadaten of the Ziel-Blocks.
    col : int
        Spalten-Index (0-basiert) the gesuchten Mode in the Block.
    include_hydrogen : bool
        hydrogen-Atome einlesen.

    Returns
    -------
    centers : list of int
        Gaussian-Center-Nummern the gelesenen Atome.
    evg : ndarray of shape (n_atoms, 3)
        Eigenvektor-Komponenten (dx, dy, dz) per Atom.

    Notes
    -----
    Bugfix B12: ``ci_idx < 1`` or ``ci_idx > 3`` bewirkt ein ``break``
    statt ``xyz = -1``.
    """
    centers: List[int] = []
    dxyz: Dict[int, List[float]] = {}
    pat = re.compile(
        rb"\s*(\d+)\s+(\d+)\s+(\d+)((?:\s+[-\d.]+)+)")

    with open(filepath, "rb") as fh:
        if bi.data_offset < 0:
            raise ValueError(
                f"data_offset invalid ({bi.data_offset}) "
                f"for block @ {bi.offset}")
        fh.seek(bi.data_offset)
        for line in fh:
            m = pat.match(line)
            if not m:
                break
            ci_idx = int(m.group(1))
            if ci_idx < 1 or ci_idx > 3:   # B12
                break
            ctr  = int(m.group(2))
            anum = int(m.group(3))
            if anum == 1 and not include_hydrogen:
                continue
            xyz = ci_idx - 1
            try:
                val = float(m.group(4).split()[col])
            except (IndexError, ValueError):
                continue
            if ctr not in dxyz:
                dxyz[ctr]  = [0.0, 0.0, 0.0]
                centers.append(ctr)
            dxyz[ctr][xyz] = val

    if not centers:
        return [], np.zeros((0, 3))
    return centers, np.array([dxyz[c] for c in centers], dtype=float)


def read_std_eigvec(
        filepath: str,
        bi: BlockInfo,
        col: int,
        include_hydrogen: bool,
) -> Tuple[List[int], np.ndarray]:
    """Reads einen Eigenvector from a Standard-Format-Block.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    bi : BlockInfo
        Metadaten of the Ziel-Blocks.
    col : int
        Spalten-Index (0-basiert) the gesuchten Mode in the Block.
    include_hydrogen : bool
        hydrogen-Atome einlesen.

    Returns
    -------
    centers : list of int
        Gaussian-Center-Nummern the gelesenen Atome.
    evg : ndarray of shape (n_atoms, 3)
        Eigenvektor-Komponenten (dx, dy, dz) per Atom.

    Notes
    -----
    Bugfix B9: Flexible Regex; ``continue`` statt ``break`` for kurzen Zeilen.
    Bugfix B11: Doppelte Center are per ``dxyz``-Dict dedupliziert.
    """
    centers: List[int] = []
    dxyz: Dict[int, List[float]] = {}
    nm  = len(bi.mode_nums)
    pat = re.compile(rb"\s*(\d+)\s+(\d+)((?:\s+[-\d.]+)+)")

    with open(filepath, "rb") as fh:
        if bi.data_offset < 0:
            raise ValueError(
                f"data_offset invalid ({bi.data_offset}) "
                f"for block @ {bi.offset}")
        fh.seek(bi.data_offset)
        for line in fh:
            m = pat.match(line)
            if not m:
                break
            ctr  = int(m.group(1))
            anum = int(m.group(2))
            if anum == 1 and not include_hydrogen:
                continue
            vals = [float(x) for x in m.group(3).split()]
            if len(vals) < nm * 3:
                continue   # B9: to kurze Zeile skip
            dx = vals[col * 3 + 0]
            dy = vals[col * 3 + 1]
            dz = vals[col * 3 + 2]
            if ctr not in dxyz:
                centers.append(ctr)
            dxyz[ctr] = [dx, dy, dz]   # B11: no Duplikat

    if not centers:
        return [], np.zeros((0, 3))
    return centers, np.array([dxyz[c] for c in centers], dtype=float)


def get_eigvec(
        filepath: str,
        bi: BlockInfo,
        col: int,
        include_hydrogen: bool,
) -> Tuple[List[int], np.ndarray]:
    """Reads the eigenvector of a single mode.

    Selects automatically HP or standard reader and validates
    ``data_offset``.

    Parameters
    ----------
    filepath : str
        Path to the Gaussian ``.log``-file.
    bi : BlockInfo
        Metadaten of the Ziel-Blocks.
    col : int
        Spalten-Index (0-basiert) the gesuchten Mode in the Block.
    include_hydrogen : bool
        hydrogen-Atome einlesen.

    Returns
    -------
    centers : list of int
        Gaussian-Center-Nummern the gelesenen Atome.
    evg : ndarray of shape (n_atoms, 3)
        Eigenvektor-Komponenten (dx, dy, dz) per Atom.

    Raises
    ------
    ValueError
        If ``bi.data_offset < 0`` (Block-Format not erkannt).
    """
    if bi.data_offset < 0:
        raise ValueError(
            f"data_offset not determined for block @ {bi.offset}. "
            f"Eigenvektor-Format not erkannt.")
    if bi.is_hp:
        return read_hp_eigvec(filepath, bi, col, include_hydrogen)
    return read_std_eigvec(filepath, bi, col, include_hydrogen)


# ===========================================================================
# PDB-Parser
# ===========================================================================

def _is_hydrogen(aname: str, elem_raw: str) -> bool:
    """Checks whether an atom is hydrogen."""
    if elem_raw.strip().upper() in ("H", "D"):
        return True
    s = aname.strip()
    if s.startswith(("H", "D")):
        return True
    if s and s[0].isdigit() and "H" in s.upper():
        return True
    return False


def _element_from_name(aname: str) -> str:
    """Leitet the elementsymbol from the PDB atomnamen ab.

    Parameters
    ----------
    aname : str
        PDB atomname (4 Zeichen, e.g. ``" FE "``).

    Returns
    -------
    str
        Grossbuchstaben-elementsymbol (e.g. ``"FE"``).

    Notes
    -----
    Bugfix B13: Zwei-Zeichen-elements (``FE``, ``ZN``, ``MG``, …)
    are korrekt erkannt.
    """
    s = aname.lstrip("0123456789 ").rstrip("0123456789 ")
    if len(s) >= 2 and s[:2].upper() in _TWO_CHAR_ELEM:
        return s[:2].upper()
    return s[0].upper() if s else "?"


def parse_pdb(pdb_path: str, chain_filter: str = "") -> Dict:
    """Reads a PDB-file complete ein.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB-file.
    chain_filter : str, optional
        Ketten-ID for ATOM-Records (e.g. ``"A"``). Leer = alle Ketten.
        HETATM-Records are immer geladen (FES kann in anderer Kette
        liegen, Bugfix B8).

    Returns
    -------
    dict
        Enthaelt folgende Schlussel:

        * ``atoms``    - ATOM-Records without H
        * ``atoms_h``  - ATOM-Records with H
        * ``hetatm``   - HETATM-Records without H
        * ``hetatm_h`` - HETATM-Records with H
        * ``all_h``    - ATOM + HETATM with H (fuer Kabsch + Koordination)
        * ``sse_elements`` - HELIX- and SHEET-Records als List of Dicts

    Notes
    -----
    Bugfix B8: Chain-Filter gilt only for ATOM-Records. HETATM (FES-Cluster)
    is immer read in, so that the Kabsch-Alignment also dann funktioniert,
    if FES in a anderen Kette liegt als the Protein.

    Bugfix B14: Nur the erste MODEL-Record is processed.
    """
    atoms:    List[Dict] = []
    atoms_h:  List[Dict] = []
    hetatm:   List[Dict] = []
    hetatm_h: List[Dict] = []
    sse_elems: List[Dict] = []

    in_model   = False
    model_done = False

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()

            # B14: only erstes Modell
            if rec == "MODEL":
                if model_done:
                    break
                in_model = True
                continue
            if rec == "ENDMDL":
                model_done = True; in_model = False; continue
            if rec in ("END", "TER"):
                if in_model:
                    model_done = True
                continue

            # HELIX / SHEET
            if rec == "HELIX":
                try:
                    ch   = line[19] if len(line) > 19 else "A"
                    r1   = int(line[21:25])
                    r2   = int(line[33:37])
                    if r2 - r1 + 1 >= 4:
                        sse_elems.append({
                            "type": "helix", "chain": ch,
                            "res_start": r1, "res_end": r2,
                            "name": f"Helix_{ch}_{r1}_{r2}"})
                except (ValueError, IndexError):
                    pass
                continue
            if rec == "SHEET":
                try:
                    ch   = line[21] if len(line) > 21 else "A"
                    r1   = int(line[22:26])
                    r2   = int(line[33:37])
                    sid  = line[11:14].strip()
                    if r2 - r1 + 1 >= 4:
                        sse_elems.append({
                            "type": "sheet", "chain": ch,
                            "res_start": r1, "res_end": r2,
                            "name": f"Sheet_{sid}_{ch}_{r1}"})
                except (ValueError, IndexError):
                    pass
                continue

            is_atom   = rec == "ATOM"
            is_hetatm = rec == "HETATM"
            if not (is_atom or is_hetatm):
                continue

            alt = line[16] if len(line) > 16 else " "
            if alt not in (" ", "A"):
                continue

            aname = line[12:16].strip() if len(line) > 15 else ""
            rname = line[17:20].strip() if len(line) > 19 else ""
            chain = line[21]            if len(line) > 21 else " "

            # B8: Chain-Filter NUR for ATOM-Records.
            # v3.7.5+: additionally are FES-Cluster-Atome (Fe and das
            # FES-Residue with Br\"ucken-S) immer durchgelassen, also wenn
            # sie als ATOM-Record (statt HETATM) markiert are and in einer
            # nicht-passenden Kette stehen. Ohne diese Ausnahme schluckt der
            # Default-pdb_chain='A'-Filter z.\,B.\ ORCA-QM/MM-Files, in denen
            # the cluster als ATOM in Kette B steht.
            aname_pre = line[12:16].strip().upper() if len(line) > 15 else ""
            elem_pre  = line[76:78].strip().upper() if len(line) > 77 else ""
            is_cluster_atom = (
                rname.upper() == "FES" or
                elem_pre == "FE" or
                aname_pre.startswith("FE")
            )
            if (is_atom and chain_filter and chain != chain_filter
                    and not is_cluster_atom):
                continue

            try:
                rnum = int(line[22:26])
            except ValueError:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                x = y = z = 0.0

            elem_raw = line[76:78].strip() if len(line) > 76 else ""
            elem     = elem_raw.upper() if elem_raw else _element_from_name(aname)
            is_h     = _is_hydrogen(aname, elem_raw)

            entry: Dict = {
                "record": rec,  "aname": aname, "rname": rname,
                "chain":  chain,"rnum":  rnum,
                "x": x,  "y": y, "z": z,
                "element": elem, "is_h": is_h,
            }

            if is_atom:
                atoms_h.append(entry)
                if not is_h:
                    atoms.append(entry)
            else:
                hetatm_h.append(entry)
                if not is_h:
                    hetatm.append(entry)

    all_h = (
        [dict(a, _list="atom",   _list_idx=i) for i, a in enumerate(atoms_h)] +
        [dict(a, _list="hetatm", _list_idx=i) for i, a in enumerate(hetatm_h)]
    )

    return {
        "atoms":       atoms,
        "atoms_h":     atoms_h,
        "hetatm":      hetatm,
        "hetatm_h":    hetatm_h,
        "all_h":       all_h,
        "sse_elements": sse_elems,
    }

__version__ = "1.4"  # modenanalyse v1.4
