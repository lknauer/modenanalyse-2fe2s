# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

"""
geometry.py
========================
Geometrie-Analyse for [2Fe-2S]-Cluster and Proteinstruktur.

Oeffentliche functions
-----------------------
find_cluster
    Identifiziert the [2Fe-2S]-Cluster in the Gaussian-atom list.
cluster_normal
    Computes the normal vector the clusterebene.
compute_dist_ref
    Computes Gleichgewichtsabstaende in the Cluster.
kabsch_align
    Kabsch-Alignment PDB → Gaussian-Koordinaten.
find_coordinating_residues
    Erkennt koordinierende amino acidn automatically from the PDB-Geometrie.
build_ss_center_map
    Creates the Mapping SS-element → Gaussian-Center-Nummern.
get_calpha_centers
    Returns Gaussian-Center-Nummern aller Cα-Atome .

classes
-------
LigandInfo
    Informationen over einen Fe-ligands-Kontakt.
CoordInfo
    Vollstaendige Koordinations-Information of the cluster.

Bugfixes (gegenvia Vorversion)
---------------------------------
B2  SS-Analyse:          PDB listnindizes ≠ Gaussian-Center-Reihenfolge.
                         Behoben through ``pdb_to_center``-Mapping
                         (PDB-Index → Gaussian-Center).
B7  Fe-N/S/O-Analyse:    Das koordinierende Fe is jetzt aus
                         ``LigandInfo.fe_center`` entnommen statt immer
                         ``fe_c[0]`` to verwenden.
"""
from __future__ import annotations

import math
import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .config import Config
from .logio import RunLog, _HIS, _CYS, parse_pdb



# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def _dist(a: Dict, b: Dict) -> float:
    """Euklidischer Abstand zweier Atom-Dicts."""
    return math.sqrt(
        (a["x"] - b["x"])**2 +
        (a["y"] - b["y"])**2 +
        (a["z"] - b["z"])**2)


def _pos(a: Dict) -> np.ndarray:
    """Koordinaten-Array a Atom-Dicts."""
    return np.array([a["x"], a["y"], a["z"]])


# ===========================================================================
# Cluster-Erkennung
# ===========================================================================

def find_all_clusters(
        atoms: List[Dict],
        cfg: Config,
) -> List[Tuple[List[int], List[int], Dict[str, float]]]:
    """Finds ALLE [2Fe-2S]-Cluster in the atom list (Hardening v3.0 #10).

    Identifiziert mehrere Cluster in Multi-Iron-Systemen (Dimere, Proteine
    with zwei or mehr [2Fe-2S]-Clustern wie Glutaredoxin-Dimer etc.).
    Verwendet Greedy-Assignment: jedes Atom kann only to a Cluster
    gehoeren, so that z. B. ein einzelnes Fe not in zwei differenten
    Clustern auftaucht.

    Algorithmus:
      1. Sortiere alle Fe-Fe-Paare below ``fe_fe_cutoff`` after Abstand.
      2. Greedy: the engste Paar is to ersten Cluster. Beide Fe atoms
         are als verbraucht markiert.
      3. Suche zwei brueckende S-Atome, the to beiden Fe < ``fe_s_cutoff``
         haben and still not verbraucht sind.
      4. Wiederhole 2-3 with the verbleibenthe atoms, bis no Fe-Fe-Paar
         mehr in Reichweite ist.

    Parameters
    ----------
    atoms : list of dict
        Gaussian-atom list (ohne hydrogen).
    cfg : Config
        Konfiguration; required ``fe_fe_cutoff`` and ``fe_s_cutoff``.

    Returns
    -------
    clusters : list of (fe_centers, s_centers, geom_info)
        Pro Cluster ein 3-Tuple:

          * ``fe_centers`` — Center-Nummern the zwei Fe atoms
          * ``s_centers``  — Center-Nummern the zwei brueckenden S-Atome
          * ``geom_info``  — Dict with ``fe_fe``, ``fe_s_min``, ``fe_s_max``
            in Angstrom (to Diagnostik / Auswahl)

        Sortierung: after Fe-Fe-Abstand (engstes Paar zuerst).
    """
    fe_atoms = [a for a in atoms if a["atomic_num"] == 26]
    s_atoms  = [a for a in atoms if a["atomic_num"] == 16]

    # Alle Fe-Fe-Paare below Cutoff, after Abstand sortiert.
    # Rundungs-stabiler Tiebreaker: Distanzen are on 6 Dezimalen
    # (10^-6 A = 1 fm) gerundet before sortiert wird; for numerisch
    # gleichen Distanzen entscheidet the Atom-Index-Reihenfolge (i, j).
    # So that is the Cluster-Nummerierung unabhaengig von Floating-Point-
    # Rundung deterministisch.
    pairs = sorted(
        [(i, j, _dist(fe_atoms[i], fe_atoms[j]))
         for i, j in combinations(range(len(fe_atoms)), 2)
         if _dist(fe_atoms[i], fe_atoms[j]) < cfg.fe_fe_cutoff],
        key=lambda x: (round(x[2], 6), x[0], x[1]),
    )

    clusters: List[Tuple[List[int], List[int], Dict[str, float]]] = []
    used_fe: set = set()
    used_s:  set = set()

    for fi, fj, d_fefe in pairs:
        if fi in used_fe or fj in used_fe:
            continue
        fe2 = [fe_atoms[fi], fe_atoms[fj]]

        # Brueckende S-Atome: nahe BEIDEN Fe and not schon vergeben
        cl_s = []
        for si, s in enumerate(s_atoms):
            if si in used_s:
                continue
            d1, d2 = _dist(s, fe2[0]), _dist(s, fe2[1])
            if d1 < cfg.fe_s_cutoff and d2 < cfg.fe_s_cutoff:
                cl_s.append((si, s, max(d1, d2)))
        cl_s.sort(key=lambda x: x[2])

        if len(cl_s) >= 2:
            si1, s1, d_fs1 = cl_s[0]
            si2, s2, d_fs2 = cl_s[1]
            cl_s_atoms = [s1, s2]
            used_s.update([si1, si2])
        elif len(cl_s) == 1:
            # Nur ein brueckendes S — Fallback: zweites S als naechst-
            # liegendes also for groesserem Abstand suchen
            si1, s1, d_fs1 = cl_s[0]
            remaining = [(si, s, min(_dist(s, fe2[0]), _dist(s, fe2[1])))
                         for si, s in enumerate(s_atoms)
                         if si != si1 and si not in used_s]
            if remaining:
                remaining.sort(key=lambda x: x[2])
                si2, s2, _ = remaining[0]
                cl_s_atoms = [s1, s2]
                used_s.update([si1, si2])
                d_fs2 = max(_dist(s2, fe2[0]), _dist(s2, fe2[1]))
            else:
                continue   # not genug S-Atome - dieses Paar discard
        else:
            continue       # gar no S in the Naehe - no Cluster

        used_fe.update([fi, fj])
        geom_info = {
            "fe_fe":    d_fefe,
            "fe_s_min": min(d_fs1, d_fs2),
            "fe_s_max": max(d_fs1, d_fs2),
        }
        clusters.append((
            [a["center"] for a in fe2],
            [a["center"] for a in cl_s_atoms],
            geom_info,
        ))

    return clusters


def find_cluster(
        atoms: List[Dict],
        cfg: Config,
) -> Tuple[List[int], List[int]]:
    """Identifiziert the primaere [2Fe-2S]-Cluster in the atom list.

    RunLogft :func:`find_all_clusters`, um alle presenten Cluster zu
    identifizieren, and liefert the through ``cfg.cluster_index`` auschosene
    Cluster (default 0 = engstes Fe-Fe-Paar). For Multi-Cluster-Systemen
    (z. B. Glutaredoxin-Dimer, e.g. ein Dimer with zwei [2Fe-2S]) gibt the einen
    nachvollziehbaren Mechanismus, the gewuenschten Cluster to choose,
    statt a willkuerliche Wahl to treffen.

    Parameters
    ----------
    atoms : list of dict
        Gaussian-atom list (ohne hydrogen).
    cfg : Config
        Konfiguration; required ``fe_fe_cutoff``, ``fe_s_cutoff`` und
        optional ``cluster_index``.

    Returns
    -------
    fe_c : list of int
        Gaussian-Center-Nummern the zwei Fe atoms of the chosenen Clusters.
    s_c : list of int
        Gaussian-Center-Nummern the zwei brueckenden S-Atome.

    Raises
    ------
    ValueError
        If no Cluster found wurde, or ``cluster_index`` ausserhalb
        of the found Cluster liegt.

    Notes
    -----
    Anschliessend geometrische Plausibilitaetspruefung with UserWarning bei
    Auffaelligkeiten (kein Programmabbruch).
    """
    clusters = find_all_clusters(atoms, cfg)
    n_total = len(clusters)

    if n_total == 0:
        n_fe = sum(1 for a in atoms if a["atomic_num"] == 26)
        raise ValueError(
            f"No [2Fe-2S] cluster found "
            f"(Fe-Fe-Cutoff = {cfg.fe_fe_cutoff} A, "
            f"Fe-S cutoff = {cfg.fe_s_cutoff} A; {n_fe} Fe atoms in system).")

    idx = getattr(cfg, "cluster_index", 0)
    if not (0 <= idx < n_total):
        # Diagnostik presenter Cluster for the User
        diag = "; ".join(
            f"#{i}: Fe-Fe={c[2]['fe_fe']:.2f}A, Fe-S={c[2]['fe_s_min']:.2f}-"
            f"{c[2]['fe_s_max']:.2f}A"
            for i, c in enumerate(clusters))
        raise ValueError(
            f"cluster_index={idx} outside the {n_total} found "
            f"Cluster. Verfuegbar: {diag}")

    # Multi-cluster-Hinweis (UserWarning, is in the main script in den
    # Runlog umgeleitet)
    if n_total > 1:
        import warnings as _w
        diag = "; ".join(
            f"#{i}: Fe-Fe={c[2]['fe_fe']:.2f}A"
            for i, c in enumerate(clusters))
        _w.warn(
            f"find_cluster: {n_total} clusters found in system ({diag}). "
            f"Analysiert is Cluster #{idx} (cluster_index={idx}). "
            f"Fuer einen anderen Cluster: cluster_index in the Konfiguration "
            f"setzen.",
            UserWarning, stacklevel=2,
        )

    fe_c, s_c, _geom = clusters[idx]

    # Geometrische Plausibilitaetspruefung
    fe_atoms = [a for a in atoms if a["atomic_num"] == 26]
    s_atoms  = [a for a in atoms if a["atomic_num"] == 16]
    fe2 = [a for a in fe_atoms if a["center"] in fe_c]
    s2  = [a for a in s_atoms  if a["center"] in s_c]
    _check_cluster_geometry(fe2, s2)

    return fe_c, s_c


def _check_cluster_geometry(
        fe2: List[Dict],
        s2:  List[Dict],
) -> None:
    """Checks whether the erkannte Cluster geometrisch plausibel ist.

    Returns UserWarnings from for Auffaelligkeiten — no Programmabbruch.

    Typische [2Fe-2S]-Referenzwerte (Literatur):
    - Fe-Fe:   2.69–2.77 Å (Rieske), 2.72–2.80 Å (Ferredoxin)
    - Fe-S:    2.15–2.35 Å
    - Fe-S-Fe: 70–80°
    """
    import warnings as _w
    msgs: List[str] = []

    if len(fe2) < 2 or len(s2) < 2:
        _w.warn(f"find_cluster: uncompleteer Kern ({len(fe2)} Fe, {len(s2)} S)",
                UserWarning, stacklevel=4)
        return

    d_fefe = _dist(fe2[0], fe2[1])
    if not (2.40 < d_fefe < 3.10):
        msgs.append(
            f"Fe-Fe = {d_fefe:.3f} Å (typical 2.69–2.80 Å for [2Fe-2S])")

    for si, s in enumerate(s2):
        for fi, fe in enumerate(fe2):
            d = _dist(s, fe)
            if d > 2.60:
                msgs.append(
                    f"Fe{fi+1}-S{si+1} = {d:.3f} Å > 2.60 Å: "
                    f"kein Brueckensulfid?")
            elif not (2.00 < d < 2.60):
                msgs.append(
                    f"Fe{fi+1}-S{si+1} = {d:.3f} Å (typical 2.15–2.35 Å)")

    # Fe-S-Fe-Winkel
    p1 = _pos(fe2[0]); p2 = _pos(fe2[1])
    for si, s in enumerate(s2):
        ps = _pos(s)
        v1 = p1 - ps; v2 = p2 - ps
        n1 = np.linalg.norm(v1); n2 = np.linalg.norm(v2)
        if n1 > 0 and n2 > 0:
            ang = float(np.degrees(
                np.arccos(np.clip(np.dot(v1, v2)/(n1*n2), -1., 1.))))
            if not (60. < ang < 95.):
                msgs.append(
                    f"Fe-S{si+1}-Fe = {ang:.1f}° (typical 70–80°)")

    for msg in msgs:
        _w.warn(f"find_cluster Geometrie: {msg}", UserWarning, stacklevel=4)




def cluster_normal(
        atoms: List[Dict],
        idx_map: Dict[int, int],
        fe_c: List[int],
        s_c: List[int],
) -> np.ndarray:
    """Computes the normal vector the Fe-S-Clusterebene.

    Parameters
    ----------
    atoms : list of dict
        Gaussian-atom list.
    idx_map : dict of {int: int}
        Center → Atom-Index.
    fe_c : list of int
        Center-Nummern the Fe atoms.
    s_c : list of int
        Center-Nummern the S-Atome.

    Returns
    -------
    ndarray of shape (3,)
        Einheits-normal vector (kleinster Singularvektor the SVD).
        Fallback: ``[0, 0, 1]`` if weniger als 3 Atome present.
    """
    pts = np.array([
        [atoms[idx_map[c]]["x"], atoms[idx_map[c]]["y"], atoms[idx_map[c]]["z"]]
        for c in fe_c + s_c if c in idx_map
    ])
    if pts.shape[0] < 3:
        return np.array([0., 0., 1.])
    ctr = pts.mean(0)
    _, _, Vt = np.linalg.svd(pts - ctr)
    return Vt[-1]


def compute_dist_ref(
        atoms: List[Dict],
        idx_map: Dict[int, int],
        fe_c: List[int],
        s_c: List[int],
        cfg: Config,
) -> Dict[str, Tuple[float, float]]:
    """Computes Gleichgewichtsabstaende in the Cluster with Unsicherheiten.

    Parameters
    ----------
    atoms : list of dict
        Gaussian-atom list.
    idx_map : dict of {int: int}
        Center → Atom-Index.
    fe_c : list of int
        Center-Nummern the Fe atoms.
    s_c : list of int
        Center-Nummern the S-Atome.
    cfg : Config
        Konfiguration; required ``sigma_coord``.

    Returns
    -------
    dict of {str: (float, float)}
        Mapping Abstandsname → ``(value_in_Angstrom, Sigma_in_Angstrom)``.
        Schlussel: ``"Fe-Fe"``, ``"Fe1-S1"``, ``"Fe2-S1"``, etc.
    """
    s_coord  = float(np.sqrt(2) * cfg.sigma_coord)
    dist_ref: Dict[str, Tuple[float, float]] = {}

    def _add(key: str, c1: int, c2: int) -> None:
        """Computes and speichert einen Clusterabstand if beide Center in the Index sind."""
        if c1 in idx_map and c2 in idx_map:
            d = float(np.linalg.norm(_pos(atoms[idx_map[c1]])
                                     - _pos(atoms[idx_map[c2]])))
            dist_ref[key] = (d, s_coord)

    if len(fe_c) >= 2:
        _add("Fe-Fe", fe_c[0], fe_c[1])
    for si, sc in enumerate(s_c[:2]):
        for fi, fc in enumerate(fe_c[:2]):
            _add(f"Fe{fi+1}-S{si+1}", fc, sc)

    return dist_ref


# ===========================================================================
# Kabsch-Alignment
# ===========================================================================

def kabsch_align(
        pdb_data: Dict,
        atoms: List[Dict],
        idx_map: Dict[int, int],
        fe_c: List[int],
        s_c: List[int],
        runlog: RunLog,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], float, bool]:
    """Kabsch-Alignment: PDB-Koordinaten → Gaussian-Koordinaten.

    Verwendet the vier Cluster-Anker (Fe1, Fe2, S1, S2) als Referenz.

    Parameters
    ----------
    pdb_data : dict
        Ausgabe von ``parse_pdb``; required ``all_h``.
    atoms : list of dict
        Gaussian-atom list (ohne hydrogen).
    idx_map : dict of {int: int}
        Center → Atom-Index.
    fe_c : list of int
        Center-Nummern the Gaussian-Fe atoms.
    s_c : list of int
        Center-Nummern the Gaussian-S-Atome.
    runlog : RunLog
        Fuer Warnmeldungen.

    Returns
    -------
    R : ndarray of shape (3, 3) or None
        Rotationsmatrix. ``None`` for Error.
    t : ndarray of shape (3,) or None
        Translationsvektor. ``None`` for Error.
    rmsd : float
        Kabsch-RMSD in Angstrom. ``nan`` for Error.
    ok : bool
        ``True`` if the Alignment erfolgreich war.

    Notes
    -----
    Transformation: ``gaus_xyz = R @ pdb_xyz + t``.
    """
    if len(fe_c) < 2 or len(s_c) < 2:
        return None, None, float("nan"), False

    log_cl: List[List[float]] = []
    for c in (fe_c + s_c)[:4]:
        if c not in idx_map:
            return None, None, float("nan"), False
        a = atoms[idx_map[c]]
        log_cl.append([a["x"], a["y"], a["z"]])
    log_cl_arr = np.array(log_cl, dtype=float)

    # Gaussian-Fe-Fe-Distanz dieses Clusters (Multi-Cluster-Anker)
    log_fefe = float(np.linalg.norm(log_cl_arr[0] - log_cl_arr[1]))

    # PDB-Fe atoms
    pdb_fe = [(a["x"], a["y"], a["z"]) for a in pdb_data["all_h"]
              if a["element"] == "FE"]
    if not pdb_fe:
        pdb_fe = [(a["x"], a["y"], a["z"]) for a in pdb_data["all_h"]
                  if a["aname"].upper().startswith("FE")]
    if len(pdb_fe) < 2:
        runlog.warn("Kabsch: Less than 2 Fe atoms in PDB.")
        return None, None, float("nan"), False

    # Multi-cluster-Fix: for mehreren Fe atomsn is dasjenige Fe-Paar
    # chosen, dessen Fe-Fe-Distanz at the naechsten an the Gaussian-Fe-Fe-
    # Distanz liegt UND the eng beieinander liegt (Cluster-Distanz < 4 A,
    # over ``fe_fe_cutoff`` hinaus with Sicherheits-Marge). For only 2 Fe
    # in PDB faellt the on the alte Logik .
    if len(pdb_fe) > 2:
        from itertools import combinations
        best_pair: Optional[Tuple[int, int]] = None
        best_score = float("inf")
        for i, j in combinations(range(len(pdb_fe)), 2):
            d = float(np.linalg.norm(np.array(pdb_fe[i]) - np.array(pdb_fe[j])))
            # Nur Paare with plausibler Cluster-Distanz (< 4 A)
            if d > 4.0:
                continue
            # Score: Abweichung von the Gaussian-Fe-Fe-Distanz
            score = abs(d - log_fefe)
            if score < best_score:
                best_score = score
                best_pair = (i, j)
        if best_pair is None:
            runlog.warn(f"Kabsch: no PDB Fe pair with cluster distance "
                        f"< 4.0 A found ({len(pdb_fe)} Fe atoms). "
                        f"Using first two Fe.")
            fe1p = np.array(pdb_fe[0]); fe2p = np.array(pdb_fe[1])
        else:
            i, j = best_pair
            fe1p = np.array(pdb_fe[i]); fe2p = np.array(pdb_fe[j])
            runlog.info(f"Kabsch: PDB Fe pair ({i}, {j}) chosen, "
                        f"Fe-Fe = {best_score + log_fefe:.3f} A "
                        f"(Gauss reference {log_fefe:.3f} A).")
            # Fuer pdb_cl: the chosene Paar als pdb_fe[0,1] umsortieren
            pdb_fe = [pdb_fe[i], pdb_fe[j]]
    else:
        fe1p = np.array(pdb_fe[0]); fe2p = np.array(pdb_fe[1])

    # PDB-S-Atome of the cluster
    pdb_s_all = [(a["x"], a["y"], a["z"]) for a in pdb_data["all_h"]
                 if a["element"] == "S"]
    pdb_s_cl  = [sv for sv in pdb_s_all
                 if (np.linalg.norm(np.array(sv) - fe1p) < 3.0 and
                     np.linalg.norm(np.array(sv) - fe2p) < 3.0)]
    if len(pdb_s_cl) < 2:
        pdb_s_cl = sorted(
            pdb_s_all,
            key=lambda s: min(np.linalg.norm(np.array(s) - fe1p),
                              np.linalg.norm(np.array(s) - fe2p)))[:2]
    if len(pdb_s_cl) < 2:
        runlog.warn("Kabsch: Less than 2 cluster S atoms in PDB.")
        return None, None, float("nan"), False

    try:
        pdb_cl = np.array(
            [pdb_fe[0], pdb_fe[1], pdb_s_cl[0], pdb_s_cl[1]], dtype=float)

        rough_t = log_cl_arr.mean(0) - pdb_cl.mean(0)
        pdb_sh  = pdb_cl + rough_t
        order, used = [], set()
        for pi in range(4):
            best = min(
                ((np.linalg.norm(pdb_sh[pi] - log_cl_arr[li]), li)
                 for li in range(4) if li not in used))
            used.add(best[1]); order.append(best[1])
        pdb_m = pdb_cl[order]

        P    = pdb_m      - pdb_m.mean(0)
        Q    = log_cl_arr - log_cl_arr.mean(0)
        U, _, Vt = np.linalg.svd(P.T @ Q)
        ds   = np.linalg.det(Vt.T @ U.T)
        R    = Vt.T @ np.diag([1, 1, ds]) @ U.T
        t    = log_cl_arr.mean(0) - R @ pdb_m.mean(0)

        P_t  = (R @ pdb_m.T).T + t
        rmsd = float(np.sqrt(np.mean(np.sum((P_t - log_cl_arr)**2, axis=1))))
        return R, t, rmsd, True

    except Exception as exc:
        runlog.warn(f"Kabsch alignment failed: {exc}")
        return None, None, float("nan"), False


# ===========================================================================
# Koordinations-Erkennung
# ===========================================================================

@dataclass
class LigandInfo:
    """Informationen over einen Fe-ligands-Kontakt.

    Attributes
    ----------
    fe_idx : int
        Index in ``fe_c`` (0 = Fe1, 1 = Fe2).
    fe_center : int
        Gaussian-Center-Nummer of the Fe atoms.
    lig_center : int
        Gaussian-Center-Nummer of the koordinierenden Atoms.
    lig_element : str
        element of the ligands-Atoms (``"N"``, ``"S"``, ``"O"``).
    lig_aname : str
        PDB atomname (e.g. ``"ND1"``).
    res_num : int
        Residuennummer in the PDB.
    res_name : str
        Residuenname (e.g. ``"HIS"``).
    res_label : str
        Anzeigebezeichnung (e.g. ``"His 255"``).
    bond_vec : ndarray of shape (3,)
        Einheitsvektor Fe → Ligand in Gaussian-Koordinaten.
    bond_len : float
        equilibrium bond length in Angstrom.
    h_center : int or None
        Gaussian-Center of the H-Atoms (nur protoniertes His).
    hn_vec : ndarray of shape (3,) or None
        Einheitsvektor N → H.
    hn_len : float or None
        N-H bond length in Angstrom.
    his_protonated : bool
        ``True`` if His protoniert is and H-Atom in the Gaussian enthalten.
    his_hn_center : int or None
        Gaussian-Center of the N-Atoms the H-bond contributes.
        Kann von ``lig_center`` (Fe-koordinierendes N) abweichen:
        e.g. NE2 if the His through ND1 koordiniert (Rieske-Typ).
    coord_n_aname : str or None
        Koordinierendes N-Atom (``"ND1"`` or ``"NE2"`` for His).
    """

    fe_idx:        int
    fe_center:     int
    lig_center:    int
    lig_element:   str
    lig_aname:     str
    res_num:       int
    res_name:      str
    res_label:     str
    bond_vec:      np.ndarray
    bond_len:      float
    h_center:      Optional[int]        = None
    hn_vec:        Optional[np.ndarray] = None
    hn_len:        Optional[float]      = None
    his_protonated: bool                = False
    his_hn_center: Optional[int]        = None   # N with H (if applicable NE2 ≠ ND1)
    coord_n_aname: Optional[str]        = None


@dataclass
class CoordInfo:
    """Vollstaendige Koordinations-Information of the [2Fe-2S]-Clusters.

    Attributes
    ----------
    ligands : list of LigandInfo
        Alle Fe-ligands-Kontakte (automatically erkannt).
    group_map : dict of {str: list of int}
        Mapping Residuen-Label → Gaussian-Center-Nummern (nur heavy atoms).
        Wird for OOP/INP-Gruppen-Amplituden used.
    pdb_to_center : dict of {int: int}
        Mapping PDB atom listn-Index → Gaussian-Center-Nummer.
        Basis of the korrekten SS-Index-Mappings (Bugfix B2).
    his_ligand_labels : list of str
        Labels aller His-ligands (fuer bedingtes His_HN-Sheet).
    pcet_info : object or None
        Vorcomputede PCET-Atom-Info (``pcet_et.PcetAtomInfo``).
        Wird in the main script einmalig after the ligands-Erkennung gefuellt.
    et_info : object or None
        Vorcomputede ET-Atom-Info (``pcet_et.EtAtomInfo``).
    atoms_h : list of dict or None
        atom list *mit* hydrogen (fuer PCET-Geometrie).
    idx_map_h : dict of {int: int} or None
        Center -> Index in ``atoms_h``.
    """

    ligands:          List[LigandInfo]
    group_map:        Dict[str, List[int]] = field(default_factory=dict)
    pdb_to_center:    Dict[int, int]       = field(default_factory=dict)
    his_ligand_labels: List[str]           = field(default_factory=list)
    # PCET/ET (Hardening v3.0 #8) — optional, is in the main script gefuellt
    pcet_info:        Optional[object]     = None
    et_info:          Optional[object]     = None
    atoms_h:          Optional[List[Dict]] = None
    idx_map_h:        Optional[Dict[int, int]] = None


def _build_pdb_to_center_map(
        pdb_data: Dict,
        atoms: List[Dict],
        idx_map: Dict[int, int],
        R: np.ndarray,
        t: np.ndarray,
        cfg: Config,
) -> Tuple[Dict[int, int], List[float], int]:
    """Creates the Mapping PDB atom listn-Index → Gaussian-Center-Nummer.

    Parameters
    ----------
    pdb_data : dict
        Ausgabe von ``parse_pdb``.
    atoms : list of dict
        Gaussian-atom list (nur heavy atoms).
    idx_map : dict of {int: int}
        Center → Atom-Index.
    R : ndarray of shape (3, 3)
        Kabsch-Rotationsmatrix (PDB → Gaussian).
    t : ndarray of shape (3,)
        Kabsch-Translationsvektor.
    cfg : Config
        Benecessaryt ``coord_match_tol``.

    Returns
    -------
    pdb_to_center : dict of {int: int}
        Mapping PDB atom-Index (0-basiert, only heavy atoms) →
        Gaussian-Center-Nummer.
    match_dists : list of float
        Euklidische Abstände [Å] aller erfolgreichen Zuordnungen;
        leer if no Zuordnung möglich.
    n_ambiguous : int
        Number of PDB atome, für the mehr als ein Gaussian-Kandidat innerhalb
        the tolerance lag (nur erstes Atom is zugeordnet).

    Notes
    -----
    Bugfix B2: Stellt the korrekte Basis for the SS-Index-Mapping bereit.
    Die drei Rückgabewerte are in ``find_coordinating_residues`` direkt
    entpackt: ``pdb_to_center, match_dists, n_ambiguous = _build_pdb_to_center_map(...)``.
    """
    log_by_elem: Dict[str, List[Tuple]] = defaultdict(list)
    for a in atoms:
        log_by_elem[a["symbol"].upper()].append(
            (a["center"], a["x"], a["y"], a["z"]))

    pdb_to_center: Dict[int, int] = {}
    used: Set[int] = set()
    # tolerance: 2x coord_match_tol. Robustheit before Strenge.
    tol2 = (cfg.coord_match_tol * 2) ** 2

    pdb_heavy = [a for a in pdb_data["atoms_h"] if not a["is_h"]]
    match_dists: List[float] = []
    n_ambiguous = 0

    for pi, pa in enumerate(pdb_heavy):
        gxyz  = R @ np.array([pa["x"], pa["y"], pa["z"]]) + t
        elem  = pa["element"]
        best_d2 = tol2
        best_c  = None
        n_cand  = 0
        for ctr, lx, ly, lz in log_by_elem.get(elem, []):
            if ctr in used:
                continue
            d2 = float(np.sum((gxyz - np.array([lx, ly, lz]))**2))
            if d2 < tol2:
                n_cand += 1
            if d2 < best_d2:
                best_d2 = d2; best_c = ctr
        if best_c is not None:
            pdb_to_center[pi] = best_c
            used.add(best_c)
            match_dists.append(float(np.sqrt(best_d2)))
            if n_cand > 1:
                n_ambiguous += 1

    return pdb_to_center, match_dists, n_ambiguous


def find_coordinating_residues(
        pdb_data:     Dict,
        atoms_all:    List[Dict],
        idx_map_all:  Dict[int, int],
        atoms_heavy:  List[Dict],
        idx_map_heavy: Dict[int, int],
        fe_c:         List[int],
        s_c:          List[int],
        R:            np.ndarray,
        t:            np.ndarray,
        cfg:          Config,
        runlog:       RunLog,
) -> CoordInfo:
    """Erkennt automatically alle koordinierenden amino acidn.

    Verwendet Kabsch-transformierte PDB-Koordinaten and the in
    ``cfg.fe_coord_cutoffs`` defineden Abstands-Schwellwerte.
    H-Atome are only for the His-Protonierungserkennung used
    and fliessen not in the Feature-Matrix ein.

    Parameters
    ----------
    pdb_data : dict
        Ausgabe von ``parse_pdb``.
    atoms_all : list of dict
        Gaussian-atom list MIT hydrogen.
    idx_map_all : dict of {int: int}
        Center → Index in ``atoms_all``.
    atoms_heavy : list of dict
        Gaussian-atom list OHNE hydrogen.
    idx_map_heavy : dict of {int: int}
        Center → Index in ``atoms_heavy``.
    fe_c : list of int
        Center-Nummern the Fe atoms.
    s_c : list of int
        Center-Nummern the S-Atome.
    R : ndarray of shape (3, 3)
        Kabsch-Rotationsmatrix.
    t : ndarray of shape (3,)
        Kabsch-Translationsvektor.
    cfg : Config
        Benecessaryt ``fe_coord_cutoffs``, ``coord_match_tol``.
    runlog : RunLog
        Fuer Warnmeldungen and Statistiken.

    Returns
    -------
    CoordInfo
        Vollstaendige Koordinations-Information with ``ligands``,
        ``group_map``, ``pdb_to_center`` and ``his_ligand_labels``.
    """
    cluster_centers: Set[int] = set(fe_c) | set(s_c)

    fe_gaus: List[Optional[np.ndarray]] = []
    for fc in fe_c:
        if fc in idx_map_all:
            a = atoms_all[idx_map_all[fc]]
            fe_gaus.append(np.array([a["x"], a["y"], a["z"]]))
        else:
            fe_gaus.append(None)

    # Spatial Index for Gaussian-Atome with H after element
    log_by_elem: Dict[str, List[Tuple]] = defaultdict(list)
    for a in atoms_all:
        log_by_elem[a["symbol"].upper()].append(
            (a["center"], a["x"], a["y"], a["z"]))

    # PDB atom-Index → Gaussian-Center (mit H, for ligands-Matching)
    pdb_to_gaus_h: Dict[int, int] = {}
    used_h: Set[int] = set()
    for pi, pa in enumerate(pdb_data["all_h"]):
        gxyz = R @ np.array([pa["x"], pa["y"], pa["z"]]) + t
        elem = pa["element"] if pa["element"] else "?"
        best_d = cfg.coord_match_tol * 3.0
        best_c = None
        for ctr, lx, ly, lz in log_by_elem.get(elem, []):
            if ctr in used_h:
                continue
            d = float(np.linalg.norm(gxyz - np.array([lx, ly, lz])))
            if d < best_d:
                best_d = d; best_c = ctr
        if best_c is not None:
            pdb_to_gaus_h[pi] = best_c
            used_h.add(best_c)

    # Koordinierende Atome suchen
    ligands:  List[LigandInfo]         = []
    found_res: Set[Tuple]              = set()

    for fe_idx, fc in enumerate(fe_c):
        if fe_gaus[fe_idx] is None:
            continue
        fe_xyz_gaus = fe_gaus[fe_idx]
        fe_xyz_pdb  = R.T @ (fe_xyz_gaus - t)

        for pi, pa in enumerate(pdb_data["all_h"]):
            if pa["is_h"]:
                continue
            elem    = pa["element"]
            cutoff  = cfg.fe_coord_cutoffs.get(elem, 0.0)
            if cutoff == 0.0:
                continue

            pdb_xyz  = np.array([pa["x"], pa["y"], pa["z"]])
            dist_pdb = float(np.linalg.norm(pdb_xyz - fe_xyz_pdb))
            if dist_pdb >= cutoff:
                continue

            lig_center = pdb_to_gaus_h.get(pi)
            if lig_center is None or lig_center in cluster_centers:
                continue

            rnum  = pa["rnum"]
            rname = pa["rname"]
            aname = pa["aname"]
            key   = (fe_idx, rnum, rname)

            if key in found_res:
                for lig in ligands:
                    if (lig.fe_idx == fe_idx and
                            lig.res_num == rnum and
                            lig.res_name == rname and
                            dist_pdb < lig.bond_len):
                        gaus_lig = R @ pdb_xyz + t
                        bvec     = gaus_lig - fe_xyz_gaus
                        bl       = float(np.linalg.norm(bvec))
                        if bl > 1e-10:
                            lig.lig_center  = lig_center
                            lig.lig_element = elem
                            lig.lig_aname   = aname
                            lig.bond_vec    = bvec / bl
                            lig.bond_len    = bl
                continue

            found_res.add(key)
            gaus_lig = R @ pdb_xyz + t
            bvec     = gaus_lig - fe_xyz_gaus
            bl       = float(np.linalg.norm(bvec))
            if bl < 1e-10:
                continue

            lig = LigandInfo(
                fe_idx      = fe_idx,
                fe_center   = fc,
                lig_center  = lig_center,
                lig_element = elem,
                lig_aname   = aname,
                res_num     = rnum,
                res_name    = rname,
                res_label   = f"{rname.capitalize()} {rnum}",
                bond_vec    = bvec / bl,
                bond_len    = bl,
            )
            if rname.upper() in _HIS and elem == "N":
                lig.coord_n_aname = aname
            ligands.append(lig)

    _add_his_hn_info(ligands, pdb_data["all_h"], pdb_to_gaus_h,
                     atoms_all, idx_map_all, R, t, runlog)

    for lig in sorted(ligands, key=lambda l: (l.fe_idx, l.res_num)):
        prot = " [prot. His, H-N verfuegbar]" if lig.his_protonated else ""
        _lig_msg = (f"Fe{lig.fe_idx+1} ← {lig.res_label} "
                    f"({lig.lig_element}, {lig.lig_aname}, "
                    f"d={lig.bond_len:.3f} A){prot}")
        runlog.info(_lig_msg)
        print(f"    {_lig_msg}")

    if not ligands:
        runlog.warn(
            "No coordinating amino acids detected. "
            "PDB-Kette, Cutoffs and Kabsch-RMSD pruefen.")

    group_map     = _build_group_map(ligands, pdb_data, pdb_to_gaus_h, runlog)
    pdb_to_center, match_dists, n_ambiguous = _build_pdb_to_center_map(
        pdb_data, atoms_heavy, idx_map_heavy, R, t, cfg)
    his_labels    = [l.res_label for l in ligands
                     if l.res_name.upper() in _HIS and l.lig_element == "N"]

    # Matching-Qualitaet in the RunLog festhalten (Hardening #4)
    n_pdb_heavy = sum(1 for a in pdb_data["atoms_h"] if not a["is_h"])
    n_matched   = len(pdb_to_center)
    match_pct   = 100.0 * n_matched / n_pdb_heavy if n_pdb_heavy else 0.0
    mean_d = float(np.mean(match_dists)) if match_dists else 0.0
    max_d  = float(np.max(match_dists))  if match_dists else 0.0
    # Matching-Statistik strukturiert speichern (Hardening #4)
    runlog.match_stats = {
        "n_matched":   n_matched,
        "n_total":     n_pdb_heavy,
        "mean_d":      mean_d,
        "max_d":       max_d,
        "n_ambiguous": n_ambiguous,
    }
    runlog.info(
        f"PDB-Gaussian matching: {n_matched}/{n_pdb_heavy} atoms "
        f"({match_pct:.1f}%), mean distance {mean_d:.4f} A, "
        f"max {max_d:.4f} A, ambiguous {n_ambiguous}")
    if match_pct < 80.0:
        runlog.warn(
            f"PDB matching: only {match_pct:.1f}% of heavy PDB atoms "
            f"zugeordnet. SS- and Gruppen-Amplituden koennen uncomplete sein.")

    runlog.group_match["ligands"] = {
        l.res_label: {
            "fe_idx":   l.fe_idx,
            "element":  l.lig_element,
            "res_name": l.res_name,
            "bond_len": l.bond_len,
        } for l in ligands}

    return CoordInfo(
        ligands          = ligands,
        group_map        = group_map,
        pdb_to_center    = pdb_to_center,
        his_ligand_labels = list(set(his_labels)),
    )


def _add_his_hn_info(
        ligands:       List[LigandInfo],
        all_pdb_h:     List[Dict],
        pdb_to_gaus_h: Dict[int, int],
        atoms_all:     List[Dict],
        idx_map_all:   Dict[int, int],
        R: np.ndarray,
        t: np.ndarray,
        runlog=None,
) -> None:
    """Ergaenzt His-ligands with H-N-bondinformation (in-place).

    Finds for each N-koordinierten His-ligands the gebundene H-Atom
    and contributes ``h_center``, ``his_hn_center``, ``hn_vec``, ``hn_len``
    sowie ``his_protonated = True`` in the ``LigandInfo``-Objekt ein.

    Strategie (beide Schritte are versucht):

    1. PDB-basiert: Abgleich over ``pdb_to_gaus_h``-Mapping.
    2. Gaussian-basiert (Fallback): Finds direkt in ``atoms_all`` nach
       H-Atomen nahe the koordinierenden N (ND1) **und** the anderen
       Ring-N (NE2, ~2.3 Å entfernt). Deckt the Rieske-Fall ab, wo
       ND1 Fe koordiniert and H on NE2 sitzt.
    """
    HN_CUT   = 1.20  # N-H bond length (typical ~1.01 Å)
    RING_CUT = 3.00  # ND1-NE2 Abstand in the Imidazolring (PDB-value: ~2.21 Å)

    # Inverse Karte: list-index → Gaussian-Center-Nummer
    # Gaussian-Atom-Dicts: "atomic_num"(int), "symbol"(str), "x/y/z"(float)
    idx_to_ctr = {v: k for k, v in idx_map_all.items()}

    # Alle Gaussian-H-Atome einmal aufbauen (atomic_num==1, KEIN "is_h"!)
    gaus_h = [
        (ai, np.array([a["x"], a["y"], a["z"]]), idx_to_ctr[ai])
        for ai, a in enumerate(atoms_all)
        if a.get("atomic_num") == 1 and ai in idx_to_ctr
    ]

    # PDB-H: welche are per pdb_to_gaus_h gemappt?
    pdb_h_ctrs: set = {
        pdb_to_gaus_h[pi]
        for pi, a in enumerate(all_pdb_h)
        if a.get("is_h") and pi in pdb_to_gaus_h
    }

    for lig in ligands:
        if lig.res_name.upper() not in _HIS or lig.lig_element != "N":
            continue
        if lig.lig_center not in idx_map_all:
            continue

        nd1_idx = idx_map_all[lig.lig_center]
        nd1_pos = np.array([atoms_all[nd1_idx]["x"],
                            atoms_all[nd1_idx]["y"],
                            atoms_all[nd1_idx]["z"]])

        # Alle Ring-N of the Imidazols: ND1 + alle N within RING_CUT
        # (ND1-NE2 in the echten His: ~2.21 Å < RING_CUT=3.0)
        ring_ns: List[Tuple[int, np.ndarray, int]] = [
            (nd1_idx, nd1_pos, lig.lig_center)
        ]
        for ai, a in enumerate(atoms_all):
            if ai == nd1_idx or a.get("atomic_num") != 7:
                continue
            a_pos = np.array([a["x"], a["y"], a["z"]])
            if float(np.linalg.norm(a_pos - nd1_pos)) < RING_CUT:
                ring_ns.append((ai, a_pos, idx_to_ctr[ai]))

        # Fuer jeden Ring-N: suche gebundenes H-Atom in Gaussian-Koordinaten
        # path A: only PDB-identifizierte H (Prioritaet, pruefe zuerst)
        # path B: alle Gaussian-H (Fallback, if PDB no H hat)
        for use_pdb_only in (True, False):
            if lig.his_protonated:
                break
            for n_ai, n_pos, n_ctr in ring_ns:
                for h_ai, h_pos, h_ctr in gaus_h:
                    if use_pdb_only and h_ctr not in pdb_h_ctrs:
                        continue
                    dist = float(np.linalg.norm(h_pos - n_pos))
                    if dist < HN_CUT:
                        hn = h_pos - n_pos
                        lig.his_protonated = True
                        lig.h_center       = h_ctr
                        lig.his_hn_center  = n_ctr
                        lig.hn_vec         = hn / max(dist, 1e-10)
                        lig.hn_len         = dist
                        via = "PDB" if use_pdb_only else "Gaussian-Fallback"
                        _hn_msg = (f"His_HN {via}: {lig.res_label}: "
                                   f"N({n_ctr})-H({h_ctr}), d={dist:.3f} A")
                        if runlog is not None: runlog.info(_hn_msg)
                        print(f"    {_hn_msg}")
                        break
                if lig.his_protonated:
                    break

            if not lig.his_protonated:
                _hn_miss = (f"His_HN NICHT erkannt: {lig.res_label} "
                            f"(Ring-N found: {len(ring_ns)}, "
                            f"Gaussian-H gesamt: {len(gaus_h)}, "
                            f"PDB-H gemappt: {len(pdb_h_ctrs)})")
                if runlog is not None: runlog.warn(_hn_miss)
                print(f"    {_hn_miss}")


def _build_group_map(
        ligands:       List[LigandInfo],
        pdb_data:      Dict,
        pdb_to_gaus_h: Dict[int, int],
        runlog=None,
) -> Dict[str, List[int]]:
    """Creates group_map: Residuen-Label → Gaussian-Center list (nur heavy atoms).

    Wird for the Gruppen-OOP-Analyse (Cys, His, Backbone) used.
    Jedes Residuum taucht only einmal auf, also if es mehrere ligands hat.
    """
    group_map: Dict[str, List[int]] = {}
    pdb_heavy = [a for a in pdb_data["atoms_h"] if not a["is_h"]]
    done: Set[str] = set()

    for lig in ligands:
        lb = lig.res_label
        if lb in done:
            continue
        done.add(lb)
        centers = [
            pdb_to_gaus_h[pi]
            for pi, a in enumerate(pdb_heavy)
            if a["rnum"] == lig.res_num and pi in pdb_to_gaus_h
        ]
        if centers:
            group_map[lb] = centers
            if runlog is not None: runlog.info(f"Group '{lb}': {len(centers)} atoms assigned")
            print(f"    {lb}: {len(centers)} atoms assigned")
        else:
            import warnings as _w
            _w.warn(
                f"Group '{lb}': no Gaussian-Center zugeordnet "
                f"(Kabsch-Matching uncomplete?)", UserWarning)

    return group_map


# ===========================================================================
# SS-Index-Mapping  (Bugfix B2)
# ===========================================================================

def build_ss_center_map(
        ss_elements:   List[Dict],
        pdb_data:      Dict,
        pdb_to_center: Dict[int, int],
        chain_filter:  str = "",
) -> Dict[str, List[int]]:
    """Creates the Mapping SS-element-Name → Gaussian-Center-Nummern.

    Parameters
    ----------
    ss_elements : list of dict
        HELIX/SHEET-Records from ``parse_pdb``.
    pdb_data : dict
        Ausgabe von ``parse_pdb``.
    pdb_to_center : dict of {int: int}
        PDB atom listn-Index → Gaussian-Center (aus ``find_coordinating_residues``).
    chain_filter : str, optional
        Ketten-ID; leer = alle Ketten.

    Returns
    -------
    dict of {str: list of int}
        Mapping SS-element-Name → Gaussian-Center-Nummern.

    Notes
    -----
    Bugfix B2: Verwendet ``pdb_to_center`` statt einfacher PDB listnindizes.
    Damit are SS-elements korrekt on the Gaussian-Eigenvektor-Zeilen
    abgebildet.
    """
    pdb_heavy     = [a for a in pdb_data["atoms_h"] if not a["is_h"]]
    ss_center_map: Dict[str, List[int]] = {}

    for elem in ss_elements:
        if chain_filter and elem.get("chain", "") != chain_filter:
            continue
        centers = [
            pdb_to_center[pi]
            for pi, a in enumerate(pdb_heavy)
            if (not chain_filter or a["chain"] == chain_filter) and
               elem["res_start"] <= a["rnum"] <= elem["res_end"] and
               pi in pdb_to_center
        ]
        if centers:
            ss_center_map[elem["name"]] = centers

    return ss_center_map


def get_calpha_centers(
        pdb_data:      Dict,
        atoms:         List[Dict],
        pdb_to_center: Dict[int, int],
) -> List[Tuple[int, int]]:
    """Returns Gaussian-Center-Nummern aller Cα-Atome .

    Parameters
    ----------
    pdb_data : dict
        Ausgabe von ``parse_pdb``.
    atoms : list of dict
        Gaussian-atom list (nur heavy atoms).
    pdb_to_center : dict of {int: int}
        PDB atom listn-Index → Gaussian-Center.

    Returns
    -------
    list of (int, int)
        List of ``(gaussian_center, residue_num)`` for each Cα.
    """
    pdb_heavy = [a for a in pdb_data["atoms_h"] if not a["is_h"]]
    return [
        (pdb_to_center[pi], a["rnum"])
        for pi, a in enumerate(pdb_heavy)
        if a["aname"] == "CA" and pi in pdb_to_center
    ]


def detect_ss_dssp(
        pdb_data: Dict,
        chain_filter: str = "A",
        min_helix: int = 4,
        min_sheet: int = 4,
) -> List[Dict]:
    """Erkennt Sekundaerstrukturelemente after the DSSP-Algorithmus
    (Kabsch & Sander 1983, [11]).

    Computes for each Residuenpaar the H-bond energy zwischen
    Amid-NH (donor) and Carbonyl-C=O (acceptor) after the DSSP-Formel:

        E = 0.084 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) * 332  [kcal/mol]

    A H-Bruecke liegt before if E < -0.5 kcal/mol. Die Klassifikation
    folgt the DSSP-Mustern:

    - Alpha-Helix (H): H-Bruecke between i and i+4 (sowie 3_10: i+3,
      Pi: i+5)
    - Beta-Strang (E): H-Bruecken to nicht-benachbarten Residuen (|i-j| > 5)

    Fehlende Amide-H are from the N-C(i-1)-Richtung geschaetzt
    (DSSP-Naeherung, bond length 1.01 Angstrom).

    Wird als primaerer Fallback called if the PDB-file keine
    HELIX/SHEET-Records enthaelt. For leerer Ausgabe greift
    ``detect_ss_phipsi`` als sekundaerer Fallback.

    Parameters
    ----------
    pdb_data : dict
        Ausgabe von ``parse_pdb``. Benecessaryt ``"atoms_h"`` (heavy atoms
        and H-Atome).
    chain_filter : str
        Ketten-ID; leer = alle Ketten. Standard: ``"A"``.
    min_helix : int
        Minimale length (Residuen) for a Helix. Standard: 4.
    min_sheet : int
        Minimale length (Residuen) for einen Beta-Strang. Standard: 4.

    Returns
    -------
    list of dict
        SS-elements with ``"auto_detected": True``.

    Notes
    -----
    Referenz: Kabsch, W. & Sander, C. (1983). Dictionary of protein
    secondary structure. *Biopolymers*, **22**(12), 2577-2637.
    doi:10.1002/bip.360221211 [11]

    Laufzeit: O(N^2) in the Number of Residuen; typical < 0.5 s fuer
    Proteine with ~300 Residuen.
    """

    # 1. Backbone-Atome per Residue sammeln (atoms_h enthaelt H-Atome)
    residues: Dict[int, Dict[str, np.ndarray]] = {}
    for a in pdb_data.get("atoms_h", []):
        if chain_filter and a.get("chain", "") != chain_filter:
            continue
        aname = a["aname"].strip()
        rnum  = a["rnum"]
        if rnum not in residues:
            residues[rnum] = {}
        coord = np.array([a["x"], a["y"], a["z"]], dtype=float)
        if   aname == "N":                                    residues[rnum]["N"]  = coord
        elif aname == "CA":                                   residues[rnum]["CA"] = coord
        elif aname == "C":                                    residues[rnum]["C"]  = coord
        elif aname == "O":                                    residues[rnum]["O"]  = coord
        elif aname in ("H", "HN", "1H", "H1") and "H" not in residues[rnum]:
            residues[rnum]["H"] = coord   # Amide-H (erstes foundes)

    rnums = sorted(residues.keys())
    if not rnums:
        return []

    # 2. Fehlende Amide-H-Positionen schaetzen (DSSP-Naeherung)
    for i, rnum in enumerate(rnums):
        res = residues[rnum]
        if "H" not in res and "N" in res and i > 0:
            prev = residues.get(rnums[i - 1], {})
            if "C" in prev:
                direction = res["N"] - prev["C"]
                norm = float(np.linalg.norm(direction))
                if norm > 1e-6:
                    res["H"] = res["N"] + direction / norm * 1.01

    # 3. H-bond energyn berechnen
    # Nur Paare with r_OH <= 3.5 A are computed (Performance-Cutoff)
    hbonds: set = set()   # (rnum_j, rnum_i): j spendet H an O von i

    for j_idx, rnum_j in enumerate(rnums):
        res_j = residues[rnum_j]
        if "N" not in res_j or "H" not in res_j:
            continue
        N_j = res_j["N"]
        H_j = res_j["H"]

        for i_idx, rnum_i in enumerate(rnums):
            if abs(i_idx - j_idx) < 2:
                continue                # Benachbarte Residuen skip
            res_i = residues[rnum_i]
            if "O" not in res_i or "C" not in res_i:
                continue

            O_i = res_i["O"]
            C_i = res_i["C"]

            r_OH = float(np.linalg.norm(O_i - H_j))
            if r_OH > 3.5:
                continue               # Ausserhalb Cutoff

            r_ON = float(np.linalg.norm(O_i - N_j))
            r_CH = float(np.linalg.norm(C_i - H_j))
            r_CN = float(np.linalg.norm(C_i - N_j))

            if min(r_OH, r_ON, r_CH, r_CN) < 0.5:
                continue               # Sterische Kollision

            E = 0.084 * (1.0/r_ON + 1.0/r_CH - 1.0/r_OH - 1.0/r_CN) * 332.0
            if E < -0.5:
                hbonds.add((rnum_j, rnum_i))

    # 4. SS-Klassifikation after DSSP-Mustern
    ss_type: Dict[int, str] = {rnum: "C" for rnum in rnums}
    rnum_to_idx: Dict[int, int] = {rnum: i for i, rnum in enumerate(rnums)}

    # Alpha-Helix (i→i+4), 3_10-Helix (i→i+3), Pi-Helix (i→i+5)
    for j_idx, rnum_j in enumerate(rnums):
        for offset in (3, 4, 5):
            if j_idx >= offset:
                rnum_i = rnums[j_idx - offset]
                if (rnum_j, rnum_i) in hbonds:
                    for k in range(max(0, j_idx - offset + 1), j_idx + 1):
                        ss_type[rnums[k]] = "H"

    # Beta-Strang: H-Bruecken to entfernten Residuen (|i-j| > 5)
    for rnum_j, rnum_i in hbonds:
        if abs(rnum_to_idx[rnum_j] - rnum_to_idx[rnum_i]) > 5:
            if ss_type.get(rnum_j, "C") == "C":
                ss_type[rnum_j] = "E"
            if ss_type.get(rnum_i, "C") == "C":
                ss_type[rnum_i] = "E"

    # 5. Konsekutive gleichartige Residuen to SS-elementsn gruppieren
    ss_elements: List[Dict] = []
    i = 0
    while i < len(rnums):
        t = ss_type.get(rnums[i], "C")
        if t == "C":
            i += 1
            continue
        j = i + 1
        while (j < len(rnums) and
               ss_type.get(rnums[j], "C") == t and
               rnums[j] - rnums[j - 1] <= 2):
            j += 1
        r1  = rnums[i]
        r2  = rnums[j - 1]
        ch  = chain_filter or "A"
        if t == "H" and (r2 - r1 + 1) >= min_helix:
            ss_elements.append({
                "type":          "helix",
                "chain":         ch,
                "res_start":     r1,
                "res_end":       r2,
                "name":          f"Helix_{ch}_{r1}_{r2}",
                "auto_detected": True,
            })
        elif t == "E" and (r2 - r1 + 1) >= min_sheet:
            ss_elements.append({
                "type":          "sheet",
                "chain":         ch,
                "res_start":     r1,
                "res_end":       r2,
                "name":          f"Sheet_S1_{ch}_{r1}",
                "auto_detected": True,
            })
        i = j

    return ss_elements


def detect_ss_phipsi(
        pdb_data: Dict,
        chain_filter: str = "A",
        min_helix: int = 4,
        min_sheet: int = 4,
) -> List[Dict]:
    """Erkennt Sekundaerstrukturelemente from Backbone-Dihedralwinkeln phi/psi.

    Sekundaerer Fallback, the only called is if ``detect_ss_dssp``
    no elements findet (e.g. because H-Atome complete fehlen).

    Verwendet vereinfachte Ramachandran-Regionen (Kabsch & Sander 1983, [11]):

    - Alpha-Helix: phi in [-90, -30] Grad, psi in [-70, -10] Grad
    - Beta-Strang: phi in [-160, -60] Grad, psi in [90, 180] oder
      [-180, -120] Grad

    Parameters
    ----------
    pdb_data : dict
        Ausgabe von ``parse_pdb``.
    chain_filter : str
        Ketten-ID; leer = alle Ketten. Standard: ``"A"``.
    min_helix : int
        Minimale length (Residuen). Standard: 4.
    min_sheet : int
        Minimale length (Residuen). Standard: 4.

    Returns
    -------
    list of dict
        SS-elements with ``"auto_detected": True``.
    """

    def _dihedral(p0: np.ndarray, p1: np.ndarray,
                  p2: np.ndarray, p3: np.ndarray) -> float:
        b0 = p0 - p1
        b1 = p2 - p1
        b2 = p3 - p2
        n  = float(np.linalg.norm(b1))
        if n < 1e-10:
            return 0.0
        b1n = b1 / n
        v   = b0 - np.dot(b0, b1n) * b1n
        w   = b2 - np.dot(b2, b1n) * b1n
        return float(np.degrees(np.arctan2(
            float(np.dot(np.cross(b1n, v), w)),
            float(np.dot(v, w)))))

    residues: Dict[int, Dict[str, np.ndarray]] = {}
    for a in pdb_data.get("atoms", []):
        if chain_filter and a.get("chain", "") != chain_filter:
            continue
        aname = a["aname"].strip()
        if aname not in ("N", "CA", "C", "O"):
            continue
        rnum = a["rnum"]
        if rnum not in residues:
            residues[rnum] = {}
        residues[rnum][aname] = np.array([a["x"], a["y"], a["z"]], dtype=float)

    rnums = sorted(residues.keys())
    ss_type: Dict[int, str] = {}
    _n_fail_phipsi = 0  # Residuen without berechenbare Dihedralwinkel

    for i, rnum in enumerate(rnums):
        res = residues[rnum]
        if not all(k in res for k in ("N", "CA", "C")):
            ss_type[rnum] = "C"
            continue

        phi: Optional[float] = None
        psi: Optional[float] = None

        if i > 0:
            prev = residues.get(rnums[i - 1], {})
            if "C" in prev:
                try:
                    phi = _dihedral(prev["C"], res["N"], res["CA"], res["C"])
                except Exception:
                    _n_fail_phipsi += 1  # phi undefined (terminal/Pro/bad coords)

        if i < len(rnums) - 1:
            nxt = residues.get(rnums[i + 1], {})
            if "N" in nxt:
                try:
                    psi = _dihedral(res["N"], res["CA"], res["C"], nxt["N"])
                except Exception:
                    _n_fail_phipsi += 1  # psi undefined (terminal/bad coords)

        if phi is None or psi is None:
            ss_type[rnum] = "C"
            _n_fail_phipsi += 1
        elif -90.0 <= phi <= -30.0 and -70.0 <= psi <= -10.0:
            ss_type[rnum] = "H"
        elif (-160.0 <= phi <= -60.0 and
              (90.0 <= psi <= 180.0 or -180.0 <= psi <= -120.0)):
            ss_type[rnum] = "E"
        else:
            ss_type[rnum] = "C"

    ss_elements: List[Dict] = []
    i = 0
    while i < len(rnums):
        t = ss_type.get(rnums[i], "C")
        if t == "C":
            i += 1
            continue
        j = i + 1
        while (j < len(rnums) and
               ss_type.get(rnums[j], "C") == t and
               rnums[j] - rnums[j - 1] <= 2):
            j += 1
        r1  = rnums[i]
        r2  = rnums[j - 1]
        ch  = chain_filter or "A"
        if t == "H" and (r2 - r1 + 1) >= min_helix:
            ss_elements.append({
                "type": "helix", "chain": ch,
                "res_start": r1, "res_end": r2,
                "name": f"Helix_{ch}_{r1}_{r2}",
                "auto_detected": True,
            })
        elif t == "E" and (r2 - r1 + 1) >= min_sheet:
            ss_elements.append({
                "type": "sheet", "chain": ch,
                "res_start": r1, "res_end": r2,
                "name": f"Sheet_S1_{ch}_{r1}",
                "auto_detected": True,
            })
        i = j
    if _n_fail_phipsi > 0:
        import warnings as _w
        _w.warn(
            f"detect_ss_phipsi: {_n_fail_phipsi} von {len(rnums)} Residuen "
            f"ohne berechenbare Dihedralwinkel (phi or psi None) -- "
            f"als Coil klassifiziert. Geometrie the PDB pruefen.",
            UserWarning, stacklevel=2)

    return ss_elements


__version__ = "1.4"  # modenanalyse v1.4
