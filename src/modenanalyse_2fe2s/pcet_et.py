# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

# -*- coding: utf-8 -*-
"""
pcet_et.py — v3.7 geometry adapter
====================================

Diese file stellt the **Atom-Setup-functions** for the Marcus-Hush-
Reorganisations-Pipeline bereit. Die eigentliche Reorg-Berechnung
occurs in :mod:`reorganization`.

Oeffentliche API
----------------

* :class:`PcetAtomInfo` -- Datenstruktur for protonierte His-ligands
  and ihre H-acceptor-Paare.
* :class:`EtAtomInfo` -- Datenstruktur for Cluster-Atome and Fe-Fe-Achse.
* :func:`build_pcet_info` -- baut PcetAtomInfo from ligands list +
  atom list with H. Is einmal per Lauf called.
* :func:`build_et_info` -- baut EtAtomInfo from the cluster-Geometrie.
* :func:`find_hbond_acceptors_for_h` -- sucht to a H-Atom alle
  acceptor-Atome in the Umkreis (3.5 A Standard-Cutoff).

Wissenschaftlicher Hintergrund
------------------------------

Marcus, J. Chem. Phys. 24, 966 (1956); ders., Rev. Mod. Phys. 65, 599
(1993): Marcus-Theorie of the Elektronentransfers; reorganization energy
als zentrale Aktivierungs-Groesse.

Hammes-Schiffer, Soudackov, J. Phys. Chem. B 112, 14108 (2008):
Theoretischer Rahmen for PCET in Loesung, Proteinen and Elektroden.
PCET als Mehr-Mode-Phaenomen.

Bergner, Dechert, Demeshko, Kupper, Mayer, Meyer, JACS 139, 701 (2017):
mitoNEET-Modell zeigt CPET to TEMPO; experimentelle Bestaetigung
der reorganization energy als zentraler PCET-Effizienz-Faktor.

Saouma, Pinney, Mayer, Inorg. Chem. 53, 3153 (2014):
Erste Charakterisierung von CPET an synthetischen [2Fe-2S]-Modellen.
"""

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ===========================================================================
# Datenstrukturen
# ===========================================================================

@dataclass
class PcetAtomInfo:
    """Vorcomputede Info over PCET-aktive Atome.

    Diese Struktur is einmal per Lauf gebaut (nach
    :func:`find_coordinating_residues`), and dann per Mode used,
    without neuerliche Geometriesuche. Das spart per Mode den
    H-Bond-acceptor-Scan.

    Attributes
    ----------
    cluster_centers : list of int
        Gaussian-Center the vier Cluster-Atome (Fe1, Fe2, S1, S2).
    his_n_centers : list of int
        Center the His-N-Atome, the H tragen (kann von dem
        Fe-koordinierenden N abweichen — z. B. Nε if the Fe an Nδ
        gebunden ist).
    his_h_centers : list of int
        Center the dazugehoerigen H-Atome (eines per N).
    acceptor_centers_per_h : list of list of int
        Pro His-H a List of acceptor-Center in
        ``cfg.pcet_hbond_cutoff_a`` Reichweite of the H-Atoms.
    eq_distances_per_h : list of list of float
        Pro His-H the Gleichgewichtsdistanzen H...acceptor in Angstrom
        (gleiche Reihenfolge wie ``acceptor_centers_per_h``).
    n_his : int
        Number of found His-ligands with protoniertem N-Atom.
    enabled : bool
        ``True`` if ``n_his > 0`` and mindestens ein H...acceptor-Paar
        found wurde. For ``False`` liefert :func:`compute_pcet_score`
        ``NaN``.
    diagnose : str
        Kurzes Diagnose-String for the REPORT (Number of His, acceptoren).
    """

    cluster_centers:        List[int]
    his_n_centers:          List[int]            = field(default_factory=list)
    his_h_centers:          List[int]            = field(default_factory=list)
    acceptor_centers_per_h: List[List[int]]      = field(default_factory=list)
    eq_distances_per_h:     List[List[float]]    = field(default_factory=list)
    n_his:                  int                  = 0
    enabled:                bool                 = False
    diagnose:               str                  = ""


@dataclass
class EtAtomInfo:
    """Vorcomputede Info over ET-aktive Atome.

    Attributes
    ----------
    cluster_centers : list of int
        Gaussian-Center the vier Cluster-Atome.
    ligand_centers : list of int
        Center aller Fe-koordinierenden ligands-Atome (Nδ von His,
        Sγ von Cys, OD von Asp/Glu, ...).
    fe_centers : tuple of (int, int)
        Center the zwei Fe atoms (fuer dr_FeFe-Berechnung).
    """

    cluster_centers: List[int]
    ligand_centers:  List[int]            = field(default_factory=list)
    fe_centers:      Tuple[int, int]      = (0, 0)


# H-Bond-acceptor-elements (Atom-Nummer)
_ACCEPTOR_ELEMENTS: set = {7, 8, 9, 16}     # N, O, F, S


def _atom_pos(atom: Dict) -> np.ndarray:
    """Position of an atom in Angstrom, as ndarray(3,)."""
    return np.array([atom["x"], atom["y"], atom["z"]], dtype=float)


# ===========================================================================
# Pre-Computation (einmal per Lauf)
# ===========================================================================

def find_hbond_acceptors_for_h(
        h_pos: np.ndarray,
        atoms_h: List[Dict],
        h_center: int,
        n_donor_center: int,
        cluster_centers: List[int],
        ligand_centers: List[int],
        cutoff_a: float = 4.0,
) -> List[Tuple[int, float]]:
    r"""Finds H-Bond-acceptor-Atome in the Naehe a H-Atoms.

    acceptoren are Atome with element N, O, F or S within von
    ``cutoff_a`` Angstrom vom H-Atom. Folgende Ausschluesse:

      * H-Atom selbst and donor-N (``h_center``, ``n_donor_center``)
      * Cluster-Atome (Fe1, Fe2, S1, S2)
      * Andere Fe-koordinierende ligands-Atome (sie are not freie
        acceptoren, but rather Teil of the bereits-koordinierten ligands-
        Geflechts; insbesondere for the donor-Imidazol the andere
        N-Atom of the Rings)

    Parameters
    ----------
    h_pos : ndarray of shape (3,)
        Position of the H-Atoms in Angstrom.
    atoms_h : list of dict
        atom list *mit* H-Atomen.
    h_center : int
        Gaussian-Center of the H-Atoms (to Selbstausschluss).
    n_donor_center : int
        Gaussian-Center of the donor-N (to Ausschluss from acceptoren).
    cluster_centers : list of int
        Cluster-Atome (Fe1, Fe2, S1, S2) — are ausgeschlossen.
    ligand_centers : list of int
        Alle Fe-koordinierenden ligands-Atome — are ausgeschlossen.
    cutoff_a : float, optional
        Maximaler H...acceptor-Abstand in Angstrom. default 4.0.

    Returns
    -------
    list of (int, float)
        Sortierte Liste ``[(center, distance), ...]`` of the found
        acceptoren, aufsteigend after Distanz.
    """
    excl: set = (set(cluster_centers) | set(ligand_centers)
                 | {h_center, n_donor_center})
    out: List[Tuple[int, float]] = []
    for atom in atoms_h:
        if atom["center"] in excl:
            continue
        if atom["atomic_num"] not in _ACCEPTOR_ELEMENTS:
            continue
        d = float(np.linalg.norm(_atom_pos(atom) - h_pos))
        if d <= cutoff_a:
            out.append((atom["center"], d))
    out.sort(key=lambda x: x[1])
    return out


def build_pcet_info(
        coord_info,
        atoms_h: List[Dict],
        idx_map_h: Dict[int, int],
        fe_c: List[int],
        s_c: List[int],
        cfg,
) -> PcetAtomInfo:
    r"""Baut the PCET-Atom-Information einmal before the modenschleife.

    Geht through ``coord_info.ligands``, identifiziert protonierte
    histidinee (``LigandInfo.his_protonated == True``), and sucht fuer
    jedes His-N--H the H-Bond-acceptoren in
    ``cfg.pcet_hbond_cutoff_a`` Reichweite.

    For a System without His-ligands (z. B. rein Cys-koordinierte
    Cluster) is ``enabled=False`` gesetzt; PCET is physikalisch
    not possible.

    Parameters
    ----------
    coord_info : CoordInfo
        Aus ``find_coordinating_residues``.
    atoms_h : list of dict
        atom list with H.
    idx_map_h : dict of {int: int}
        Center -> Index in ``atoms_h``.
    fe_c, s_c : list of int
        Cluster-Center.
    cfg : Config
        Config; required ``pcet_hbond_cutoff_a``.

    Returns
    -------
    PcetAtomInfo
    """
    cluster_centers = list(fe_c) + list(s_c)
    info = PcetAtomInfo(cluster_centers=cluster_centers)

    cutoff_a = float(getattr(cfg, "pcet_hbond_cutoff_a", 4.0))

    # Alle Fe-koordinierenden Atome (ligands) als acceptor-Ausschluss
    all_lig_centers = [lig.lig_center for lig in coord_info.ligands]

    n_his_total = 0
    n_pairs_total = 0
    for lig in coord_info.ligands:
        # Nur protonierte His-ligands with identifiziertem H
        if lig.res_name != "HIS":
            continue
        if not lig.his_protonated:
            continue
        if lig.h_center is None or lig.his_hn_center is None:
            continue

        n_his_total += 1
        h_idx = idx_map_h.get(lig.h_center)
        if h_idx is None:
            continue
        h_pos = _atom_pos(atoms_h[h_idx])

        accs = find_hbond_acceptors_for_h(
            h_pos, atoms_h,
            h_center=lig.h_center,
            n_donor_center=lig.his_hn_center,
            cluster_centers=cluster_centers,
            ligand_centers=all_lig_centers,
            cutoff_a=cutoff_a,
        )

        info.his_n_centers.append(lig.his_hn_center)
        info.his_h_centers.append(lig.h_center)
        info.acceptor_centers_per_h.append([c for c, _d in accs])
        info.eq_distances_per_h.append([d for _c, d in accs])
        n_pairs_total += len(accs)

    info.n_his = n_his_total
    info.enabled = (n_his_total > 0 and n_pairs_total > 0)
    info.diagnose = (
        f"{n_his_total} His-ligands, "
        f"{n_pairs_total} H-Bond-Paare (Cutoff {cutoff_a:.1f} A)"
    )
    return info


def build_et_info(
        coord_info,
        fe_c: List[int],
        s_c: List[int],
) -> EtAtomInfo:
    r"""Baut the ET-Atom-Information.

    Cluster + alle Fe-koordinierenden ligands-Atome — His-unabhaengig.

    Parameters
    ----------
    coord_info : CoordInfo
        Aus ``find_coordinating_residues``.
    fe_c, s_c : list of int
        Cluster-Center.

    Returns
    -------
    EtAtomInfo
    """
    cluster_centers = list(fe_c) + list(s_c)
    ligand_centers = [lig.lig_center for lig in coord_info.ligands]
    fe_pair: Tuple[int, int] = (
        fe_c[0] if len(fe_c) > 0 else 0,
        fe_c[1] if len(fe_c) > 1 else 0,
    )
    return EtAtomInfo(
        cluster_centers=cluster_centers,
        ligand_centers=ligand_centers,
        fe_centers=fe_pair,
    )


__version__ = "3.0"
