# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit-Tests fuer pcet_et.find_hbond_acceptors_for_h.

Testet die H-Bond-Akzeptor-Suche, die im PCET-Setup verwendet wird,
um pro His-H die naheliegenden Akzeptor-Atome (N, O, F, S) zu finden.

Wichtige Aspekte:
* Cutoff-Filterung
* Ausschluss von H, donor-N, Cluster-Atomen, anderen Liganden
* Element-Filter (nur N, O, F, S)
* Sortierung nach Distanz

Aufruf::

    python -m pytest tests/test_pcet_et.py -v
    python tests/test_pcet_et.py
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "src"))

from modenanalyse_2fe2s.pcet_et import find_hbond_acceptors_for_h


def _atom(center: int, atomic_num: int, x: float, y: float, z: float) -> dict:
    """Hilfsfunktion fuer Mock-Atom-Dicts."""
    return {
        "center":     center,
        "atomic_num": atomic_num,
        "x":          x,
        "y":          y,
        "z":          z,
    }


# =============================================================================
# Grundfunktionalitaet
# =============================================================================

def test_find_acceptor_within_cutoff():
    """Ein Wasser-O bei 2.8 A Distanz wird gefunden."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1,  1, 0.0, 0.0, 0.0),    # H selbst
        _atom(2,  7, 0.0, 0.0, 1.0),    # donor-N
        _atom(3,  8, 0.0, 0.0, 2.8),    # H-Bond-Akzeptor (Wasser-O)
        _atom(4,  6, 5.0, 0.0, 0.0),    # weit entferntes C
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos=h_pos,
        atoms_h=atoms_h,
        h_center=1,
        n_donor_center=2,
        cluster_centers=[],
        ligand_centers=[],
        cutoff_a=4.0,
    )
    assert len(accs) == 1
    assert accs[0][0] == 3
    assert abs(accs[0][1] - 2.8) < 1e-9


def test_acceptor_outside_cutoff_skipped():
    """Akzeptor jenseits Cutoff wird ignoriert."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1, 0.0, 0.0, 0.0),
        _atom(2, 7, 0.0, 0.0, 1.0),
        _atom(3, 8, 0.0, 0.0, 5.0),    # 5 A jenseits 4-A-Cutoff
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    assert accs == []


def test_acceptor_at_exact_cutoff_included():
    """Akzeptor genau am Cutoff (<=) wird mitgenommen."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1, 0.0, 0.0, 0.0),
        _atom(2, 7, 0.0, 0.0, 1.0),
        _atom(3, 8, 0.0, 0.0, 4.0),    # exakt 4.0 A
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    assert len(accs) == 1


def test_donor_n_excluded():
    """Das donor-N selbst darf nicht als Akzeptor gewertet werden."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1, 0.0, 0.0, 0.0),
        _atom(2, 7, 0.0, 0.0, 1.0),    # donor-N
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    assert accs == []


def test_h_self_excluded():
    """H-Atom (atomic_num=1) ist sowieso kein Akzeptor; Test verifiziert
    konsistent das Verhalten."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1, 0.0, 0.0, 0.0),    # H selbst
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=99,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    assert accs == []


def test_cluster_atoms_excluded():
    """Fe1, Fe2, S1, S2 (Cluster-Atome) werden nicht als Akzeptoren
    gewertet, auch wenn sie geometrisch in Reichweite sind."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1,  0.0, 0.0, 0.0),
        _atom(2, 7,  0.0, 0.0, 1.0),
        _atom(10, 26, 0.0, 0.0, 2.5),  # Fe1 (atomic_num=26)
        _atom(11, 16, 0.0, 0.0, 3.0),  # Cluster-S (S = atomic_num=16)
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[10, 11],
        ligand_centers=[], cutoff_a=4.0,
    )
    assert accs == []


def test_ligand_atoms_excluded():
    """Andere Fe-koordinierende Liganden-Atome (z.B. anderes His-N,
    Cys-S) werden ausgeschlossen, weil sie zum koordinierten Geflecht
    gehoeren."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1,  1, 0.0, 0.0, 0.0),
        _atom(2,  7, 0.0, 0.0, 1.0),
        _atom(20, 7, 0.0, 0.0, 3.0),    # anderes His-N
        _atom(21, 16, 1.0, 0.0, 2.5),   # Cys-S
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[20, 21],
        cutoff_a=4.0,
    )
    assert accs == []


def test_carbon_not_acceptor():
    """C (atomic_num=6) ist KEIN H-Bond-Akzeptor (nur N, O, F, S
    erlaubt)."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1, 0.0, 0.0, 0.0),
        _atom(2, 7, 0.0, 0.0, 1.0),
        _atom(3, 6, 0.0, 0.0, 2.5),     # C in Reichweite (atomic_num=6)
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    assert accs == []


def test_iron_not_acceptor():
    """Fe (atomic_num=26) ist kein Akzeptor (auch wenn nicht im
    cluster_centers gelistet)."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1, 0.0, 0.0, 0.0),
        _atom(2, 7, 0.0, 0.0, 1.0),
        _atom(3, 26, 0.0, 0.0, 2.5),    # einsames Fe
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    assert accs == []


# =============================================================================
# Element-Filter: alle erlaubten Akzeptor-Elemente
# =============================================================================

def test_all_acceptor_elements_recognized():
    """Alle vier erlaubten Elemente (N, O, F, S) werden als Akzeptoren
    erkannt."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1,  1, 0.0, 0.0, 0.0),    # H
        _atom(2,  7, 0.0, 0.0, 1.0),    # donor-N
        _atom(3,  7, 0.0, 0.0, 2.5),    # N-Akzeptor
        _atom(4,  8, 0.0, 1.0, 2.5),    # O-Akzeptor
        _atom(5,  9, 1.0, 0.0, 2.5),    # F-Akzeptor
        _atom(6, 16, 1.0, 1.0, 2.5),    # S-Akzeptor
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    centers = sorted(c for c, _ in accs)
    assert centers == [3, 4, 5, 6]


# =============================================================================
# Sortierung
# =============================================================================

def test_acceptors_sorted_by_distance():
    """Mehrere Akzeptoren werden nach Distanz sortiert (kleinster
    Abstand zuerst)."""
    h_pos = np.array([0.0, 0.0, 0.0])
    atoms_h = [
        _atom(1, 1,  0.0, 0.0, 0.0),
        _atom(2, 7,  0.0, 0.0, 1.0),
        _atom(3, 8,  0.0, 0.0, 3.5),    # 3.5 A
        _atom(4, 8,  0.0, 0.0, 2.0),    # 2.0 A (am naechsten)
        _atom(5, 8,  0.0, 0.0, 3.0),    # 3.0 A
    ]
    accs = find_hbond_acceptors_for_h(
        h_pos, atoms_h, h_center=1, n_donor_center=2,
        cluster_centers=[], ligand_centers=[], cutoff_a=4.0,
    )
    distances = [d for _, d in accs]
    assert distances == sorted(distances), (
        f"Distanzen nicht sortiert: {distances}")
    # Erster Eintrag: das naechste Atom (Center 4 bei 2.0 A)
    assert accs[0][0] == 4
    assert abs(accs[0][1] - 2.0) < 1e-9


# =============================================================================
# Standalone-Aufruf
# =============================================================================

if __name__ == "__main__":
    import inspect
    test_funcs = [
        (n, fn) for n, fn in inspect.getmembers(sys.modules[__name__])
        if n.startswith("test_") and callable(fn)
    ]
    print(f"Running {len(test_funcs)} tests in {os.path.basename(__file__)}\n")
    n_pass = n_fail = 0
    for n, fn in test_funcs:
        try:
            fn()
            print(f"  [OK]   {n}")
            n_pass += 1
        except Exception as e:
            print(f"  [FAIL] {n}: {e!r}")
            n_fail += 1
    print(f"\n{n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
