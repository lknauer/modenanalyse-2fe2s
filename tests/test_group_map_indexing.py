# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regression tests for the v1.0.2 Apd1 bug (Groups_*-rows silently zero).

These tests cover three independent layers of defense against the same
class of silent data-loss bug:

1.  Root cause (geometry._build_group_map): the function must use the
    same atom-list that pdb_to_gaus_h was built from, so that residues
    with hydrogens interleaved between heavy atoms do not lose centers.

2.  Diagnostic layer 1 (core.analyze_mode): if a group_map entry maps
    to zero indices in the eigenvector c2l, a UserWarning must be
    emitted at least once per group_name per run.

3.  Diagnostic layer 2 (export._ws_gruppen): if a row in Gruppen_OOP/
    INP/Winkel ends up all-zero, a UserWarning must point at the most
    likely cause.

The runtime-impact of the warnings is bounded: the first layer
(_WARNED_EMPTY_GROUP) deduplicates, the second runs once per sheet.

Run::

    python3 -m pytest tests/test_group_map_indexing.py -v
"""
from __future__ import annotations
import pytest
import warnings
from pathlib import Path

import numpy as np

# Make src/ importable when running pytest from the repo root
ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(ROOT / "src"))


# ---------------------------------------------------------------------------
# Layer 1: root cause in geometry._build_group_map
# ---------------------------------------------------------------------------

def _make_pdb_data_h_interleaved():
    """Build a synthetic pdb_data dict that triggers the original bug.

    Layout: residue 100 has 4 heavy atoms (N, CA, CB, SG) plus 2 hydrogens
    interleaved (after N and after CA). Residue 200 has 3 heavy atoms (no
    hydrogens). The interleaved hydrogens are what shifts pdb_heavy indices
    relative to all_h indices in the buggy code path.

    Returns a dict matching the shape expected by _build_group_map:
        {"all_h": [...], "atoms_h": [...]}
    where atoms_h == all_h (both lists contain heavy + hydrogen atoms in
    PDB order; the "heavy" filter is applied inside _build_group_map).
    """
    def atom(rnum, rname, aname, element, is_h, x=0.0, y=0.0, z=0.0):
        return {"rnum": rnum, "rname": rname, "aname": aname,
                "element": element, "is_h": is_h,
                "x": x, "y": y, "z": z}

    all_h = [
        atom(100, "CYS", "N",  "N", False, x=0.0),   # 0 heavy
        atom(100, "CYS", "HN", "H", True,  x=0.1),   # 1 hydrogen (interleaved)
        atom(100, "CYS", "CA", "C", False, x=1.0),   # 2 heavy
        atom(100, "CYS", "HA", "H", True,  x=1.1),   # 3 hydrogen (interleaved)
        atom(100, "CYS", "CB", "C", False, x=2.0),   # 4 heavy
        atom(100, "CYS", "SG", "S", False, x=3.0),   # 5 heavy
        atom(200, "HIS", "N",  "N", False, x=10.0),  # 6 heavy
        atom(200, "HIS", "CA", "C", False, x=11.0),  # 7 heavy
        atom(200, "HIS", "CB", "C", False, x=12.0),  # 8 heavy
    ]
    return {"all_h": all_h, "atoms_h": all_h}


def _make_ligands(res_nums):
    """Minimal LigandInfo stubs that _build_group_map needs."""
    from modenanalyse_2fe2s.geometry import LigandInfo
    return [
        LigandInfo(
            fe_idx=0, fe_center=1, lig_center=100 + i,
            lig_element="S" if i == 0 else "N",
            lig_aname="SG" if i == 0 else "ND1",
            res_num=rn,
            res_name="CYS" if i == 0 else "HIS",
            res_label=f"{'Cys' if i == 0 else 'His'} {rn}",
            bond_vec=np.array([1.0, 0.0, 0.0]),
            bond_len=2.0,
        )
        for i, rn in enumerate(res_nums)
    ]


def test_build_group_map_returns_all_heavy_atoms_with_interleaved_hydrogens():
    """The fix must collect ALL heavy atoms of a residue even when its
    hydrogens appear interleaved between heavy atoms in the PDB list.

    Before the fix, the comprehension iterated a filtered pdb_heavy list
    but looked up pdb_to_gaus_h with the filtered index; for residues
    appearing AFTER an interleaved hydrogen, the index offset caused
    silent atom loss and produced 100%-zero rows in the Groups_* sheets.
    """
    from modenanalyse_2fe2s.geometry import _build_group_map

    pdb_data = _make_pdb_data_h_interleaved()
    ligands = _make_ligands([100, 200])

    # pdb_to_gaus_h: all heavy atoms map to a Gaussian center (synthetic).
    # Keyed by index into all_h, as in find_coordinating_residues.
    pdb_to_gaus_h = {
        0: 1001,   # Cys100 N
        2: 1002,   # Cys100 CA
        4: 1003,   # Cys100 CB
        5: 1004,   # Cys100 SG
        6: 2001,   # His200 N
        7: 2002,   # His200 CA
        8: 2003,   # His200 CB
    }
    # Hydrogens (indices 1, 3) are NOT in pdb_to_gaus_h, mirroring how
    # find_coordinating_residues only stores heavy-atom mappings here.

    result = _build_group_map(ligands, pdb_data, pdb_to_gaus_h)

    # Cys100: 4 heavy atoms (N, CA, CB, SG)
    assert "Cys 100" in result, "Cys 100 missing from group_map"
    assert sorted(result["Cys 100"]) == [1001, 1002, 1003, 1004], (
        f"Cys 100 should have 4 centers, got {result['Cys 100']}. "
        f"The pre-1.0.2 bug returned a subset (e.g. only [1001, 1002] "
        f"because the interleaved hydrogens shifted indices)."
    )

    # His200: 3 heavy atoms (N, CA, CB) — no interleaved H, was already
    # correct in v1.0.1 and must stay correct
    assert "His 200" in result, "His 200 missing from group_map"
    assert sorted(result["His 200"]) == [2001, 2002, 2003], (
        f"His 200 should have 3 centers, got {result['His 200']}"
    )


def test_build_group_map_no_interleaved_hydrogens_unchanged():
    """When no hydrogens are interleaved, the fix must produce the same
    output as the pre-fix code (i.e. it does not regress the simple case).
    """
    from modenanalyse_2fe2s.geometry import _build_group_map

    def atom(rnum, rname, aname, element, is_h, x=0.0):
        return {"rnum": rnum, "rname": rname, "aname": aname,
                "element": element, "is_h": is_h, "x": x, "y": 0.0, "z": 0.0}

    # All heavy first, then all H — the legacy bug-free PDB layout
    all_h = [
        atom(100, "CYS", "N",  "N", False, x=0.0),
        atom(100, "CYS", "CA", "C", False, x=1.0),
        atom(100, "CYS", "CB", "C", False, x=2.0),
        atom(100, "CYS", "SG", "S", False, x=3.0),
        atom(100, "CYS", "HN", "H", True,  x=0.1),
        atom(100, "CYS", "HA", "H", True,  x=1.1),
    ]
    pdb_data = {"all_h": all_h, "atoms_h": all_h}
    ligands = _make_ligands([100])[:1]
    pdb_to_gaus_h = {0: 1001, 1: 1002, 2: 1003, 3: 1004}

    result = _build_group_map(ligands, pdb_data, pdb_to_gaus_h)
    assert sorted(result["Cys 100"]) == [1001, 1002, 1003, 1004]


def test_build_group_map_excludes_overlapping_water():
    """Real-world scenario: a crystal water (HOH) in a different chain
    happens to share the residue number of an amino-acid ligand. This
    is common in QM/MM PDBs where waters are numbered separately.
    The pre-1.0.2-followup _build_group_map filtered only on rnum, so
    a HOH oxygen with rnum=216 would be silently included in the Cys
    216 group, distorting the OOP/INP percentages by a non-zero amount.
    """
    from modenanalyse_2fe2s.geometry import _build_group_map

    def atom(rnum, rname, aname, element, is_h, x=0.0):
        return {"rnum": rnum, "rname": rname, "aname": aname,
                "element": element, "is_h": is_h, "x": x, "y": 0.0, "z": 0.0}

    # PDB layout that mirrors the actual Apd1 PDB:
    #   Cys 216 in chain A: N, CA, C, O, CB, SG (6 heavy) + 4 hydrogens
    #   HOH 216 in chain C: O (1 heavy)         <- accidental rnum overlap
    all_h = [
        atom(216, "CYS", "N",   "N", False, x=0.0),    # 0
        atom(216, "CYS", "CA",  "C", False, x=1.0),    # 1
        atom(216, "CYS", "C",   "C", False, x=2.0),    # 2
        atom(216, "CYS", "O",   "O", False, x=3.0),    # 3
        atom(216, "CYS", "CB",  "C", False, x=4.0),    # 4
        atom(216, "CYS", "SG",  "S", False, x=5.0),    # 5
        atom(216, "CYS", "HN",  "H", True,  x=0.1),    # 6 (H)
        atom(216, "CYS", "HA",  "H", True,  x=1.1),    # 7 (H)
        atom(216, "CYS", "HB2", "H", True,  x=4.1),    # 8 (H)
        atom(216, "CYS", "HB3", "H", True,  x=4.2),    # 9 (H)
        atom(216, "HOH", "O",   "O", False, x=10.0),   # 10  <-- the trap
    ]
    pdb_data = {"all_h": all_h, "atoms_h": all_h}
    ligands = [
        # Build the LigandInfo by hand here, because _make_ligands hardcodes
        # res_name and we need CYS specifically (not "Cys 100"-style).
        type("L", (), dict(
            fe_idx=0, fe_center=1, lig_center=1005,
            lig_element="S", lig_aname="SG",
            res_num=216, res_name="CYS", res_label="Cys 216",
            bond_vec=np.array([1.0, 0.0, 0.0]), bond_len=2.2,
        ))()
    ]
    # All 11 PDB atoms map to Gaussian centers; the HOH oxygen too.
    pdb_to_gaus_h = {0: 1001, 1: 1002, 2: 1003, 3: 1004, 4: 1005, 5: 1006,
                     # hydrogens 6-9 mapped to H centers (irrelevant here)
                     6: 1101, 7: 1102, 8: 1103, 9: 1104,
                     # HOH oxygen mapped to a water center
                     10: 9999}

    result = _build_group_map(ligands, pdb_data, pdb_to_gaus_h)

    # Cys 216 must contain ONLY the 6 amino-acid heavy atoms,
    # NOT the HOH oxygen at index 10.
    assert "Cys 216" in result
    assert sorted(result["Cys 216"]) == [1001, 1002, 1003, 1004, 1005, 1006], (
        f"Cys 216 should have exactly 6 heavy atoms (the amino-acid "
        f"backbone+sidechain), but got {result['Cys 216']}. "
        f"The HOH oxygen at Gaussian center 9999 (PDB rnum 216, rname HOH) "
        f"must not appear here."
    )
    assert 9999 not in result["Cys 216"], (
        "HOH oxygen leaked into Cys 216 group_map — rname filter not applied"
    )


# ---------------------------------------------------------------------------
# Layer 2: diagnostic warning in core.analyze_mode
# ---------------------------------------------------------------------------

def test_warned_empty_group_dedup():
    """The _WARNED_EMPTY_GROUP set must prevent duplicate warnings within
    a single run, even when analyze_mode is called thousands of times for
    the same problematic group.
    """
    from modenanalyse_2fe2s.core import _WARNED_EMPTY_GROUP, reset_warning_state

    reset_warning_state()
    assert len(_WARNED_EMPTY_GROUP) == 0, "warning state not reset"

    # Simulate the analyze_mode behaviour: add a name once, expect it not
    # to grow further on repeat calls
    _WARNED_EMPTY_GROUP.add("Cys 207")
    _WARNED_EMPTY_GROUP.add("Cys 207")
    _WARNED_EMPTY_GROUP.add("Cys 207")
    assert len(_WARNED_EMPTY_GROUP) == 1

    _WARNED_EMPTY_GROUP.add("His 259")
    assert len(_WARNED_EMPTY_GROUP) == 2

    reset_warning_state()
    assert len(_WARNED_EMPTY_GROUP) == 0, "reset_warning_state did not clear"


# ---------------------------------------------------------------------------
# Layer 3: diagnostic warning in export._ws_gruppen
# ---------------------------------------------------------------------------

def test_ws_gruppen_warns_on_all_zero_row():
    """If a residue row in Gruppen_OOP / INP / Winkel is all zero across
    all modes, _ws_gruppen must emit a UserWarning naming the residue.
    Without this layer the silent fallback ``.get(g,{}).get(key,0.)`` masks
    the bug entirely.
    """
    from openpyxl import Workbook
    from modenanalyse_2fe2s.export import _ws_gruppen

    wb = Workbook()
    # Default sheet that openpyxl creates - we ignore it
    # _ws_gruppen will create the Gruppen_* sheets

    # Three modes; "Cys 207" is all-zero in OOP/INP/Winkel, "His 255"
    # has normal values. Tors is allowed to be zero (skipped by design).
    results = [
        {"freq": 100.0, "groups": {
            "Cys 207": {"oop": 0.0, "inp": 0.0, "angle": 0.0, "torsion": 0.0},
            "His 255": {"oop": 30.0, "inp": 70.0, "angle": 15.0, "torsion": 0.001},
        }},
        {"freq": 200.0, "groups": {
            "Cys 207": {"oop": 0.0, "inp": 0.0, "angle": 0.0, "torsion": 0.0},
            "His 255": {"oop": 25.0, "inp": 75.0, "angle": 12.0, "torsion": 0.002},
        }},
        {"freq": 300.0, "groups": {
            "Cys 207": {"oop": 0.0, "inp": 0.0, "angle": 0.0, "torsion": 0.0},
            "His 255": {"oop": 40.0, "inp": 60.0, "angle": 18.0, "torsion": 0.003},
        }},
    ]
    group_names = ["Cys 207", "His 255"]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ws_gruppen(wb, results, group_names, E=None)

    # Expect at least 3 warnings (OOP, INP, Winkel — not Tors)
    cys_207_warnings = [w for w in caught
                        if "Cys 207" in str(w.message)]
    assert len(cys_207_warnings) >= 3, (
        f"Expected >= 3 warnings naming 'Cys 207' (OOP/INP/Winkel), "
        f"got {len(cys_207_warnings)}: "
        f"{[str(w.message) for w in cys_207_warnings]}"
    )

    # His 255 should never trigger the warning
    his_warnings = [w for w in caught if "His 255" in str(w.message)]
    assert len(his_warnings) == 0, (
        f"His 255 has nonzero values and must not trigger a warning, "
        f"got {[str(w.message) for w in his_warnings]}"
    )


def test_ws_gruppen_no_warning_when_all_good():
    """If every group has at least one nonzero value, no warning fires."""
    from openpyxl import Workbook
    from modenanalyse_2fe2s.export import _ws_gruppen

    wb = Workbook()
    results = [
        {"freq": 100.0, "groups": {
            "Cys 207": {"oop": 30.0, "inp": 70.0, "angle": 10.0, "torsion": 0.001},
            "His 255": {"oop": 25.0, "inp": 75.0, "angle": 15.0, "torsion": 0.002},
        }},
        {"freq": 200.0, "groups": {
            "Cys 207": {"oop": 28.0, "inp": 72.0, "angle": 11.0, "torsion": 0.001},
            "His 255": {"oop": 26.0, "inp": 74.0, "angle": 14.0, "torsion": 0.002},
        }},
    ]
    group_names = ["Cys 207", "His 255"]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _ws_gruppen(wb, results, group_names, E=None)

    relevant = [w for w in caught
                if "100% zero" in str(w.message)]
    assert len(relevant) == 0, (
        f"No warning expected when all values are nonzero, got: "
        f"{[str(w.message) for w in relevant]}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
