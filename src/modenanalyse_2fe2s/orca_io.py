# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Part of modenanalyse_2fe2s — see LICENSE in the project root.

"""
ORCA format adapter (standalone, no external dependencies beyond numpy).

This module implements its own ORCA ``.hess`` parser plus the bridge
functions that convert ORCA data into the runner-internal format.

Format specification (ORCA 5.x, 6.x)
------------------------------------
A ``.hess`` file besteht from with Dollarzeichen gekennzeichneten
blocksn, beendet through ``$end``. Wir parsen folgende blocks:

* ``$atoms`` -- atom list (element, Masse in amu, Position in Bohr).
* ``$vibrational_frequencies`` -- Frequenzen in cm$^{-1}$.
* ``$normal_modes`` -- 3N x M Matrix the Cartesian-unit-Eigenvektoren,
  ausgiven in 5er-Spaltenbloecken.

Eigenvektor-Konvention
----------------------
ORCA ``.hess`` liefert Cartesian-unit-Eigenvectors ($\\sum_j |l_j|^2 = 1$),
identical to Gaussian-hpmodes-Konvention. Sie are also direkt
verwendbar.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

from .logio import BlockInfo

# Bohr -> Angstroem (CODATA 2018: a0 = 5.29177210903e-11 m)
_BOHR_TO_ANG = 0.5291772109030

# element-Symbol -> Atomnummer (deckt alle in [2Fe-2S]-Proteinen
# and gaengigen Co-Faktoren vorkommenden elements ab)
_SYM_TO_Z = {
    "H": 1,  "He": 2,  "Li": 3,  "Be": 4,  "B": 5,   "C": 6,   "N": 7,
    "O": 8,  "F": 9,   "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14,
    "P": 15, "S": 16,  "Cl": 17, "Ar": 18, "K": 19,  "Ca": 20, "Sc": 21,
    "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
    "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35,
    "Mo": 42, "Cd": 48, "I": 53,  "W": 74,
}


@dataclass
class OrcaHessResult:
    """Vereinfachtes Resultat-Container for ORCA ``.hess``-Daten."""
    atoms:           List[Dict]
    frequencies_cm1: np.ndarray
    eigenvectors:    np.ndarray
    n_atoms:         int
    n_modes:         int
    source_path:     str = ""


def is_orca_input(path: str) -> bool:
    """Checks whether ``path`` als ORCA ``.hess`` file behandelt are soll."""
    p = Path(path)
    ext = p.suffix.lower()
    if ext == ".hess":
        return True
    if ext in {".log", ".out", ".dat"}:
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as fh:
                head = fh.read(8192)
            if "$orca_hessian_file" in head or "$vibrational_frequencies" in head:
                return True
        except OSError:
            pass
    return False


def _split_blocks(text: str) -> Dict[str, str]:
    """Spaltet ``.hess`` text in seine ``$``-blocks."""
    blocks: Dict[str, str] = {}
    current_name: Optional[str] = None
    current_body: List[str] = []
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("$"):
            if current_name is not None:
                blocks[current_name] = "\n".join(current_body)
            tag = stripped.split()[0]
            if tag == "$end":
                current_name = None
                current_body = []
            else:
                current_name = tag
                current_body = []
        elif current_name is not None:
            current_body.append(line)
    if current_name is not None:
        blocks[current_name] = "\n".join(current_body)
    return blocks


def _parse_atoms_block(body: str) -> List[Dict]:
    lines = [ln.strip() for ln in body.strip().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("ORCA .hess: $atoms block empty")
    try:
        n_atoms = int(lines[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"ORCA .hess: first line of $atoms is not a number: "
                         f"{lines[0]!r}") from exc
    if len(lines) < 1 + n_atoms:
        raise ValueError(f"ORCA .hess: $atoms states {n_atoms} atoms, "
                         f"but only {len(lines)-1} data rows.")
    atoms = []
    for i in range(n_atoms):
        toks = lines[1 + i].split()
        if len(toks) < 5:
            raise ValueError(f"ORCA .hess: atom line {i} incomplete: "
                             f"{lines[1+i]!r}")
        elem = toks[0]
        try:
            mass = float(toks[1])
            xb = float(toks[2]); yb = float(toks[3]); zb = float(toks[4])
        except ValueError as exc:
            raise ValueError(f"ORCA .hess: could not parse atom line {i}: "
                             f"{lines[1+i]!r}") from exc
        atoms.append({
            "element":    elem,
            "mass_amu":   mass,
            "x_ang":      xb * _BOHR_TO_ANG,
            "y_ang":      yb * _BOHR_TO_ANG,
            "z_ang":      zb * _BOHR_TO_ANG,
            "atomic_num": _SYM_TO_Z.get(elem, 0),
        })
    return atoms


def _parse_frequencies_block(body: str) -> np.ndarray:
    lines = [ln.strip() for ln in body.strip().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("ORCA .hess: $vibrational_frequencies empty")
    try:
        n_modes = int(lines[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"ORCA .hess: first line of "
                         f"$vibrational_frequencies is not a number: "
                         f"{lines[0]!r}") from exc
    freqs = np.zeros(n_modes, dtype=float)
    for i in range(n_modes):
        toks = lines[1 + i].split()
        if len(toks) < 2:
            raise ValueError(f"ORCA .hess: Frequenzzeile {i} uncomplete: "
                             f"{lines[1+i]!r}")
        try:
            freqs[i] = float(toks[1])
        except ValueError as exc:
            raise ValueError(f"ORCA .hess: could not parse frequency {i}: "
                             f"{lines[1+i]!r}") from exc
    return freqs


def _parse_normal_modes_block(body: str) -> np.ndarray:
    """Parses the 5-column block matrix in ``$normal_modes``."""
    lines = [ln.rstrip() for ln in body.strip().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("ORCA .hess: $normal_modes empty")
    toks = lines[0].split()
    if len(toks) < 2:
        raise ValueError(f"ORCA .hess: $normal_modes Header uncomplete: "
                         f"{lines[0]!r}")
    try:
        n_rows = int(toks[0])
        n_cols = int(toks[1])
    except ValueError as exc:
        raise ValueError(f"ORCA .hess: $normal_modes header not numbers: "
                         f"{lines[0]!r}") from exc
    matrix = np.zeros((n_rows, n_cols), dtype=float)

    cursor = 1
    while cursor < len(lines):
        col_header_toks = lines[cursor].split()
        try:
            col_indices = [int(t) for t in col_header_toks]
        except ValueError:
            break
        cursor += 1
        for r in range(n_rows):
            if cursor >= len(lines):
                raise ValueError(f"ORCA .hess: $normal_modes vorzeitig "
                                 f"ended at row {r}")
            data_toks = lines[cursor].split()
            cursor += 1
            if len(data_toks) < 1 + len(col_indices):
                raise ValueError(f"ORCA .hess: $normal_modes row {r} "
                                 f"uncomplete: {lines[cursor-1]!r}")
            try:
                row_idx = int(data_toks[0])
            except ValueError:
                row_idx = r
            for k, c in enumerate(col_indices):
                if c < n_cols and row_idx < n_rows:
                    matrix[row_idx, c] = float(data_toks[1 + k])
    return matrix


def load_orca_hess(filepath: str) -> OrcaHessResult:
    """Reads an ORCA ``.hess`` file in ein ``OrcaHessResult``."""
    p = Path(filepath)
    if not p.is_file():
        raise FileNotFoundError(f"ORCA .hess file not found: {p}")

    text = p.read_text(encoding="utf-8", errors="replace")
    blocks = _split_blocks(text)

    for required in ("$atoms", "$vibrational_frequencies", "$normal_modes"):
        if required not in blocks:
            raise ValueError(
                f"ORCA .hess {str(p)!r}: required block {required!r} fehlt. "
                f"Gefundene blocks: {sorted(blocks.keys())}")

    atoms = _parse_atoms_block(blocks["$atoms"])
    freqs = _parse_frequencies_block(blocks["$vibrational_frequencies"])
    eigvecs = _parse_normal_modes_block(blocks["$normal_modes"])

    n_atoms = len(atoms)
    n_modes = freqs.size
    if eigvecs.shape != (3 * n_atoms, n_modes):
        raise ValueError(
            f"ORCA .hess {str(p)!r}: Eigenvektor-Shape {eigvecs.shape} "
            f"does not match (3*N, M) = ({3*n_atoms}, {n_modes})")

    return OrcaHessResult(
        atoms           = atoms,
        frequencies_cm1 = freqs,
        eigenvectors    = eigvecs,
        n_atoms         = n_atoms,
        n_modes         = n_modes,
        source_path     = str(p),
    )


# ---------------------------------------------------------------------------
# Adapter to runner-internen Format
# ---------------------------------------------------------------------------

def parseresult_to_atoms(orca_res: OrcaHessResult,
                         include_hydrogen: bool = True
                         ) -> Tuple[List[Dict], Dict[int, int]]:
    """Converts ``OrcaHessResult.atoms`` ins runner-interne Atomformat."""
    atoms_out: List[Dict] = []
    idx_map: Dict[int, int] = {}
    for original_idx, a in enumerate(orca_res.atoms):
        z = a["atomic_num"]
        if not include_hydrogen and z == 1:
            continue
        center = original_idx + 1   # 1-basiert wie Gaussian
        elem_str = str(a["element"])
        atoms_out.append({
            "center":     center,
            "atomic_num": z,
            "x":          float(a["x_ang"]),
            "y":          float(a["y_ang"]),
            "z":          float(a["z_ang"]),
            "mass":       float(a["mass_amu"]),
            "element":    elem_str,
            "symbol":     elem_str,    # Alias for Gauss-Kompatibilitaet
        })
        idx_map[center] = len(atoms_out) - 1
    return atoms_out, idx_map


def parseresult_to_blocks(orca_res: OrcaHessResult
                          ) -> Tuple[List[BlockInfo],
                                     Dict[int, BlockInfo],
                                     Dict[int, List[BlockInfo]]]:
    """Synthetisiert einen ORCA-Pseudoblock with allen Modes."""
    n_modes = orca_res.n_modes
    freqs = list(map(float, orca_res.frequencies_cm1))
    block = BlockInfo(
        offset      = 0,
        mode_nums   = list(range(1, n_modes + 1)),
        freqs       = freqs,
        red_masses  = [1.0] * n_modes,
        frc_consts  = [0.0] * n_modes,
        syms        = ["A"] * n_modes,
        is_hp       = True,
        data_offset = 0,
    )
    all_blocks: List[BlockInfo] = [block]
    best_block: Dict[int, BlockInfo] = {mn: block for mn in block.mode_nums}
    cand_map:   Dict[int, List[BlockInfo]] = {mn: [block] for mn in block.mode_nums}
    return all_blocks, best_block, cand_map


def get_eigvec_orca(orca_res: OrcaHessResult,
                    mode_num: int,
                    atoms: List[Dict],
                    idx_map: Dict[int, int],
                    include_hydrogen: bool = True
                    ) -> Tuple[List[int], np.ndarray]:
    """Reads the eigenvector of a single mode from the ORCA result."""
    col = mode_num - 1
    if not (0 <= col < orca_res.n_modes):
        raise IndexError(f"mode_num {mode_num} outside "
                         f"[1, {orca_res.n_modes}]")
    ev_full = orca_res.eigenvectors[:, col]
    ev_3d = ev_full.reshape(orca_res.n_atoms, 3)

    n_use = len(atoms)
    centers: List[int] = []
    evg = np.zeros((n_use, 3), dtype=float)
    for i, a in enumerate(atoms):
        orig_idx = a["center"] - 1
        if 0 <= orig_idx < orca_res.n_atoms:
            evg[i] = ev_3d[orig_idx]
        centers.append(a["center"])
    return centers, evg
