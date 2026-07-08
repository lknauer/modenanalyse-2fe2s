# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regressionstest fuer die ORCA-Reduzierte-Massen-Rekonstruktion [fix K2].

Hintergrund
-----------
ORCA-``$normal_modes`` sind einheits-kartesisch (sum_i |l_i|^2 = 1, empirisch
an realen QM/MM-.hess-Dateien bestaetigt), identisch zur Gaussian-hpmodes-
Konvention. Fuer diese Konvention ist die effektive/reduzierte Masse der
Normalkoordinate mu_k = sum_i m_i * l_{i,k}^2 -- genau der Wert, den Gaussian
im Log als "Red. mass" ausgibt. Frueher setzte ``parseresult_to_blocks``
faelschlich 1.0 amu, wodurch jede absolute ORCA-Amplitude um sqrt(mu) daneben
lag.

Dieser Test baut ein kleines synthetisches ``OrcaHessResult`` mit bekannten
Massen und Modenvektoren und prueft, dass die rekonstruierten ``red_masses``
exakt sum_i m_i l_i^2 ergeben -- ohne die (mehrere hundert MB grossen) realen
.hess-Dateien zu benoetigen.
"""

from __future__ import annotations
import os
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "src"))

from modenanalyse_2fe2s import orca_io


def _make_result(masses, eigvecs, freqs):
    """Baut ein OrcaHessResult mit vorgegebenen Massen/Moden."""
    n_atoms = len(masses)
    atoms = [{"element": "X", "atomic_num": 6,
              "mass_amu": float(masses[i]),
              "x_ang": 0.0, "y_ang": 0.0, "z_ang": 0.0}
             for i in range(n_atoms)]
    return orca_io.OrcaHessResult(
        atoms           = atoms,
        frequencies_cm1 = np.asarray(freqs, dtype=float),
        eigenvectors    = np.asarray(eigvecs, dtype=float),
        n_atoms         = n_atoms,
        n_modes         = len(freqs),
        source_path     = "<synthetic>",
    )


def test_orca_reduced_mass_equals_sum_m_lsq():
    """red_masses[k] == sum_i m_i * l_{i,k}^2 fuer einheits-normierte Moden."""
    rng = np.random.default_rng(0)
    n_atoms = 5
    masses = np.array([1.008, 12.011, 15.999, 32.06, 55.845])  # H,C,O,S,Fe
    m_dof = np.repeat(masses, 3)

    # drei zufaellige, einheits-kartesisch normierte Moden
    M = np.zeros((3 * n_atoms, 3))
    for k in range(3):
        v = rng.standard_normal(3 * n_atoms)
        M[:, k] = v / np.linalg.norm(v)      # sum |l|^2 = 1

    res = _make_result(masses, M, freqs=[150.0, 800.0, 2900.0])
    _, best, _ = orca_io.parseresult_to_blocks(res)
    red = best[1].red_masses

    expected = (m_dof[:, None] * M ** 2).sum(axis=0)
    assert np.allclose(red, expected, rtol=1e-12, atol=1e-12), \
        f"red_masses {red} != sum m|l|^2 {expected.tolist()}"

    # Physikalische Plausibilitaet: jede reduzierte Masse liegt zwischen der
    # kleinsten und groessten Atommasse.
    for r in red:
        assert masses.min() - 1e-9 <= r <= masses.max() + 1e-9


def test_orca_pure_translation_reduced_mass_is_mean_mass():
    """Reine Translation (alle Atome gleich, +x): mu = sum m_i * l_i^2.

    Fuer l_i = (1/sqrt(N)) e_x pro Atom ist mu = (1/N) sum m_i = mittlere
    Atommasse -- eine analytisch bekannte Referenz."""
    n_atoms = 4
    masses = np.array([12.0, 12.0, 16.0, 16.0])
    N = n_atoms
    v = np.zeros(3 * n_atoms)
    v[0::3] = 1.0 / np.sqrt(N)          # x-Komponente jedes Atoms
    res = _make_result(masses, v.reshape(-1, 1), freqs=[0.0])
    _, best, _ = orca_io.parseresult_to_blocks(res)
    mu = best[1].red_masses[0]
    assert mu == pytest.approx(masses.mean(), rel=1e-12)
