# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit-Tests fuer die v3.5-Helfer in core.py.

Testet:
* classify_significance + classify_difference_significance
* _oop_ring_metrics
* _bend_split (inkl. NaN-Guard fuer Prozente bei kleinem Signal)

Aufruf::

    python3 tests/test_helpers_v35.py

Erfolgreich, wenn am Ende "ALL TESTS PASSED" gedruckt wird.
"""
from __future__ import annotations
import math
import sys
from pathlib import Path

# Pfad-Setup: src/ ins sys.path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from modenanalyse_2fe2s.core import (  # noqa: E402
    classify_significance,
    classify_difference_significance,
    _oop_ring_metrics,
    _bend_split,
)


def _assert(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


# ---------------------------------------------------------------------------
# classify_significance
# ---------------------------------------------------------------------------
def test_classify_significance() -> None:
    _assert(classify_significance(0.5, 1.0) == "trivial",
            "0.5 sigma sollte trivial sein")
    _assert(classify_significance(2.0, 1.0) == "significant",
            "2 sigma sollte signifikant sein")
    _assert(classify_significance(5.0, 1.0) == "high",
            "5 sigma sollte hoch sein")
    _assert(classify_significance(0.0, 1.0) == "trivial",
            "Null-Wert ist trivial")
    _assert(classify_significance(float("nan"), 1.0) == "trivial",
            "NaN-Wert ist trivial")
    _assert(classify_significance(1.5, 0.0) == "high",
            "Sigma=0 mit Wert>0 sollte hoch sein (deterministisch)")
    print("[OK] classify_significance")


# ---------------------------------------------------------------------------
# classify_difference_significance
# ---------------------------------------------------------------------------
def test_classify_difference_significance() -> None:
    # Diff = 0.0033, sigma_diff = sqrt(2)*0.0008 = 0.001131
    # ratio = 0.0033/0.001131 = 2.92 -> 'significant'
    res = classify_difference_significance(0.0042, 0.0008, 0.0009, 0.0008)
    _assert(res == "significant", f"erwartet 'significant', got {res}")

    # Diff = 0.0035, ratio = 3.09 -> 'high'
    res = classify_difference_significance(0.0044, 0.0008, 0.0009, 0.0008)
    _assert(res == "high", f"erwartet 'high', got {res}")

    # Identische Werte -> trivial
    res = classify_difference_significance(0.001, 0.0005, 0.001, 0.0005)
    _assert(res == "trivial", f"erwartet 'trivial', got {res}")
    print("[OK] classify_difference_significance")


# ---------------------------------------------------------------------------
# _oop_ring_metrics
# ---------------------------------------------------------------------------
def test_oop_ring_metrics() -> None:
    n_hat = np.array([0.0, 0.0, 1.0])

    # Reine OOP-Bewegung
    evg = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.8]])
    m = _oop_ring_metrics(evg, n_hat, sigma_ev=0.001)
    _assert(m["oop_pct"] == 100.0, f"OOP=100, got {m['oop_pct']}")
    _assert(m["inp_pct"] == 0.0, f"INP=0, got {m['inp_pct']}")
    _assert(m["n_atoms"] == 3, f"n=3, got {m['n_atoms']}")

    # Reine INP-Bewegung
    evg = np.array([[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.7, -0.3, 0.0]])
    m = _oop_ring_metrics(evg, n_hat, sigma_ev=0.001)
    _assert(abs(m["oop_pct"]) < 1e-10, f"OOP=0, got {m['oop_pct']}")

    # Leerer Ring
    m = _oop_ring_metrics(np.zeros((0, 3)), n_hat, sigma_ev=0.001)
    _assert(m["n_atoms"] == 0, "leerer Ring: n_atoms=0")
    _assert(m["oop_pct"] == 0.0, "leerer Ring: oop=0")
    print("[OK] _oop_ring_metrics")


# ---------------------------------------------------------------------------
# _bend_split
# ---------------------------------------------------------------------------
def test_bend_split() -> None:
    bhat = np.array([1.0, 0.0, 0.0])      # Fe->Lig in x
    n_hat = np.array([0.0, 0.0, 1.0])     # Cluster-Normale in z
    d_fe = np.array([0.0, 0.0, 0.0])
    sigma_ev = 0.001

    # Reine OOP-Biegung: Lig bewegt sich in z, senkrecht zu bhat
    d_lig = np.array([0.0, 0.0, 0.05])
    s = _bend_split(d_lig, d_fe, bhat, n_hat, sigma_ev)
    _assert(abs(s["stretch"]) < 1e-10, "stretch=0 erwartet")
    _assert(abs(s["bend_inp"]) < 1e-10, f"bend_inp=0, got {s['bend_inp']}")
    _assert(abs(s["bend_oop"] - 0.05) < 1e-10, f"bend_oop=0.05, got {s['bend_oop']}")
    _assert(abs(s["bend_oop_pct"] - 100.0) < 1e-6, "100% OOP")

    # Reine INP-Biegung: Lig bewegt sich in y (senkrecht zu bhat=x und n=z)
    d_lig = np.array([0.0, 0.04, 0.0])
    s = _bend_split(d_lig, d_fe, bhat, n_hat, sigma_ev)
    _assert(abs(s["bend_oop"]) < 1e-10, f"bend_oop=0, got {s['bend_oop']}")
    _assert(abs(s["bend_inp"] - 0.04) < 1e-10, f"bend_inp=0.04, got {s['bend_inp']}")
    _assert(abs(s["bend_inp_pct"] - 100.0) < 1e-6, "100% INP")

    # Pythagoras-Mischfall: bend_inp=0.03, bend_oop=0.04, bend_tot=0.05.
    # Prozente sind Energie-Anteile (Quadrat-Anteile, konsistent mit oop_pct):
    #   bend_inp_pct = 100 * bend_inp^2 / bend_tot^2 = 100 * 0.09/0.25 = 36 %
    #   bend_oop_pct = 100 * bend_oop^2 / bend_tot^2 = 100 * 0.16/0.25 = 64 %
    # Summe ist exakt 100 % (Pythagoras), waehrend Laengen-Anteile das
    # nicht erfuellen wuerden.
    d_lig = np.array([0.0, 0.03, 0.04])
    s = _bend_split(d_lig, d_fe, bhat, n_hat, sigma_ev)
    _assert(abs(s["bend_inp"] - 0.03) < 1e-10, f"bend_inp=0.03, got {s['bend_inp']}")
    _assert(abs(s["bend_oop"] - 0.04) < 1e-10, f"bend_oop=0.04, got {s['bend_oop']}")
    _assert(abs(s["bend_inp_pct"] - 36.0) < 1e-6, f"36% inp, got {s['bend_inp_pct']}")
    _assert(abs(s["bend_oop_pct"] - 64.0) < 1e-6, f"64% oop, got {s['bend_oop_pct']}")
    _assert(abs(s["bend_inp_pct"] + s["bend_oop_pct"] - 100.0) < 1e-6,
            "Summe der Anteile = 100% (Pythagoras)")

    # NaN bei zu kleinem Signal: bend << 2*sqrt(2)*sigma_ev
    d_lig = np.array([0.0, 1e-4, 1e-4])
    s = _bend_split(d_lig, d_fe, bhat, n_hat, sigma_ev)
    _assert(math.isnan(s["bend_inp_pct"]), "NaN bei kleinem Signal")
    _assert(math.isnan(s["bend_oop_pct"]), "NaN bei kleinem Signal")

    # Stretch + Biege-Mix
    d_lig = np.array([0.02, 0.0, 0.03])  # stretch=0.02 (in bhat), oop=0.03
    s = _bend_split(d_lig, d_fe, bhat, n_hat, sigma_ev)
    _assert(abs(s["stretch"] - 0.02) < 1e-10, f"stretch=0.02, got {s['stretch']}")
    _assert(abs(s["bend_oop"] - 0.03) < 1e-10, "bend_oop=0.03")
    _assert(abs(s["bend_inp"]) < 1e-10, "bend_inp=0")

    # Sigma-Konvention: sigma der Differenz e_lig-e_fe ist sigma_ev*sqrt(2)
    expected_sigma = sigma_ev * math.sqrt(2.0)
    _assert(abs(s["sigma_bend_inp"] - expected_sigma) < 1e-15, "sigma_bend_inp")
    _assert(abs(s["sigma_bend_oop"] - expected_sigma) < 1e-15, "sigma_bend_oop")
    _assert(abs(s["sigma_stretch"] - expected_sigma) < 1e-15, "sigma_stretch")

    print("[OK] _bend_split (incl. NaN-guard, sigmas)")


# ---------------------------------------------------------------------------
# Konsistenz: bend_inp^2 + bend_oop^2 = |bend|^2 (Pythagoras)
# ---------------------------------------------------------------------------
def test_bend_split_pythagoras() -> None:
    """Verifiziere die exakte Pythagoras-Beziehung fuer beliebige Eingabe."""
    rng = np.random.default_rng(42)
    n_hat = np.array([0.0, 0.0, 1.0])
    bhat = np.array([1.0, 0.0, 0.0])
    sigma_ev = 0.001

    for _ in range(50):
        d_fe = rng.normal(0, 0.02, 3)
        d_lig = rng.normal(0, 0.02, 3)
        s = _bend_split(d_lig, d_fe, bhat, n_hat, sigma_ev)
        bend_total_from_split = math.sqrt(
            s["bend_inp"] ** 2 + s["bend_oop"] ** 2
        )
        # Direkter Vergleich: bend = |rel - (rel*bhat)bhat|
        rel = d_lig - d_fe
        rel_perp = rel - (rel @ bhat) * bhat
        bend_direct = float(np.linalg.norm(rel_perp))
        _assert(
            abs(bend_total_from_split - bend_direct) < 1e-12,
            f"Pythagoras verletzt: split={bend_total_from_split}, "
            f"direkt={bend_direct}",
        )
    print("[OK] _bend_split Pythagoras (50 Zufallsvektoren)")


# ---------------------------------------------------------------------------
def main() -> int:
    test_classify_significance()
    test_classify_difference_significance()
    test_oop_ring_metrics()
    test_bend_split()
    test_bend_split_pythagoras()
    print()
    print("ALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
