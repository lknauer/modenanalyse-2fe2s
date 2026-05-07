# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit-Tests fuer die zentralen physikalischen Funktionen in core.py.

Testet:
* compute_thermal_amplitude (QHO-Formel mit coth-Faktor):
  - klassisches und Quanten-Limit
  - Edge-Cases (T=None, T=0, freq=0, mass=0)
  - bekannte Referenzwerte
* classify_oop_inp (binaere und 7-Stufen-Klassifikation):
  - Symmetrie um 50 %
  - Schwellen-Edge-Cases
  - Default-Schwellen vs. Custom-Schwellen

Aufruf::

    python -m pytest tests/test_core_physics.py -v
    python tests/test_core_physics.py
"""
import math
import os
import sys

import numpy as np
import pytest

# Pfad-Setup fuer direkten Aufruf
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "src"))

from modenanalyse_2fe2s.core import (
    compute_thermal_amplitude,
    classify_oop_inp,
)


# Physikalische Konstanten fuer manuelle Referenz-Rechnungen
_HBAR = 1.054571817e-34   # J*s
_KB   = 1.380649e-23      # J/K
_C    = 2.99792458e10     # cm/s
_AMU  = 1.66053906660e-27 # kg


# =============================================================================
# compute_thermal_amplitude - Edge-Cases
# =============================================================================

def test_thermal_amplitude_temp_none_returns_fallback():
    """temp_k=None aktiviert den klassischen Modus mit Fallback-Amplitude."""
    u = compute_thermal_amplitude(freq_cm1=200.0, red_mass_amu=20.0,
                                   temp_k=None, amplitude=0.5)
    assert u == 0.5


def test_thermal_amplitude_temp_zero_returns_fallback():
    """temp_k=0 zaehlt als 'kein QHO' und liefert Fallback."""
    u = compute_thermal_amplitude(freq_cm1=200.0, red_mass_amu=20.0,
                                   temp_k=0.0, amplitude=0.7)
    assert u == 0.7


def test_thermal_amplitude_temp_negative_returns_fallback():
    """Negative Temperatur (unphysikalisch) liefert Fallback."""
    u = compute_thermal_amplitude(freq_cm1=200.0, red_mass_amu=20.0,
                                   temp_k=-100.0, amplitude=1.5)
    assert u == 1.5


def test_thermal_amplitude_zero_freq_returns_fallback():
    """freq_cm1=0 (translations-/rotations-aehnliche Mode) -> Fallback."""
    u = compute_thermal_amplitude(freq_cm1=0.0, red_mass_amu=20.0,
                                   temp_k=300.0, amplitude=2.0)
    assert u == 2.0


def test_thermal_amplitude_zero_mass_returns_fallback():
    """red_mass_amu=0 -> Fallback (waere physikalisch sinnlos)."""
    u = compute_thermal_amplitude(freq_cm1=200.0, red_mass_amu=0.0,
                                   temp_k=300.0, amplitude=3.0)
    assert u == 3.0


# =============================================================================
# compute_thermal_amplitude - Quantum-Limit (T -> 0 mit hoher Frequenz)
# =============================================================================

def test_thermal_amplitude_quantum_limit_zero_point():
    """Bei sehr niedriger Temperatur und hoher Frequenz (hbar*omega >>
    k_B*T): coth -> 1, also u_rms^2 = hbar/(2*m*omega).

    Test: omega = 1000 cm^-1, m = 1 amu, T = 1 K.
    hbar*omega/(2*k_B*T) ~ 720 -> deutlich im Quanten-Limit.
    """
    freq = 1000.0
    mass = 1.0
    T    = 1.0
    u = compute_thermal_amplitude(freq, mass, T)

    # Erwartet: u_rms = sqrt(hbar/(2*m*omega))
    omega_si = 2.0 * math.pi * _C * freq
    m_si     = mass * _AMU
    u_expected_m = math.sqrt(_HBAR / (2.0 * m_si * omega_si))
    u_expected_a = u_expected_m * 1e10

    assert abs(u - u_expected_a) / u_expected_a < 1e-6, (
        f"Quantum-Limit: u={u}, erwartet {u_expected_a}")


def test_thermal_amplitude_quantum_limit_does_not_diverge():
    """Sehr hohe Frequenz, sehr niedrige Temperatur (x > 500):
    Code muss coth = 1 setzen statt zu ueberlaufen. Resultat finite."""
    u = compute_thermal_amplitude(freq_cm1=4000.0, red_mass_amu=1.0,
                                   temp_k=1.0)
    assert math.isfinite(u)
    assert u > 0
    # Plausibilitaets-Bound: u_rms einer 4000 cm^-1 Mode bei 1 amu
    # ist O(0.07 A); kein Grund auch nur in die Naehe von 1 A zu kommen
    assert u < 0.2


# =============================================================================
# compute_thermal_amplitude - Klassisches Limit (T sehr gross)
# =============================================================================

def test_thermal_amplitude_classical_limit():
    """Bei hoher Temperatur (k_B*T >> hbar*omega): coth(x) -> 1/x,
    also u_rms^2 -> k_B*T/(m*omega^2). Das ist das klassische
    Aequipartitions-Limit.

    Test: omega = 100 cm^-1, m = 50 amu, T = 5000 K.
    hbar*omega/(2*k_B*T) ~ 0.014 -> deutlich im klassischen Limit.
    """
    freq = 100.0
    mass = 50.0
    T    = 5000.0
    u = compute_thermal_amplitude(freq, mass, T)

    omega_si = 2.0 * math.pi * _C * freq
    m_si     = mass * _AMU
    u_classical_sq = _KB * T / (m_si * omega_si**2)
    u_classical_m = math.sqrt(u_classical_sq)
    u_classical_a = u_classical_m * 1e10

    # Im klassischen Limit ist die QHO-Formel innerhalb weniger Prozent
    assert abs(u - u_classical_a) / u_classical_a < 0.05, (
        f"Klassisches Limit: u={u:.4f}, erwartet ~{u_classical_a:.4f} A")


def test_thermal_amplitude_increases_with_temperature():
    """Monotonie: u_rms steigt mit T (oberhalb der ZPE-Schwelle)."""
    u_low  = compute_thermal_amplitude(200.0, 20.0, 50.0)
    u_mid  = compute_thermal_amplitude(200.0, 20.0, 300.0)
    u_high = compute_thermal_amplitude(200.0, 20.0, 1000.0)
    assert u_low <= u_mid <= u_high


def test_thermal_amplitude_decreases_with_mass():
    """Monotonie: u_rms ~ 1/sqrt(m), also schwerere Atome haben
    kleinere thermische Amplitude bei gleicher Frequenz."""
    u_light = compute_thermal_amplitude(200.0, 1.0,  300.0)
    u_heavy = compute_thermal_amplitude(200.0, 50.0, 300.0)
    assert u_heavy < u_light
    # Speziell: u_rms^2 ~ 1/m -> Verhaeltnis ~ sqrt(50)
    ratio = u_light / u_heavy
    assert 6.5 < ratio < 7.5, f"Verhaeltnis {ratio}, erwartet ~sqrt(50)=7.07"


def test_thermal_amplitude_decreases_with_freq_in_quantum_limit():
    """Im Quantum-Limit: u_rms^2 = hbar/(2m*omega), also u_rms ~
    1/sqrt(omega). Hoehere Frequenz -> kleinere Amplitude."""
    u_low_omega  = compute_thermal_amplitude(100.0,  1.0, 1.0)
    u_high_omega = compute_thermal_amplitude(1000.0, 1.0, 1.0)
    assert u_high_omega < u_low_omega
    # Verhaeltnis ~ sqrt(10)
    ratio = u_low_omega / u_high_omega
    assert 2.8 < ratio < 3.4, f"Verhaeltnis {ratio}, erwartet ~sqrt(10)=3.16"


# =============================================================================
# compute_thermal_amplitude - Anwendungsnahe Werte (FeFe, NRVS bei 5K)
# =============================================================================

def test_thermal_amplitude_fefe_5K():
    """FeFe-Atmungsmode (~250 cm^-1) bei mu_FeFe ~28 amu, T=5K
    (Tieftemperatur-NRVS): typisch u_rms ~ 0.04 A.
    """
    mu_fefe = (55.845 * 55.845) / (2 * 55.845)  # ~27.92 amu
    u = compute_thermal_amplitude(250.0, mu_fefe, 5.0)
    # Nicht zu klein, nicht zu gross fuer FeFe-Atmung bei Tieftemp.
    assert 0.02 < u < 0.08, f"u={u} A erscheint unphysikalisch"


def test_thermal_amplitude_h_low_freq():
    """H-Atom in einer Low-Frequency-Mode (~200 cm^-1, mu=1 amu, 5K).
    Da H sehr leicht ist, sind Amplituden hier substantieller."""
    u = compute_thermal_amplitude(200.0, 1.0, 5.0)
    # H-Atom in 200 cm^-1 Mode bei 5K: u_rms ~ 0.15-0.3 A
    assert 0.10 < u < 0.40, f"u={u} A erscheint unphysikalisch fuer H"


# =============================================================================
# classify_oop_inp - Edge-Cases und Symmetrie
# =============================================================================

def test_classify_pure_oop_100():
    """OOP=100%: Pur OOP, Out-of-plane."""
    broad, detail = classify_oop_inp(100.0)
    assert broad == "Out-of-plane"
    assert detail == "Pure OOP"


def test_classify_pure_inp_0():
    """OOP=0%: Pur INP, In-plane."""
    broad, detail = classify_oop_inp(0.0)
    assert broad == "In-plane"
    assert detail == "Pure INP"


def test_classify_mixed_50():
    """OOP=50%: Gemischt + Torsional/Mixed."""
    broad, detail = classify_oop_inp(50.0)
    assert broad == "Torsional/Mixed"
    assert detail == "Mixed"


def test_classify_strong_oop_80():
    """OOP=80% (zwischen mid=75 und high=90): Stark OOP, Out-of-plane."""
    broad, detail = classify_oop_inp(80.0)
    assert broad == "Out-of-plane"
    assert detail == "Strong OOP"


def test_classify_strong_inp_20():
    """OOP=20% (also INP=80%): Stark INP, In-plane (per Symmetrie)."""
    broad, detail = classify_oop_inp(20.0)
    assert broad == "In-plane"
    assert detail == "Strong INP"


def test_classify_threshold_exact_60():
    """OOP=60% (genau auf binary-Schwelle): Out-of-plane,
    Mehrheitlich OOP."""
    broad, detail = classify_oop_inp(60.0)
    assert broad == "Out-of-plane"
    assert detail == "Majority OOP"


def test_classify_just_below_threshold_59():
    """OOP=59%: knapp drunter -> Torsional/Mixed (broad), Gemischt
    (detail; weil 59 < 60 = low)."""
    broad, detail = classify_oop_inp(59.0)
    assert broad == "Torsional/Mixed"
    assert detail == "Mixed"


def test_classify_symmetry_around_50():
    """Detail-Label symmetrisch um 50 % (OOP-INP-Spiegelung)."""
    pairs = [
        (95.0, "Pure OOP",          5.0,  "Pure INP"),
        (80.0, "Strong OOP",        20.0, "Strong INP"),
        (65.0, "Majority OOP", 35.0, "Majority INP"),
    ]
    for oop, dexp, inp, iexp in pairs:
        _, det_o = classify_oop_inp(oop)
        _, det_i = classify_oop_inp(inp)
        assert det_o == dexp, f"OOP={oop}: bekam {det_o}, erwartet {dexp}"
        assert det_i == iexp, f"OOP={inp}: bekam {det_i}, erwartet {iexp}"


def test_classify_custom_thresholds():
    """Custom-Schwellen werden respektiert."""
    # Schwellen (50, 70, 85) statt (60, 75, 90)
    broad, detail = classify_oop_inp(72.0,
                                       binary_threshold=50.0,
                                       detail_thresholds=(50.0, 70.0, 85.0))
    assert broad == "Out-of-plane"  # 72% >= 50%
    assert detail == "Strong OOP"    # 70% <= 72% < 85%


def test_classify_custom_low_threshold():
    """Mit niedrigerer binary_threshold (z.B. 50%) wird OOP-Klassifikation
    weniger streng."""
    # OOP = 55%
    # Mit Default (60%): Torsional/Mixed
    # Mit binary=50%: Out-of-plane
    broad_default, _ = classify_oop_inp(55.0)
    broad_loose, _   = classify_oop_inp(55.0, binary_threshold=50.0)
    assert broad_default == "Torsional/Mixed"
    assert broad_loose   == "Out-of-plane"


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
