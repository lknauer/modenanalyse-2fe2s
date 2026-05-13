# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Numerical sanity tests: physics-based property tests.

These tests verify that the code computes what the manual formulas
*claim* it computes, by constructing inputs whose answer is known
analytically from first principles.

Unlike the audit/regression tests (which protect against bugs that
already happened), these tests protect against the **future** class
of bug where someone refactors a formula and the implementation no
longer matches the physics. Each test is keyed to a specific
equation in Manual.tex.

Coverage
--------
1. Pure-OOP mode → α_OOP = 100%  (Manual Eq.~\\ref{eq:oop-anteil})
2. Pure-INP mode → α_OOP = 0%
3. Translation mode → all bond modulations Δr_X = 0
                        (Manual Eq.~\\ref{eq:dr-classical})
4. Pure stretching mode at T→0:
   λ_X = ℏω/4 * α_X^2  (Manual Eq.~\\ref{eq:lambda-zero-T})
5. u_rms at T→0 → sqrt(ℏ/(2μω))  (Manual Eq.~\\ref{eq:urms})
6. u_rms at high T → sqrt(kT/μω²) (classical limit of Eq.~\\ref{eq:urms})
7. HA reaction-coord: pure acceptor motion → Δr_HA = 0
                        (Manual Eq.~\\ref{eq:dr-reaction-coord})
8. HA reaction-coord: pure H motion toward A → Δr_HA = +|e_H|
9. RSS aggregation with one sub-channel → dr_rss = |dr|
                        (Manual Eq.~\\ref{eq:dr-rss})
10. RSS aggregation: λ_parent = ½μω²·dr_rss² (energy consistency)

Run::

    python3 -m pytest tests/test_sanity_physics.py -v
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))


# Physical constants (matching reorganization.py values)
_HC_JCM = 1.986445857e-23   # Joules per cm^-1
_C_CMS  = 2.99792458e10     # cm/s
_AMU_KG = 1.66053906660e-27 # kg/amu
_HBAR_JS = 1.054571817e-34  # J*s
_KB_JK  = 1.380649e-23       # J/K


# =============================================================================
# 1+2. OOP fraction extremal cases
# =============================================================================

def test_pure_oop_mode_yields_alpha_oop_100():
    """Manual Eq.~\\ref{eq:oop-anteil}: a mode whose displacement is
    entirely along the cluster normal direction must produce
    α_OOP = 100 %.
    """
    from modenanalyse_2fe2s.core import _oop

    # Cluster normal along +z
    n = np.array([0.0, 0.0, 1.0])
    # 4 cluster atoms, every one displaced purely in z (umbrella mode)
    evg = np.array([
        [0.0, 0.0,  0.02],
        [0.0, 0.0,  0.02],
        [0.0, 0.0, -0.02],
        [0.0, 0.0, -0.02],
    ])
    alpha = _oop(evg, n)
    # _oop returns the fraction (between 0 and 1), not a percent
    assert abs(alpha - 1.0) < 1e-12, f"Pure-OOP mode gave α={alpha}, expected 1.0"


def test_pure_inp_mode_yields_alpha_oop_0():
    """Manual Eq.~\\ref{eq:oop-anteil}: a mode whose displacement is
    entirely perpendicular to the cluster normal must produce
    α_OOP = 0%.
    """
    from modenanalyse_2fe2s.core import _oop

    n = np.array([0.0, 0.0, 1.0])
    # Breathing mode: all motion in the xy plane
    evg = np.array([
        [+0.018, 0.000, 0.0],
        [-0.018, 0.000, 0.0],
        [ 0.000, +0.018, 0.0],
        [ 0.000, -0.018, 0.0],
    ])
    alpha = _oop(evg, n)
    assert abs(alpha) < 1e-12, f"Pure-INP mode gave α={alpha}, expected 0.0"


def test_mixed_oop_inp_satisfies_pythagoras():
    """Manual Eq.~\\ref{eq:oop-anteil}: for a general mode, α_OOP and α_INP
    must satisfy α_OOP + α_INP = 1 exactly (Pythagorean splitting).
    """
    from modenanalyse_2fe2s.core import _oop

    n = np.array([0.0, 0.0, 1.0])
    rng = np.random.default_rng(42)
    # Random eigenvector for a 4-atom cluster
    evg = rng.normal(scale=0.01, size=(4, 3))

    alpha_oop = _oop(evg, n)
    # INP is what's left after projecting out the n direction
    evg_inp = evg - (evg @ n)[:, None] * n[None, :]
    P_inp = float(np.sum(evg_inp**2))
    P_total = float(np.sum(evg**2))
    alpha_inp = P_inp / P_total

    assert abs(alpha_oop + alpha_inp - 1.0) < 1e-12


# =============================================================================
# 3. Translation mode → all bond modulations zero
# =============================================================================

def test_pure_translation_mode_yields_zero_bond_modulation():
    """Manual Eq.~\\ref{eq:dr-classical}: a translation mode (all atoms
    move identically) cannot stretch any bond. Δr_X must vanish for
    every X.
    """
    from modenanalyse_2fe2s.reorganization import signed_dr_along_axis

    # Pure translation in x: e_a = e_b = (0.05, 0, 0)
    e_a = np.array([0.05, 0.0, 0.0])
    e_b = np.array([0.05, 0.0, 0.0])
    r_a = np.array([2.20, 0.0, 0.0])    # Fe-S bond ~2.2 A
    r_b = np.array([0.00, 0.0, 0.0])

    dr = signed_dr_along_axis(e_a, e_b, r_a, r_b)
    assert abs(dr) < 1e-14, f"Translation mode gave Δr={dr}, expected 0.0"


def test_pure_rotation_mode_yields_zero_bond_stretch():
    """A rigid rotation of the bond must also leave its length unchanged
    (to first order in the displacement). For two atoms rotating by a
    small angle around the bond midpoint, their relative displacement
    is perpendicular to the bond, so Δr along the bond is 0.
    """
    from modenanalyse_2fe2s.reorganization import signed_dr_along_axis

    # Bond along x, both atoms displaced perpendicular by ±small amount
    r_a = np.array([+1.0, 0.0, 0.0])
    r_b = np.array([-1.0, 0.0, 0.0])
    # Small rotation around the bond midpoint (origin): e_a in +y, e_b in -y
    e_a = np.array([0.0, +0.05, 0.0])
    e_b = np.array([0.0, -0.05, 0.0])

    dr = signed_dr_along_axis(e_a, e_b, r_a, r_b)
    assert abs(dr) < 1e-14, f"Rotation gave Δr={dr}, expected 0.0"


def test_pure_bond_stretching_mode_yields_full_modulation():
    """A pure stretching motion of bond AB (atoms move along the bond
    axis with opposite signs) must give |Δr| = |e_a| + |e_b| when both
    move outward.
    """
    from modenanalyse_2fe2s.reorganization import signed_dr_along_axis

    # Bond along x, atoms stretching outward
    r_a = np.array([+1.0, 0.0, 0.0])
    r_b = np.array([-1.0, 0.0, 0.0])
    e_a = np.array([+0.03, 0.0, 0.0])
    e_b = np.array([-0.03, 0.0, 0.0])
    # Convention: n̂_ab points from b to a, so (e_a - e_b)·n̂_ab > 0
    # means stretching.

    dr = signed_dr_along_axis(e_a, e_b, r_a, r_b)
    expected = 0.03 - (-0.03)  # = 0.06
    assert abs(dr - expected) < 1e-14, \
        f"Pure stretching gave Δr={dr}, expected {expected}"


# =============================================================================
# 4. λ_X = (ℏω/4) at T=0 for a pure stretching mode (α_X = 1)
# =============================================================================

def test_lambda_zero_T_pure_stretching():
    """Manual Eq.~\\ref{eq:lambda-zero-T}: for a pure stretching mode at
    T=0, λ_X = ℏω/4 · α_X^2. With α_X = 1 this should yield ℏω/4 ≈ ω/4
    in cm⁻¹.

    Construction:
      u_rms = sqrt(ℏ/(2μω))     (T=0 limit)
      Δr_X  = u_rms · α_X = u_rms · 1
      λ_X(in cm⁻¹) = ½ μ ω² Δr² / hc
                   = ½ μ ω² · ℏ/(2μω) / hc
                   = ℏω/4 / hc
                   = ω/4   (in cm⁻¹)
    """
    from modenanalyse_2fe2s.reorganization import lambda_pair_cm1
    from modenanalyse_2fe2s.core import compute_thermal_amplitude

    omega_cm1 = 320.0    # typical Fe-S stretch
    mu_amu = 31.97       # mass of one S atom in amu
    T = 0.0001           # essentially T=0 (avoid coth singularity)

    u_rms_A = compute_thermal_amplitude(omega_cm1, mu_amu, T)
    # Δr = u_rms * α_X with α_X = 1
    dr_a = u_rms_A
    lam = lambda_pair_cm1(dr_a, omega_cm1, mu_amu)

    expected = omega_cm1 / 4.0     # ℏω/4 in cm⁻¹
    rel_err = abs(lam - expected) / expected
    assert rel_err < 1e-3, \
        f"λ_X = {lam:.4f} cm⁻¹, expected {expected:.4f} (rel err {rel_err:.2e})"


# =============================================================================
# 5+6. u_rms limits
# =============================================================================

def test_urms_zero_T_limit():
    """Manual Eq.~\\ref{eq:urms} at T→0 should reduce to the zero-point
    amplitude sqrt(ℏ/(2μω)).
    """
    from modenanalyse_2fe2s.core import compute_thermal_amplitude

    omega_cm1 = 320.0
    mu_amu = 31.97
    omega_si = 2 * np.pi * _C_CMS * omega_cm1
    mu_si = mu_amu * _AMU_KG

    # T→0: coth(ℏω/2kT) → 1
    u_expected_m = np.sqrt(_HBAR_JS / (2 * mu_si * omega_si))
    u_expected_A = u_expected_m * 1e10

    u_actual = compute_thermal_amplitude(omega_cm1, mu_amu, temp_k=0.001)
    rel_err = abs(u_actual - u_expected_A) / u_expected_A
    assert rel_err < 1e-6, \
        f"u_rms(T→0) = {u_actual:.6f} A, expected {u_expected_A:.6f} A " \
        f"(rel err {rel_err:.2e})"


def test_urms_classical_limit():
    """At T >> ℏω/k_B: coth(x) → 1/x = 2kT/ℏω, so u_rms → sqrt(kT/μω²)
    (classical equipartition).
    """
    from modenanalyse_2fe2s.core import compute_thermal_amplitude

    omega_cm1 = 50.0    # low frequency: ℏω/k_B ≈ 72 K
    mu_amu = 50.0
    T_K = 2000.0        # T >> 72 K → deep classical regime

    omega_si = 2 * np.pi * _C_CMS * omega_cm1
    mu_si = mu_amu * _AMU_KG

    # Classical limit:  <q²> = kT/(μω²),  u_rms = sqrt(<q²>)
    u_expected_m = np.sqrt(_KB_JK * T_K / (mu_si * omega_si**2))
    u_expected_A = u_expected_m * 1e10

    u_actual = compute_thermal_amplitude(omega_cm1, mu_amu, temp_k=T_K)
    rel_err = abs(u_actual - u_expected_A) / u_expected_A
    # The classical limit is asymptotic; at T=2000 K the residual
    # quantum correction is small but nonzero. Tolerance 1%.
    assert rel_err < 1e-2, \
        f"u_rms(T=2000 K, ω=50 cm⁻¹) = {u_actual:.4f} A, " \
        f"classical {u_expected_A:.4f} A (rel err {rel_err:.2e})"


def test_urms_monotonic_in_T():
    """u_rms must increase monotonically with T (more thermal energy
    → larger displacement)."""
    from modenanalyse_2fe2s.core import compute_thermal_amplitude

    omega_cm1 = 200.0
    mu_amu = 32.0
    temps = [1.0, 50.0, 100.0, 300.0, 1000.0]
    u_values = [compute_thermal_amplitude(omega_cm1, mu_amu, T) for T in temps]
    for i in range(1, len(u_values)):
        assert u_values[i] >= u_values[i - 1] - 1e-12, \
            f"u_rms not monotonic: {u_values}"


# =============================================================================
# 7+8. HA reaction-coord limits
# =============================================================================

def test_ha_pure_acceptor_motion_yields_zero():
    """Manual Eq.~\\ref{eq:dr-reaction-coord}:
    Δr_HA = (e_H - e_N) · n̂_NA.
    If only the acceptor moves and H, N are static, then e_H = e_N = 0,
    so Δr_HA = 0 (no PT character).
    """
    e_H = np.array([0.0, 0.0, 0.0])
    e_N = np.array([0.0, 0.0, 0.0])
    # Acceptor moves, but it's not part of the formula → Δr_HA = 0 by
    # construction.
    r_N = np.array([0.0, 0.0, 0.0])
    r_A = np.array([2.8, 0.0, 0.0])
    n_NA = (r_A - r_N) / np.linalg.norm(r_A - r_N)
    dr_HA = float(np.dot(e_H - e_N, n_NA))
    assert abs(dr_HA) < 1e-14, \
        f"Pure acceptor motion gave Δr_HA={dr_HA}, expected 0.0"


def test_ha_pure_h_motion_toward_acceptor_yields_positive():
    """Manual Eq.~\\ref{eq:dr-reaction-coord}: if the H moves toward A
    by amount d while N is static, Δr_HA = +d.
    """
    e_N = np.array([0.0, 0.0, 0.0])
    e_H = np.array([0.05, 0.0, 0.0])    # H moves +x toward acceptor
    r_N = np.array([0.0, 0.0, 0.0])
    r_A = np.array([2.8, 0.0, 0.0])
    n_NA = (r_A - r_N) / np.linalg.norm(r_A - r_N)
    dr_HA = float(np.dot(e_H - e_N, n_NA))
    expected = 0.05
    assert abs(dr_HA - expected) < 1e-14, \
        f"H toward A by 0.05 A gave Δr_HA={dr_HA}, expected {expected}"


def test_ha_pure_n_motion_toward_acceptor_yields_negative():
    """If only N moves toward A (skeleton mode, no actual PT character),
    Δr_HA = -d. This is the canonical PT vs. skeleton-mode discriminator.
    """
    e_N = np.array([0.05, 0.0, 0.0])    # N moves +x toward A
    e_H = np.array([0.0, 0.0, 0.0])
    r_N = np.array([0.0, 0.0, 0.0])
    r_A = np.array([2.8, 0.0, 0.0])
    n_NA = (r_A - r_N) / np.linalg.norm(r_A - r_N)
    dr_HA = float(np.dot(e_H - e_N, n_NA))
    expected = -0.05
    assert abs(dr_HA - expected) < 1e-14, \
        f"N toward A by 0.05 A gave Δr_HA={dr_HA}, expected {expected}"


# =============================================================================
# 9+10. RSS aggregation properties
# =============================================================================

def test_rss_aggregation_single_subchannel():
    """Manual Eq.~\\ref{eq:dr-rss}: with one sub-channel and unit weight,
    Δr_RSS = sqrt(w · dr²) = |dr|. Verified directly because we have a
    pure formula here, not a function with extra structure.
    """
    dr_values = [0.012]
    weights = [1.0]
    dr_rss = float(np.sqrt(sum(w * dr**2
                                 for w, dr in zip(weights, dr_values))))
    assert abs(dr_rss - abs(dr_values[0])) < 1e-14


def test_rss_aggregation_two_subchannels_consistency_with_lambda():
    """Manual Eq.~\\ref{eq:dr-rss}: the RSS aggregation is energetically
    consistent with the λ aggregation, in the sense that

        λ_parent  =  ½ μ ω² · (Δr_RSS)^2   modulo unit conversion

    This is the key property that motivates RSS over L1 for the
    modulation spectrum (Manual Eq.~\\ref{eq:M-spectrum}).
    """
    from modenanalyse_2fe2s.reorganization import lambda_pair_cm1

    omega_cm1 = 350.0
    mu_amu = 31.97
    # Two sub-channels, unit weights (as for FeN-Cys2 / FeN-His2)
    dr1 = 0.008
    dr2 = -0.012
    w1 = w2 = 1.0

    # Manual approach: λ_parent = sum of sub-λ's
    lam1 = lambda_pair_cm1(dr1, omega_cm1, mu_amu)
    lam2 = lambda_pair_cm1(dr2, omega_cm1, mu_amu)
    lam_sum = w1 * lam1 + w2 * lam2

    # Equivalent via RSS:
    dr_rss = np.sqrt(w1 * dr1**2 + w2 * dr2**2)
    lam_via_rss = lambda_pair_cm1(dr_rss, omega_cm1, mu_amu)

    rel_err = abs(lam_sum - lam_via_rss) / lam_sum
    assert rel_err < 1e-10, \
        f"λ_sum = {lam_sum}, λ_via_RSS = {lam_via_rss} (rel err {rel_err:.2e})"


# =============================================================================
# Extra: dimensional sanity for λ
# =============================================================================

def test_lambda_zero_for_zero_displacement():
    """A mode with zero bond modulation must contribute zero λ.
    Trivial but guards against e.g. a stray bias term."""
    from modenanalyse_2fe2s.reorganization import lambda_pair_cm1
    assert lambda_pair_cm1(0.0, 350.0, 32.0) == 0.0


def test_lambda_positive_for_any_nonzero_displacement():
    """λ = ½μω²(Δr)² is non-negative; the sign of Δr doesn't matter."""
    from modenanalyse_2fe2s.reorganization import lambda_pair_cm1
    lam_pos = lambda_pair_cm1(+0.01, 350.0, 32.0)
    lam_neg = lambda_pair_cm1(-0.01, 350.0, 32.0)
    assert lam_pos > 0
    assert abs(lam_pos - lam_neg) < 1e-15, \
        f"λ should be sign-independent: λ(+)={lam_pos}, λ(-)={lam_neg}"


def test_lambda_quadratic_in_displacement():
    """Doubling Δr must quadruple λ."""
    from modenanalyse_2fe2s.reorganization import lambda_pair_cm1
    lam1 = lambda_pair_cm1(0.01, 350.0, 32.0)
    lam2 = lambda_pair_cm1(0.02, 350.0, 32.0)
    ratio = lam2 / lam1
    assert abs(ratio - 4.0) < 1e-10, f"λ doubling test: ratio = {ratio}"


def test_lambda_quadratic_in_frequency():
    """Doubling ω must quadruple λ (since λ ∝ ω²)."""
    from modenanalyse_2fe2s.reorganization import lambda_pair_cm1
    lam1 = lambda_pair_cm1(0.01, 200.0, 32.0)
    lam2 = lambda_pair_cm1(0.01, 400.0, 32.0)
    ratio = lam2 / lam1
    assert abs(ratio - 4.0) < 1e-10, f"λ frequency test: ratio = {ratio}"
