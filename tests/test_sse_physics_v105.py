"""v1.0.5 SSE-physics fixes: rigid-body decomposition, mass weighting,
C-alpha axis, and a regression guard against the bending-feature bug.

These tests exercise ``analyze_sse_element`` directly with hand-built
displacement fields whose translation / rotation / internal content is
known a priori, so the decomposition can be checked against ground truth.
"""
import pathlib
import numpy as np
import pytest

from modenanalyse_2fe2s.core import analyze_sse_element

SRC = pathlib.Path(__file__).resolve().parents[1] / "src" / "modenanalyse_2fe2s"


def _build(coords, disps, symbols=None):
    """Assemble (evg, c2l, sse_centers, atoms, idx_map) for n atoms.

    centers are 1..n, evg rows are 0..n-1.
    """
    coords = np.asarray(coords, float)
    disps = np.asarray(disps, float)
    n = len(coords)
    if symbols is None:
        symbols = ["C"] * n
    centers = list(range(1, n + 1))
    c2l = {c: i for i, c in enumerate(centers)}
    idx_map = {c: i for i, c in enumerate(centers)}
    atoms = [{"x": coords[i, 0], "y": coords[i, 1], "z": coords[i, 2],
              "symbol": symbols[i]} for i in range(n)]
    return disps.copy(), c2l, centers, atoms, idx_map


# straight backbone trace along z (4 "C-alpha" atoms)
Z = np.array([[0, 0, -3.0], [0, 0, -1.0], [0, 0, 1.0], [0, 0, 3.0]])


def test_pure_translation():
    t = np.array([0.2, -0.1, 0.05])
    disps = np.tile(t, (4, 1))
    evg, c2l, cen, atoms, idx = _build(Z, disps)
    r = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    assert r["com_amplitude"] == pytest.approx(np.linalg.norm(t), abs=1e-9)
    assert r["tilting_angle"] == pytest.approx(0.0, abs=1e-6)
    assert r["internal_amplitude"] == pytest.approx(0.0, abs=1e-9)


def test_pure_tilt_rotation():
    # rotation about x (perpendicular to the z helix axis): u_i = omega x r_i
    omega = np.array([0.05, 0.0, 0.0])
    disps = np.cross(omega, Z)
    evg, c2l, cen, atoms, idx = _build(Z, disps)
    r = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    # tilt is recovered as degrees(|omega_perp|); COM and internal ~ 0
    assert r["tilting_angle"] == pytest.approx(np.degrees(0.05), abs=1e-3)
    assert r["com_amplitude"] == pytest.approx(0.0, abs=1e-9)
    assert r["internal_amplitude"] == pytest.approx(0.0, abs=1e-9)


def test_pure_twist_is_not_tilt():
    # rotation ABOUT the helix axis (z): should NOT register as a tilt
    omega = np.array([0.0, 0.0, 0.05])
    # give the atoms a perpendicular offset so a z-rotation actually moves them
    coords = np.array([[1, 0, -3.0], [1, 0, -1.0], [1, 0, 1.0], [1, 0, 3.0]])
    disps = np.cross(omega, coords - coords.mean(0))
    evg, c2l, cen, atoms, idx = _build(coords, disps)
    r = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    assert r["tilting_angle"] == pytest.approx(0.0, abs=1e-3)
    assert r["internal_amplitude"] == pytest.approx(0.0, abs=1e-9)


def test_pure_internal_deformation():
    # u_i = (|z_i|, 0, 0): even in z, inexpressible as translation + rotation
    disps = np.array([[abs(z), 0.0, 0.0] for z in Z[:, 2]])
    evg, c2l, cen, atoms, idx = _build(Z, disps)
    r = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    assert r["internal_amplitude"] > 0.1
    assert r["tilting_angle"] == pytest.approx(0.0, abs=1e-6)


def test_com_is_mass_weighted():
    # two atoms, very different masses; only the light one moves
    coords = np.array([[0, 0, 0.0], [0, 0, 2.0]])
    v = np.array([0.0, 0.0, 0.3])
    disps = np.array([[0, 0, 0.0], v])
    evg, c2l, cen, atoms, idx = _build(coords, disps, symbols=["FE", "H"])
    r = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    mFe, mH = 55.845, 1.008
    expected = mH * 0.3 / (mFe + mH)          # mass-weighted COM displacement
    unweighted = 0.15                          # what the old code returned
    assert r["com_amplitude"] == pytest.approx(expected, abs=1e-6)
    assert abs(r["com_amplitude"] - unweighted) > 0.1


def test_ca_axis_overrides_sidechain_contamination():
    # short helix: 3 CA along z, plus 3 far side-chain atoms spread in x.
    # All-atom SVD axis would point along x; CA-only axis must be ~z.
    coords = np.array([
        [0, 0, -1.5], [0, 0, 0.0], [0, 0, 1.5],     # CA (centers 1,2,3)
        [8, 0, -1.5], [-8, 0, 0.0], [8, 0, 1.5],    # side chain (centers 4,5,6)
    ])
    disp_z = np.tile([0, 0, 0.1], (6, 1))            # pure z-displacement
    evg, c2l, cen, atoms, idx = _build(coords, disp_z)
    ca_centers = [1, 2, 3]
    r_ca = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix",
                              axis_centers=ca_centers)
    r_all = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    # with CA axis ~ z, a z-displacement is almost fully axial
    assert r_ca["axial_amplitude"] > r_ca["lateral_amplitude"]
    # with the contaminated all-atom axis (~x), the same field looks lateral
    assert r_all["axial_amplitude"] < r_all["lateral_amplitude"]


def test_zero_and_singleton_safe():
    # empty centers -> all-zero dict, no crash
    z = analyze_sse_element(np.zeros((0, 3)), {}, [], [], {}, "helix")
    assert z["amplitude_mean"] == 0.0 and z["tilting_angle"] == 0.0
    # single atom -> finite, no crash
    evg, c2l, cen, atoms, idx = _build([[0, 0, 0.0]], [[0.1, 0, 0]])
    r = analyze_sse_element(evg, c2l, cen, atoms, idx, "helix")
    assert np.isfinite(r["com_amplitude"])


def test_sse_umap_feature_names_exist():
    """Regression guard for the v1.0.5 bending bug: every metric used by the
    SSE-UMAP must be a key actually produced by analyze_sse_element, and the
    metric list must come from the single source of truth in core."""
    from modenanalyse_2fe2s.core import SSE_UMAP_METRICS

    evg, c2l, cen, atoms, idx = _build(Z, np.cross([0.05, 0, 0], Z))
    produced = set(analyze_sse_element(evg, c2l, cen, atoms, idx, "helix").keys())

    missing = [m for m in SSE_UMAP_METRICS if m not in produced]
    assert not missing, f"SSE-UMAP references non-existent descriptor(s): {missing}"
    # the buggy names must never reappear as actual metrics
    assert "bending_std" not in SSE_UMAP_METRICS
    assert "bending_mean" not in SSE_UMAP_METRICS

    # the consumer must use the shared constant, not a local literal
    src = (SRC / "embedding.py").read_text(encoding="utf-8")
    assert "_SSE_UMAP_METRICS" in src, \
        "embedding.py should consume core.SSE_UMAP_METRICS (single source of truth)"
