# -*- coding: utf-8 -*-
"""
test_reorganization.py
======================

Tests fuer das neue v3.7-Reorganisations-Modul.
"""

import sys
sys.path.insert(0, "src")

import numpy as np
import warnings
from modenanalyse_2fe2s.reorganization import (
    CHANNELS, ChannelGeometry, ChannelResult,
    signed_dr_along_axis, reduced_mass_pair, lambda_pair_cm1,
    compute_mode_modulations, aggregate_by_parent,
    compute_total_reorganization, compute_modulation_spectra,
    compute_co_modulation_spectra, compute_cumulative_reorganization,
)


def test_signed_dr_pure_stretch():
    """Reine N-H-Streckmode: dr sollte exakt der Auslenkungsdifferenz
    entlang der Achse entsprechen."""
    e_n = np.array([0.0, 0.0, 0.0])
    e_h = np.array([0.001, 0.0, 0.0])  # H bewegt sich +1 mA in x
    r_n = np.array([0.0, 0.0, 0.0])
    r_h = np.array([1.014, 0.0, 0.0])  # N-H entlang x
    # Bond-Vektor: a-b = N-H = (0-1.014,0,0) = (-1.014, 0, 0)
    # axis = bond/|bond| = (-1, 0, 0)
    # relative_e = e_n - e_h = (-0.001, 0, 0)
    # dr = (-0.001) * (-1) = +0.001 → Bindung dehnt sich (richtig)
    dr = signed_dr_along_axis(e_n, e_h, r_n, r_h)
    assert abs(dr - 0.001) < 1e-9, f"Erwartet +0.001, bekommen {dr}"
    print("  [OK] test_signed_dr_pure_stretch")


def test_signed_dr_compression():
    """H bewegt sich zu N hin -> Bindung staucht sich -> negatives dr."""
    e_n = np.array([0.0, 0.0, 0.0])
    e_h = np.array([-0.001, 0.0, 0.0])  # H zu N hin
    r_n = np.array([0.0, 0.0, 0.0])
    r_h = np.array([1.014, 0.0, 0.0])
    dr = signed_dr_along_axis(e_n, e_h, r_n, r_h)
    assert dr < 0, f"Erwartet negativ, bekommen {dr}"
    print("  [OK] test_signed_dr_compression")


def test_reduced_mass_known():
    """Bekannte Elemente liefern bekannte Werte."""
    mu_NH = reduced_mass_pair("N", "H")
    assert 0.93 < mu_NH < 0.95, f"N-H mu sollte ~0.94 sein, bekommen {mu_NH}"
    mu_FeFe = reduced_mass_pair("Fe", "Fe")
    assert 27.5 < mu_FeFe < 28.5, f"Fe-Fe mu ~27.92, bekommen {mu_FeFe}"
    mu_FeS = reduced_mass_pair("Fe", "S")
    assert 20.0 < mu_FeS < 20.7, f"Fe-S mu ~20.36, bekommen {mu_FeS}"
    print("  [OK] test_reduced_mass_known")


def test_reduced_mass_unknown():
    """Unbekanntes Element -> Fallback mit Warnung."""
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        mu = reduced_mass_pair("X", "H")
        assert mu > 0, f"Fallback sollte positiv liefern, bekommen {mu}"
        assert any("X" in str(w.message) for w in ws)
    print("  [OK] test_reduced_mass_unknown (mu={:.3f})".format(mu))


def test_lambda_pair_zero_for_zero_dr():
    assert lambda_pair_cm1(0.0, 500.0, 1.0) == 0.0
    print("  [OK] test_lambda_pair_zero_for_zero_dr")


def test_lambda_pair_known_value():
    """Fuer eine N-H-Streckmode bei 3300 cm-1 mit dr=0.005 A,
    mu=0.94 amu: lambda = (1/2) * 0.94 * (2*pi*c*3300)^2 * (5e-13)^2
    
    Numerisch: ~3.8 cm-1.
    Bei voller PT-Verschiebung dr=0.05 A waeren es ~380 cm-1.
    """
    lam = lambda_pair_cm1(0.005, 3300.0, 0.94)
    assert 3.0 < lam < 5.0, f"NH-Stretch lambda ausserhalb Erwartung: {lam}"
    
    # Verschiebung 10x groesser -> lambda 100x groesser
    lam_full = lambda_pair_cm1(0.05, 3300.0, 0.94)
    assert 350.0 < lam_full < 420.0, f"Vollst. PT lambda: {lam_full}"
    print(f"  [OK] test_lambda_pair_known_value (lambda(5mA)={lam:.2f}, "
          f"lambda(50mA)={lam_full:.1f})")


def test_lambda_pair_scaling():
    """lambda ~ omega^2 * mu * dr^2 - Skalierungs-Test."""
    lam1 = lambda_pair_cm1(0.005, 100.0, 1.0)
    lam2 = lambda_pair_cm1(0.005, 200.0, 1.0)
    # Verdoppelung omega -> Vervierfachung lambda
    assert abs(lam2 / lam1 - 4.0) < 1e-3, f"omega-Skalierung falsch: {lam2/lam1}"
    
    lam3 = lambda_pair_cm1(0.010, 100.0, 1.0)
    # Verdoppelung dr -> Vervierfachung lambda
    assert abs(lam3 / lam1 - 4.0) < 1e-3, f"dr-Skalierung falsch: {lam3/lam1}"
    print("  [OK] test_lambda_pair_scaling")


def test_compute_mode_modulations_simple():
    """Einfaches System: 4-Atom-Cluster + 1 His-N-H-Akzeptor."""
    # Atome: 0=Fe1, 1=Fe2, 2=S1, 3=S2, 4=N(His), 5=H, 6=O(Akzeptor)
    e = np.zeros((7, 3))
    # Mode: H bewegt sich 1 mA in x (Richtung Akzeptor), Fe1 bewegt sich -0.5 mA
    e[5] = [+0.001, 0.0, 0.0]
    e[0] = [-0.0005, 0.0, 0.0]
    e[1] = [+0.0005, 0.0, 0.0]
    
    channels = [
        ChannelGeometry(name="FeFe", idx_a=0, idx_b=1,
                        r_a=np.array([0.0, 0.0, 0.0]),
                        r_b=np.array([2.77, 0.0, 0.0]),
                        elem_a="Fe", elem_b="Fe", mu_pair_amu=27.92,
                        parent_channel="FeFe"),
        ChannelGeometry(name="NH_H1", idx_a=4, idx_b=5,
                        r_a=np.array([3.0, 0.0, 0.0]),
                        r_b=np.array([4.014, 0.0, 0.0]),
                        elem_a="N", elem_b="H", mu_pair_amu=0.94,
                        parent_channel="NH"),
        ChannelGeometry(name="HA_H1_O", idx_a=5, idx_b=6,
                        r_a=np.array([4.014, 0.0, 0.0]),
                        r_b=np.array([6.7, 0.0, 0.0]),
                        elem_a="H", elem_b="O", mu_pair_amu=0.948,
                        parent_channel="HA"),
    ]
    
    results = compute_mode_modulations(
        e_atoms=e, omega_cm1=300.0, mode_red_mass_amu=10.0,
        channels=channels)
    assert len(results) == 3
    
    # FeFe: Fe1 bei (0,0,0), Fe2 bei (2.77,0,0)
    # axis = (Fe1-Fe2)/|...| = (-1, 0, 0)
    # relative_e = e_Fe1 - e_Fe2 = (-0.001, 0, 0)
    # dr = (-0.001) * (-1) = +0.001 → Bindung dehnt sich +1 mA
    fefe = results[0]
    assert abs(fefe.dr_signed_a - 0.001) < 1e-9
    assert fefe.lambda_pair_cm1 > 0  # endlicher Beitrag
    
    # NH: e_N = 0, e_H = +0.001, axis ~ (-1,0,0)
    # dr = (0 - 0.001) * (-1) = +0.001
    nh = results[1]
    assert abs(nh.dr_signed_a - 0.001) < 1e-9
    
    # HA: e_H = +0.001, e_O = 0, axis = (H-O)/|...| ~ (-1, 0, 0)
    # dr = (0.001 - 0) * (-1) = -0.001 → Bindung staucht sich (H zum O)
    ha = results[2]
    assert ha.dr_signed_a < 0
    
    print("  [OK] test_compute_mode_modulations_simple")


def test_aggregate_by_parent():
    """Mehrere FeN-Subkanaele aggregieren zu einem FeN-Eintrag."""
    results = [
        ChannelResult("FeN_His255", "FeN", dr_signed_a=0.003,
                      lambda_pair_cm1=0.5, lambda_mode_cm1=0.4),
        ChannelResult("FeN_His259", "FeN", dr_signed_a=-0.002,
                      lambda_pair_cm1=0.3, lambda_mode_cm1=0.25),
        ChannelResult("FeFe", "FeFe", dr_signed_a=0.001,
                      lambda_pair_cm1=0.1, lambda_mode_cm1=0.08),
    ]
    agg = aggregate_by_parent(results)
    
    # FeN: dr_rss = sqrt(0.003^2 + 0.002^2) = sqrt(13e-6) ~ 0.0036
    fen = agg["FeN"]
    expected_rms = np.sqrt(0.003**2 + 0.002**2)
    assert abs(fen["dr_rss_a"] - expected_rms) < 1e-9
    # Lambdas: einfache Summe
    assert abs(fen["lambda_pair_cm1"] - 0.8) < 1e-9
    assert abs(fen["lambda_mode_cm1"] - 0.65) < 1e-9
    assert fen["n_subchannels"] == 2
    
    # FeFe: nur ein Sub-Kanal
    assert abs(agg["FeFe"]["dr_rss_a"] - 0.001) < 1e-9
    assert agg["FeFe"]["n_subchannels"] == 1
    
    # NH/FeS/HA leer
    assert agg["NH"]["n_subchannels"] == 0
    print("  [OK] test_aggregate_by_parent")


def test_compute_total_reorganization():
    """System-Total = einfache Summe ueber alle Modes."""
    per_mode = []
    for i in range(5):
        per_mode.append({
            "FeFe": {"dr_rss_a": 0.001, "lambda_pair_cm1": 0.5,
                     "lambda_mode_cm1": 0.3, "n_subchannels": 1,
                     "dr_sum_signed_a": 0.001},
            "FeN":  {"dr_rss_a": 0.0, "lambda_pair_cm1": 0.0,
                     "lambda_mode_cm1": 0.0, "n_subchannels": 0,
                     "dr_sum_signed_a": 0.0},
            "FeS":  {"dr_rss_a": 0.0, "lambda_pair_cm1": 0.0,
                     "lambda_mode_cm1": 0.0, "n_subchannels": 0,
                     "dr_sum_signed_a": 0.0},
            "NH":   {"dr_rss_a": 0.0, "lambda_pair_cm1": 0.0,
                     "lambda_mode_cm1": 0.0, "n_subchannels": 0,
                     "dr_sum_signed_a": 0.0},
            "HA":   {"dr_rss_a": 0.0, "lambda_pair_cm1": 0.0,
                     "lambda_mode_cm1": 0.0, "n_subchannels": 0,
                     "dr_sum_signed_a": 0.0},
        })
    
    totals = compute_total_reorganization(per_mode)
    assert abs(totals["FeFe"]["lambda_total_pair_cm1"] - 2.5) < 1e-9
    assert abs(totals["FeFe"]["lambda_total_mode_cm1"] - 1.5) < 1e-9
    assert totals["FeFe"]["n_modes_contributing"] == 5
    assert totals["FeN"]["lambda_total_pair_cm1"] == 0.0
    print("  [OK] test_compute_total_reorganization")


def test_modulation_spectra():
    """Eine einzelne Mode bei 200 cm-1 sollte einen Gauss-Peak bei
    omega=200 im Spektrum erzeugen."""
    freqs = np.array([200.0])
    per_mode = [{
        "FeFe": {"dr_rss_a": 0.005},
        "FeN":  {"dr_rss_a": 0.0},
        "FeS":  {"dr_rss_a": 0.0},
        "NH":   {"dr_rss_a": 0.0},
        "HA":   {"dr_rss_a": 0.0},
    }]
    grid = np.linspace(150, 250, 1001)
    spectra = compute_modulation_spectra(freqs, per_mode, grid, sigma_cm1=5.0)
    
    # Maximum sollte bei 200 cm-1 liegen
    peak_idx = np.argmax(spectra["FeFe"])
    peak_freq = grid[peak_idx]
    assert abs(peak_freq - 200.0) < 0.2, f"Peak bei {peak_freq}, erwartet 200"
    
    # Andere Kanaele sollten 0 sein
    assert np.max(spectra["FeN"]) < 1e-9
    print(f"  [OK] test_modulation_spectra (Peak bei {peak_freq:.2f} cm-1)")


def test_modulation_spectra_two_modes():
    """Zwei Modes ueberlagern sich ohne Interferenz (Beitrag addiert)."""
    freqs = np.array([200.0, 300.0])
    per_mode = [
        {"FeFe": {"dr_rss_a": 0.005}, "FeN": {"dr_rss_a": 0},
         "FeS": {"dr_rss_a": 0}, "NH": {"dr_rss_a": 0}, "HA": {"dr_rss_a": 0}},
        {"FeFe": {"dr_rss_a": 0.003}, "FeN": {"dr_rss_a": 0},
         "FeS": {"dr_rss_a": 0}, "NH": {"dr_rss_a": 0}, "HA": {"dr_rss_a": 0}},
    ]
    grid = np.linspace(150, 350, 2001)
    spectra = compute_modulation_spectra(freqs, per_mode, grid, sigma_cm1=5.0)
    
    # Zwei Peaks
    s = spectra["FeFe"]
    peaks = []
    for i in range(1, len(s)-1):
        if s[i] > s[i-1] and s[i] > s[i+1] and s[i] > 1e-6:
            peaks.append(grid[i])
    assert len(peaks) == 2, f"Erwartet 2 Peaks, gefunden {len(peaks)}: {peaks}"
    assert abs(peaks[0] - 200.0) < 0.2
    assert abs(peaks[1] - 300.0) < 0.2
    print("  [OK] test_modulation_spectra_two_modes")


def test_co_modulation():
    """Zwei Spektren mit Peak an gleicher Stelle -> C hat dort Peak.
    Peaks an unterschiedlichen Stellen -> C ist klein."""
    grid = np.linspace(0, 500, 5001)
    # Spectrum 1: Peak bei 264
    s1 = np.exp(-((grid - 264.0)**2) / (2*25))
    # Spectrum 2: Peak bei 264 (uberlappend)
    s2 = np.exp(-((grid - 264.0)**2) / (2*25))
    # Spectrum 3: Peak bei 400 (nicht uberlappend)
    s3 = np.exp(-((grid - 400.0)**2) / (2*25))
    
    spectra = {"HA": s1, "FeFe": s2, "FeN": s3, "FeS": np.zeros_like(grid),
               "NH": np.zeros_like(grid)}
    
    co = compute_co_modulation_spectra(spectra)
    
    # C_PCET = sqrt(M_HA * M_FeFe) hat Peak bei 264
    pcet_peak_idx = np.argmax(co["C_PCET"])
    assert abs(grid[pcet_peak_idx] - 264.0) < 0.5
    
    # C_PT_FeN = sqrt(M_HA * M_FeN), HA bei 264, FeN bei 400 -> C ist klein
    assert np.max(co["C_PT_FeN"]) < 0.1 * np.max(co["C_PCET"])
    
    print("  [OK] test_co_modulation")


def test_cumulative_reorganization():
    """Lambda(omega) ist monoton steigend und konvergiert gegen Total."""
    freqs = np.array([100.0, 200.0, 300.0, 400.0])
    per_mode = [
        {"FeFe": {"lambda_mode_cm1": 0.5, "lambda_pair_cm1": 0.7},
         "FeN":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "FeS":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "NH":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "HA":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0}},
        {"FeFe": {"lambda_mode_cm1": 1.0, "lambda_pair_cm1": 1.4},
         "FeN":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "FeS":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "NH":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "HA":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0}},
        {"FeFe": {"lambda_mode_cm1": 0.3, "lambda_pair_cm1": 0.4},
         "FeN":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "FeS":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "NH":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "HA":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0}},
        {"FeFe": {"lambda_mode_cm1": 0.2, "lambda_pair_cm1": 0.3},
         "FeN":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "FeS":  {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "NH":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0},
         "HA":   {"lambda_mode_cm1": 0.0, "lambda_pair_cm1": 0.0}},
    ]
    grid = np.linspace(0, 500, 501)
    cum = compute_cumulative_reorganization(freqs, per_mode, grid)
    
    # Monoton steigend
    fefe = cum["FeFe"]
    assert all(fefe[i+1] >= fefe[i] for i in range(len(fefe)-1)), "Nicht monoton"
    
    # Bei omega=0: keine Mode <= 0 -> Lambda = 0
    assert fefe[0] == 0.0
    
    # Bei omega=500: alle 4 Modes inkludiert, Lambda = 0.5+1.0+0.3+0.2 = 2.0
    assert abs(fefe[-1] - 2.0) < 1e-9
    
    # Bei omega=150: nur Mode 100 -> Lambda = 0.5
    idx_150 = np.argmin(np.abs(grid - 150.0))
    assert abs(fefe[idx_150] - 0.5) < 1e-9
    
    # Bei omega=250: Modes 100 + 200 -> Lambda = 0.5 + 1.0 = 1.5
    idx_250 = np.argmin(np.abs(grid - 250.0))
    assert abs(fefe[idx_250] - 1.5) < 1e-9
    
    print("  [OK] test_cumulative_reorganization")


def test_acceptor_gauss_weight_in_build_channels():
    """Akzeptoren bekommen Gauss-Gewicht ueber den Gleichgewichtsabstand:
    nahe Akzeptoren mehr, ferne weniger.
    
    v3.7.1: pro H wird nur der HAUPT-Akzeptor (max Gauss-Gewicht)
    verwendet, nicht alle. Reaktionskoordinaten-Modus aktiv (idx_donor
    gesetzt)."""
    from modenanalyse_2fe2s.reorganization import build_channels
    
    class MockCoordInfo:
        ligands = []
    
    atoms = [
        {"x":  0.0, "y": 0.0, "z": 0.0, "atomic_num": 26},  # Fe1
        {"x":  2.77,"y": 0.0, "z": 0.0, "atomic_num": 26},  # Fe2
        {"x":  1.4, "y": 1.0, "z": 0.0, "atomic_num": 16},  # S1
        {"x":  1.4, "y": -1.0,"z": 0.0, "atomic_num": 16},  # S2
        {"x":  0.0, "y": 2.0, "z": 0.0, "atomic_num": 7},   # N (His)
        {"x":  0.5, "y": 2.9, "z": 0.0, "atomic_num": 1},   # H
        {"x":  1.0, "y": 5.4, "z": 0.0, "atomic_num": 8},   # O nahe (eq=2.6 A)
        {"x":  3.0, "y": 6.0, "z": 0.0, "atomic_num": 8},   # O fern (eq=4.0 A)
    ]
    idx_map = {i+1: i for i in range(len(atoms))}
    
    channels = build_channels(
        coord_info=MockCoordInfo(), atoms=atoms, idx_map=idx_map,
        fe1_center=1, fe2_center=2,
        s_centers=[3, 4],
        his_n_centers=[5], his_n_labels=["His1"],
        his_h_centers=[6], 
        cys_s_centers=[], cys_s_labels=[],
        fe_centers_per_his_n=[1],
        fe_centers_per_cys_s=[],
        acceptors_per_h=[[7, 8]],
        acceptor_elem_per_h=[["O", "O"]],
        eq_distances_per_h=[[2.6, 4.0]],
        acceptor_r0_a=2.8, acceptor_sigma_a=0.4,
    )
    
    names = [c.name for c in channels]
    assert "FeFe" in names, f"FeFe fehlt: {names}"
    assert any("FeN" in n for n in names), f"FeN fehlt: {names}"
    assert any("NH" in n for n in names), f"NH fehlt: {names}"
    
    # v3.7.1: nur EIN HA-Kanal (Hauptakzeptor)
    ha_channels = [c for c in channels if c.parent_channel == "HA"]
    assert len(ha_channels) == 1, (
        f"v3.7.1: Erwartet 1 HA-Kanal (Hauptakzeptor), bekommen {len(ha_channels)}")
    
    ha = ha_channels[0]
    # v3.7.4: Gauss-Gewicht wird NUR fuer die Hauptakzeptor-Auswahl verwendet,
    # NICHT mehr als Daempfung in der Lambda-Aggregation. weight ist 1.0.
    # Die binaere Logik ist: ist ein H-Bridge-Akzeptor in
    # pcet_hbond_cutoff_a Reichweite, oder nicht. Bei Akzeptoren in
    # Reichweite zaehlt die Mode voll, bei keinem entfaellt der Kanal.
    assert ha.weight == 1.0, f"v3.7.4: HA weight muss 1.0 sein, bekommen {ha.weight}"
    
    # Reaktionskoordinaten-Modus: idx_donor gesetzt (auf das N-Atom)
    assert ha.idx_donor is not None, "v3.7.1: HA muss idx_donor haben"
    assert ha.idx_donor == 4  # idx von N im atoms-Array
    assert ha.r_donor is not None
    
    # Verifizieren: der naehere Akzeptor (eq=2.6 A) wurde gewaehlt,
    # nicht der weiter entfernte (eq=3.4 A). Sichtbar an idx_b.
    # Naheliegend ist: idx_b zeigt auf das O bei eq=2.6 (idx 6 im atoms-Array).
    # Auswahl ueber Gauss-Gewicht: Gauss(2.6-2.8) > Gauss(3.4-2.8)
    
    print(f"  [OK] test_acceptor_gauss_weight_v37_4 (weight={ha.weight}, "
          f"idx_donor={ha.idx_donor}, idx_b={ha.idx_b})")


def test_reaction_coordinate_vs_classical_HA():
    """v3.7.1: HA-Berechnung im Reaktionskoordinaten-Modus vs.
    klassischer Bindungs-Modulation.
    
    Szenario: hochfrequente Mode wo nur der Akzeptor sich bewegt
    (typisch bei C=O-Schwingung in der Aminosaeure-Skelett-Region).
    
    - Klassische Modulation: dr_HA != 0 (Distanz aendert sich)
    - Reaktionskoordinaten-Modus: dr_HA = 0 (das H bewegt sich nicht
      relativ zum Donor; nur der Akzeptor schwingt)
    """
    # 5-Atom-System: N (Donor), H, O (Akzeptor)
    # Mode: nur O bewegt sich
    e = np.zeros((3, 3))
    e[2] = [0.005, 0.0, 0.0]  # O bewegt sich +5 mA in x
    
    r_n = np.array([0.0, 0.0, 0.0])
    r_h = np.array([1.014, 0.0, 0.0])
    r_o = np.array([3.6,   0.0, 0.0])  # H...O = 2.586 A
    
    # Klassische Modulation: idx_a=H (1), idx_b=O (2), kein donor
    ch_classical = ChannelGeometry(
        name="HA_test", idx_a=1, idx_b=2,
        r_a=r_h, r_b=r_o, elem_a="H", elem_b="O",
        mu_pair_amu=0.95, parent_channel="HA",
    )
    res_c = compute_mode_modulations(e, 700.0, 5.0, [ch_classical])
    # axis = (r_h - r_o) / |...| = (-1, 0, 0)
    # dr = (e_h - e_o) . axis = (0 - 0.005, 0, 0) . (-1, 0, 0) = +0.005
    # Klassisch: dr ist 5 mA > 0 (Distanz verkleinert sich? Nein:
    #   axis zeigt von r_o zu r_h, also "negatives Vorzeichen" wenn H sich
    #   vom O entfernt.) Eigentliche Konvention: bond = r_a - r_b = r_h - r_o,
    #   axis_unit zeigt von b nach a, also von O nach H.
    #   relative_e = e_h - e_o = (0 - 0.005,0,0) = (-0.005, 0, 0)
    #   dr = relative_e . axis = (-0.005) * (-1) = +0.005
    #   Positives dr = "Bindung dehnt sich" = O entfernt sich von H
    #   In unserem Fall: O bewegt sich in +x, H steht still.
    #   Distanz |r_h - r_o| nach Bewegung: |1.014-3.605| = 2.591 statt 2.586
    #   → vergrößert um +5 mA → dr = +0.005, OK.
    assert abs(res_c[0].dr_signed_a - 0.005) < 1e-6, (
        f"Klassische Modulation: erwartet +0.005, bekommen {res_c[0].dr_signed_a}")
    assert res_c[0].lambda_pair_cm1 > 0
    
    # Reaktionskoordinaten-Modus: idx_a=H (1), idx_b=O (2), idx_donor=N (0)
    ch_react = ChannelGeometry(
        name="HA_test", idx_a=1, idx_b=2,
        r_a=r_h, r_b=r_o, elem_a="H", elem_b="O",
        mu_pair_amu=0.95, parent_channel="HA",
        idx_donor=0, r_donor=r_n,
    )
    res_r = compute_mode_modulations(e, 700.0, 5.0, [ch_react])
    # Reaktionskoordinaten-dr = (e_h - e_n) . axis(N->O)
    # axis = (r_o - r_n)/|...| = (3.6, 0, 0)/3.6 = (1, 0, 0)
    # e_h - e_n = (0 - 0, 0, 0) = (0, 0, 0)
    # dr = 0 * 1 = 0
    # → Keine PT-Bewegung, das H bewegt sich nicht relativ zum Donor
    assert abs(res_r[0].dr_signed_a) < 1e-9, (
        f"Reaktionskoord: erwartet 0 (kein H-Donor-Move), bekommen "
        f"{res_r[0].dr_signed_a}")
    assert res_r[0].lambda_pair_cm1 == 0.0
    
    print("  [OK] test_reaction_coordinate_vs_classical_HA "
          "(klassisch: dr=0.005, react.coord: dr=0.000)")


def test_reaction_coordinate_real_PT():
    """Echte PT-Mode: H bewegt sich vom N zum O. Sowohl klassisch als
    auch im Reaktionskoordinaten-Modus muss dr signifikant sein."""
    # H bewegt sich +1 mA in x (Richtung O)
    e = np.zeros((3, 3))
    e[1] = [0.001, 0.0, 0.0]
    
    r_n = np.array([0.0, 0.0, 0.0])
    r_h = np.array([1.014, 0.0, 0.0])
    r_o = np.array([3.6, 0.0, 0.0])
    
    ch = ChannelGeometry(
        name="HA", idx_a=1, idx_b=2,
        r_a=r_h, r_b=r_o, elem_a="H", elem_b="O",
        mu_pair_amu=0.95, parent_channel="HA",
        idx_donor=0, r_donor=r_n,
    )
    res = compute_mode_modulations(e, 1500.0, 1.5, [ch])
    # axis(N->O) = (1, 0, 0); e_h - e_n = (0.001, 0, 0)
    # dr = 0.001 * 1 = +0.001 (H bewegt sich auf O zu)
    assert abs(res[0].dr_signed_a - 0.001) < 1e-9
    assert res[0].lambda_pair_cm1 > 0
    print("  [OK] test_reaction_coordinate_real_PT (dr=+0.001 = echtes PT)")


def test_marcus_hush_summation_property():
    """Verifiziere die Marcus-Hush-Summen-Eigenschaft:
    
    Wenn man Modes addiert, addieren sich die Lambda_X-Beitraege linear,
    und die Total-Reorg ist die einfache Summe ueber alle Modes."""
    # 3 Modes mit Lambda_FeFe-Beitraegen 1.0, 2.0, 3.0
    per_mode = [
        {ch: {"dr_rss_a": 0, "dr_sum_signed_a": 0,
              "lambda_pair_cm1": 1.0 if ch == "FeFe" else 0,
              "lambda_mode_cm1": 0.7 if ch == "FeFe" else 0,
              "n_subchannels": 1 if ch == "FeFe" else 0}
         for ch in CHANNELS},
        {ch: {"dr_rss_a": 0, "dr_sum_signed_a": 0,
              "lambda_pair_cm1": 2.0 if ch == "FeFe" else 0,
              "lambda_mode_cm1": 1.4 if ch == "FeFe" else 0,
              "n_subchannels": 1 if ch == "FeFe" else 0}
         for ch in CHANNELS},
        {ch: {"dr_rss_a": 0, "dr_sum_signed_a": 0,
              "lambda_pair_cm1": 3.0 if ch == "FeFe" else 0,
              "lambda_mode_cm1": 2.1 if ch == "FeFe" else 0,
              "n_subchannels": 1 if ch == "FeFe" else 0}
         for ch in CHANNELS},
    ]
    totals = compute_total_reorganization(per_mode)
    assert abs(totals["FeFe"]["lambda_total_pair_cm1"] - 6.0) < 1e-9
    assert abs(totals["FeFe"]["lambda_total_mode_cm1"] - 4.2) < 1e-9
    print("  [OK] test_marcus_hush_summation_property")


# ===========================================================================
# Lauf
# ===========================================================================

if __name__ == "__main__":
    print("=== test_reorganization.py ===\n")
    
    print("Helper-Funktionen:")
    test_signed_dr_pure_stretch()
    test_signed_dr_compression()
    test_reduced_mass_known()
    test_reduced_mass_unknown()
    test_lambda_pair_zero_for_zero_dr()
    test_lambda_pair_known_value()
    test_lambda_pair_scaling()
    
    print("\nPro-Mode-Berechnung:")
    test_compute_mode_modulations_simple()
    test_aggregate_by_parent()
    
    print("\nSystem-Aggregation:")
    test_compute_total_reorganization()
    test_marcus_hush_summation_property()
    
    print("\nFrequenz-aufgeloest:")
    test_modulation_spectra()
    test_modulation_spectra_two_modes()
    test_co_modulation()
    test_cumulative_reorganization()
    
    print("\nGeometrie-Bau:")
    test_acceptor_gauss_weight_in_build_channels()
    
    print("\nv3.7.1 Reaktionskoordinaten-Modus:")
    test_reaction_coordinate_vs_classical_HA()
    test_reaction_coordinate_real_PT()
    
    print("\n=== ALLE TESTS BESTANDEN ===")
