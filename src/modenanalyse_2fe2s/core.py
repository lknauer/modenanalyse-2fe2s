# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

# -*- coding: utf-8 -*-
"""
core.py
====================
Physikalische Berechnungen per normal mode.

Oeffentliche functions
-----------------------
compute_thermal_amplitude
    QM-RMS-Amplitude a harmonischen Oszillators.
analyze_mode
    Vollstaendige Analyse a einzelnen normal mode.
analyze_fe_ligand
    Fe-N/S/O Streck- and Biegebewegungen.
analyze_his_hn
    H-N bondaenderung for protonierten histidineen.
analyze_all_ss
    Sekundaerstruktur-Amplitudenanalyse aller SS-elements.
classify_kernel_mode_from_evg
    Klassifikation after D2h-Symmetriekoordinaten.
compute_scsd_for_mode_full
    SCSD-Zerlegung the clustergeometrie.

Bugfixes (gegenvia Vorversion)
---------------------------------
B1  _run_scsd:                    Versucht mehrere scsdpy-API-Formate.
B2  analyze_all_ss:               Korrekte Indizierung via ``c2l``.
B3  compute_thermal_amplitude:    Fehlender Faktor 2 in the Nenner behoben.
B7  analyze_fe_ligand:            Fe-Center dynamisch from ``LigandInfo``.
"""
from __future__ import annotations
import math, warnings
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .config import Config
from .logio     import RunLog, get_eigvec, BlockInfo, _HIS


def _get_eigvec_smart(filepath: str,
                       bi: BlockInfo,
                       col: int,
                       include_hydrogen: bool,
                       cfg=None,
                       atoms=None,
                       idx_map=None):
    """Format-agnostischer Eigenvektor-Reader.

    For Gaussian is the bestehende ``get_eigvec`` called
    (liest Bytes from the .log-File). For ORCA is the in
    ``cfg._orca_parse_result`` injizierte ``ParseResult`` used --
    the completee Eigenvector liegt schon in the RAM.

    Returns
    -------
    centers : list of int
        Atom-Center-Nummern (1-basiert).
    evg : ndarray of shape (n_centers, 3)
        Cartesian-unit-Eigenvector (\\textbf{ohne} u_rms-Skalierung).
    """
    pr = getattr(cfg, "_orca_parse_result", None) if cfg is not None else None
    if pr is None:
        # Gaussian-path
        return get_eigvec(filepath, bi, col, include_hydrogen=include_hydrogen)

    # ORCA-path: from the ParseResult lesen
    from .orca_io import parseresult_to_atoms, get_eigvec_orca
    # mode_num from the Block (col-basiert)
    mode_num = bi.mode_nums[col] if col < len(bi.mode_nums) else (col + 1)
    # Atome in the richtigen Format (mit/ohne H)
    pr_atoms, _pr_idx = parseresult_to_atoms(pr, include_hydrogen=include_hydrogen)
    # get_eigvec_orca gibt (centers, evg)-Tuple zurueck; wir nehmen
    # only the evg-Teil and nutzen pr_atoms als Centers-Quelle
    _, evg = get_eigvec_orca(pr, mode_num, pr_atoms, _pr_idx,
                              include_hydrogen=include_hydrogen)
    centers = [a["center"] for a in pr_atoms]
    return centers, evg
from .geometry import CoordInfo, LigandInfo


MAX_REASONABLE_KERNEL_D  = 1.0

# Oeffentliche Konstante - is von Embedding and Export importiert
SCORE_KEYS: List[str] = [
    "Translation", "Umbrella", "Hinge-Fe", "Hinge-S", "OOP-Twist",
    "Fe-stretching", "S-stretching", "Breathing", "Rhombus-shear",
    "Rotation-ip",
]
_SCORE_KEYS = SCORE_KEYS   # Abwaertskompatibilitaet


# Set of group names for which the "empty eigenvector" warning has already
# been issued during the current run. Used by analyze_mode() to suppress
# duplicate warnings; populated as modes are processed and reset by
# reset_warning_state() at the start of each run.
_WARNED_EMPTY_GROUP: Set[str] = set()

# v1.0.4 new: identical mechanism, but for the Fe-ligand c2l lookup that can
# silently fall back to _zero_lig() (the v1.0.4 fix targets the same class of
# silent data-loss bug as the v1.0.2 group_map fix). One warning per
# (residue label, failure reason) pair per run.
_WARNED_EMPTY_LIG:   Set[str] = set()

# v1.0.4 new: identical mechanism, but for protonated-His H-N lookups that
# can silently skip a ligand in analyze_his_hn(). Distinguishes legitimate
# deprotonation (no warning) from a true lookup error (warning).
_WARNED_HN_SKIP:     Set[str] = set()


def reset_warning_state() -> None:
    """Reset per-run warning state.

    Called by runner.run() before iterating over modes. Without this reset,
    a long-lived Python session that calls run() repeatedly would only see
    warnings from the first run.
    """
    _WARNED_EMPTY_GROUP.clear()
    _WARNED_EMPTY_LIG.clear()
    _WARNED_HN_SKIP.clear()


# ===========================================================================
# Thermische Amplitude  (Bugfix B3: Faktor 2)
# ===========================================================================

def compute_thermal_amplitude(freq_cm1:      float,
                               red_mass_amu:  float,
                               temp_k:        Optional[float],
                               amplitude:     float = 1.0,
                               ) -> float:
    """Computes the QM-RMS-Amplitude a harmonischen Oszillators.

    Parameters
    ----------
    freq_cm1 : float
        Schwingungsfrequenz in cm⁻^1. Muss > 0 sein.
    red_mass_amu : float
        Reduced mass in atomaren Masseneinheiten.
    temp_k : float or None
        Temperature in Kelvin. ``None`` aktiviert klassischen Modus.
    amplitude : float, optional
        Fallback-Amplitude in the klassischen Modus (``temp_k`` is ``None``
        or ``<= 0``). default is ``1.0``.

    Returns
    -------
    float
        RMS-Amplitude in Angstrom.

    Notes
    -----
    Formel (quantenmechanisch):

    .. math::

        u_\\text{rms} = \\sqrt{\\frac{\\hbar}{2 m_r \\omega}
                       \\coth\\!\\left(\\frac{\\hbar\\omega}{2 k_B T}\\right)}

    Bugfix B3: Vorherige Version hatte ``hbar / (m * omega)`` statt
    ``hbar / (2 * m * omega)``, was the value um Faktor ``sqrt(2)``
    ueberschaetzte.
    """
    if temp_k is None or temp_k <= 0 or freq_cm1 <= 0 or red_mass_amu <= 0:
        return amplitude

    hbar  = 1.054571817e-34      # J*s
    k_b   = 1.380649e-23         # J/K
    c     = 2.99792458e10        # cm/s
    amu   = 1.66053906660e-27    # kg

    omega = 2 * math.pi * c * freq_cm1
    m     = red_mass_amu * amu
    x     = hbar * omega / (2 * k_b * temp_k)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        coth = 1.0 / np.tanh(x) if x < 500 else 1.0

    # B3-Fix: hbar / (2 * m * omega)  statt  hbar / (m * omega)
    u_sq  = hbar / (2.0 * m * omega) * coth
    return float(np.sqrt(max(u_sq, 0.0))) * 1e10   # → Angstrom


# ===========================================================================
# error propagation
# ===========================================================================

def _oop(evg: np.ndarray, n: np.ndarray) -> float:
    """Out-of-plane Bruchteil: P/S with P=Σ(eᵢ*n)^2, S=Σ|eᵢ|^2. returns 0 for S≈0."""
    t = float(np.sum(evg**2))
    return float(np.sum((evg @ n)**2)) / t if t > 1e-30 else 0.0

def _ang(evg: np.ndarray, n: np.ndarray, amp_thr: float = 0.001) -> float:
    """Mean angle (Grad) between displacement vectors and normal vector n.

    Nur Atome with |eᵢ| > amp_thr are einbezogen. returns 0 zurueck wenn
    no Atom the Schwelle ueberschreitet.
    """
    r = []
    for v in evg:
        m = np.linalg.norm(v)
        if m > amp_thr:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                r.append(np.degrees(np.arcsin(min(1.0, abs(float(v @ n)) / m))))
    return float(np.mean(r)) if r else 0.0

def _kd(evg: np.ndarray, _=None) -> float:
    """Mittlere displacement length over alle Atome: mean(|eᵢ|)."""
    return float(np.mean(np.linalg.norm(evg, axis=1)))

def _sigma_oop(evg: np.ndarray, n_hat: np.ndarray, sigma: float) -> float:
    """Analytische error propagation for the OOP-Bruchteil P/S.

    Derivation: OOP = P/S, P = Σ(eᵢ*n̂)^2, S = Σ|eᵢ|^2
    → sigma^2(P/S) = 4sigma^2*P*(S−P)/S^3  (exakt, identical with Zentraldifferenzen).
    """
    S = float(np.sum(evg**2))
    if S < 1e-60: return 0.
    P = float(np.sum((evg @ n_hat)**2))
    return 2.0 * sigma * math.sqrt(max(P * (S - P), 0.)) / (S ** 1.5)

def _sigma_kd(evg: np.ndarray, sigma: float) -> float:
    """Analytische error propagation for the mittlere displacement length.

    Derivation: ∂(mean|eᵢ|)/∂eᵢⱼ = eᵢⱼ/(n*|eᵢ|)  →  sigma(kd) = sigma/sqrtn.
    """
    return sigma / math.sqrt(max(evg.shape[0], 1))

def _sigma_ang(evg: np.ndarray, n_hat: np.ndarray,
               amp_thr: float, sigma: float) -> float:
    """Analytische error propagation for the mittleren displacement angle.

    Derivation: sigma(arcsin(|eᵢ*n̂|/|eᵢ|)) = sigma/|eᵢ|
    → sigma(Mean) = sigma/n_v * sqrt(Σᵢ 1/|eᵢ|^2)  [in Grad].
    """
    norms = np.linalg.norm(evg, axis=1)
    valid = norms > amp_thr
    if not np.any(valid): return 0.
    n_v = int(valid.sum())
    return float(np.degrees(sigma / n_v * math.sqrt(
        float(np.sum(1.0 / np.maximum(norms[valid], 1e-30) ** 2)))))


# ===========================================================================
# Universal significance-Klassifikation (v3.5)
# ===========================================================================
# Ein einziges Konvention for alle Groessen with Unsicherheit:
#   |X| <= 1*sigma     → "trivial"      (value not von 0 unterscheidbar)
#   1*sigma < |X| <= k*sigma → "significant"  (1σ < |X|, k = high-Schwelle)
#   |X| > k*sigma      → "high"         (k-Sigma-Niveau, default k=3)
#
# Schwellen are over cfg.significance_threshold_low (default 1.0) und
# cfg.significance_threshold_high (default 3.0) konfigurierbar.

def classify_significance(value:     float,
                          sigma:     float,
                          thr_low:   float = 1.0,
                          thr_high:  float = 3.0,
                          ) -> str:
    """Klassifiziere einen value relativ to seiner Unsicherheit.

    Parameters
    ----------
    value : float
        Zu klassifizierender value (Vorzeichen is ignoriert).
    sigma : float
        1-Sigma-Unsicherheit on ``value``.
    thr_low : float, default 1.0
        Threshold for "trivial → signifikant" in multiples of sigma.
    thr_high : float, default 3.0
        Threshold for "signifikant → hoch" in multiples of sigma.

    Returns
    -------
    {"trivial", "significant", "high"}
    """
    if not math.isfinite(value) or not math.isfinite(sigma):
        return "trivial"
    if sigma <= 0:
        # Sigma=0 → value is exakt; Klassifikation after value selbst
        return "high" if abs(value) > 0 else "trivial"
    ratio = abs(value) / sigma
    if ratio > thr_high:
        return "high"
    if ratio > thr_low:
        return "significant"
    return "trivial"


def classify_difference_significance(v1:      float,
                                     s1:      float,
                                     v2:      float,
                                     s2:      float,
                                     thr_low: float = 1.0,
                                     thr_high: float = 3.0,
                                     ) -> str:
    """Klassifiziere the Signifikanz a Differenz |v1 - v2|.

    Verwendet Standard-error propagation: sigma_diff = sqrt(s1^2 + s2^2).
    Geeignet for Vergleiche wie bend_oop vs. bend_inp oder
    lig_oop_pct vs. lig_inp_pct.
    """
    diff = v1 - v2
    sigma_diff = math.sqrt(max(s1 * s1 + s2 * s2, 0.))
    return classify_significance(diff, sigma_diff, thr_low, thr_high)


# ===========================================================================
# OOP/INP-Zerlegung per Atom-Set (Ring 1/2/3)
# ===========================================================================

def _oop_ring_metrics(evg_ring:  np.ndarray,
                      n_hat:     np.ndarray,
                      sigma_ev:  float,
                      ) -> Dict[str, float]:
    """Berechne OOP/INP-Fractione + mittlere displacement for einen Atom-Ring.

    Identische Definition wie bisher for cluster core (kern_oop, kern_d),
    jetzt als generische Helper-Funktion for beliebige Ring-Atomsets.

    Parameters
    ----------
    evg_ring : ndarray of shape (n_ring, 3)
        Eigenvektor-Zeilen the Ring-Atome (bereits scaled with u_rms).
    n_hat : ndarray of shape (3,)
        cluster normal (Einheitsvektor).
    sigma_ev : float
        Eigenvektor-Unsicherheit per Komponente (scaled with u_rms).

    Returns
    -------
    dict with keyn:
        oop_pct, inp_pct, d, sigma_oop_pct, sigma_d, n_atoms
    """
    n_ring = int(evg_ring.shape[0])
    if n_ring == 0:
        return {
            "oop_pct":       0.0,
            "inp_pct":       0.0,
            "d":             0.0,
            "sigma_oop_pct": 0.0,
            "sigma_d":       0.0,
            "n_atoms":       0,
        }
    oop_frac = _oop(evg_ring, n_hat)
    d_mean   = _kd(evg_ring)
    s_oop    = _sigma_oop(evg_ring, n_hat, sigma_ev)
    s_d      = _sigma_kd(evg_ring, sigma_ev)
    return {
        "oop_pct":       100.0 * oop_frac,
        "inp_pct":       100.0 * (1.0 - oop_frac),
        "d":             d_mean,
        "sigma_oop_pct": 100.0 * s_oop,
        "sigma_d":       s_d,
        "n_atoms":       n_ring,
    }


# ===========================================================================
# ligands-Bend-Aufteilung in OOP/INP (relativ to cluster plane)
# ===========================================================================

def _bend_split(d_lig:    np.ndarray,
                d_fe:     np.ndarray,
                bhat:     np.ndarray,
                n_hat:    np.ndarray,
                sigma_ev: float,
                ) -> Dict[str, float]:
    """Zerlege the Fe-ligands-Biegekomponente in INP- and OOP-Fraction.

    Die Fe-ligands-Relativbewegung is zerplaces in:

        rel        = e_lig - e_fe
        stretch    = rel * b̂                       (entlang bond)
        rel_perp   = rel - (rel * b̂) b̂            (senkrecht to bond)
        bend_oop   = |rel_perp * n̂|                (out-of-plane bzgl. Cluster)
        bend_inp   = |rel_perp - (rel_perp*n̂) n̂|   (in-plane bzgl. Cluster)
        |bend|^2   = bend_inp^2 + bend_oop^2       (Pythagoras, exakt)

    Beide Komponenten are Skalare (nicht-negative lengthn) in Angstroem.

    Die Prozente bend_inp_pct and bend_oop_pct are **Quadrat-Fractione**
    (energy fractions), konsistent with the Konvention for oop_pct und
    lig_oop_pct:

        bend_inp_pct = 100 * bend_inp^2 / |bend|^2
        bend_oop_pct = 100 * bend_oop^2 / |bend|^2

    Damit gilt bend_inp_pct + bend_oop_pct = 100 % exakt (Pythagoras).

    Sigma values for the lengths follow from Standard-error propagation:
    jede Komponente is a Linearkombination the eigenvector-Komponenten
    von e_lig and e_fe (jebecauses with Unsicherheit sigma_ev unkorreliert),
    therefore sigma(bend_*) = sigma_ev * sqrt(2). Fuer the Prozente werden
    no Sigmas computed (waeren for kleiner Gesamt-Bend-Amplitude
    numerisch instabil; the Signifikanz lesen Anwender from den
    lengthn-valuesn and ihren Sigmas).

    Parameters
    ----------
    d_lig, d_fe : ndarray of shape (3,)
        displacement vectors von Ligand and Fe (in A, with u_rms scaled).
    bhat : ndarray of shape (3,)
        Einheitsvektor Fe→ligand bondrichtung.
    n_hat : ndarray of shape (3,)
        cluster normaln-Einheitsvektor.
    sigma_ev : float
        Eigenvektor-Unsicherheit per Komponente (scaled).

    Returns
    -------
    dict with Feldern:
        stretch (A), bend_inp (A), bend_oop (A),
        bend_inp_pct (%), bend_oop_pct (%),
        sigma_stretch, sigma_bend_inp, sigma_bend_oop
    """
    rel       = d_lig - d_fe
    stretch_s = float(rel @ bhat)              # signed scalar
    rel_perp  = rel - stretch_s * bhat
    oop_s     = float(rel_perp @ n_hat)        # signed scalar (oop component)
    inp_vec   = rel_perp - oop_s * n_hat       # in-plane vector
    inp_len   = float(np.linalg.norm(inp_vec)) # always non-negative

    stretch   = abs(stretch_s)
    bend_oop  = abs(oop_s)
    bend_inp  = inp_len
    bend_tot2 = bend_inp * bend_inp + bend_oop * bend_oop

    # Prozent-values: Quadratsummen-Fractione (konsistent with anderen *_pct
    # in the Programm, z. B. lig_oop_pct). So gilt bend_inp_pct +
    # bend_oop_pct = 100 % (oder NaN for to kleiner Gesamt-Bend).
    # Nur sinnvoll if Gesamt-Bend deutlich groesser als
    # Eigenvektor-noise. If |bend| < 2*sigma_diff (= sigma der
    # Differenz e_lig - e_fe), Prozente on NaN setzen.
    sigma_diff = sigma_ev * math.sqrt(2.0)
    bend_tot   = math.sqrt(bend_tot2)
    if bend_tot < 2.0 * sigma_diff or bend_tot2 < 1e-24:
        inp_pct  = float("nan")
        oop_pct  = float("nan")
    else:
        inp_pct  = 100.0 * (bend_inp * bend_inp) / bend_tot2
        oop_pct  = 100.0 * (bend_oop * bend_oop) / bend_tot2

    return {
        "stretch":          stretch,
        "bend_inp":         bend_inp,
        "bend_oop":         bend_oop,
        "bend_inp_pct":     inp_pct,
        "bend_oop_pct":     oop_pct,
        "sigma_stretch":    sigma_diff,
        "sigma_bend_inp":   sigma_diff,
        "sigma_bend_oop":   sigma_diff,
    }



# ===========================================================================
# Fe-ligands-Analyse  (Bugfix B7: naechstes Fe dynamisch bestimmt)
# ===========================================================================

def analyze_fe_ligand(evg:      np.ndarray,
                      c2l:      Dict[int,int],
                      coord_info: CoordInfo,
                      atoms:    List[Dict],
                      idx_map:  Dict[int,int],
                      cfg:      Config,
                      n_hat:    np.ndarray = None,
                      u_rms:    float = 1.0,
                      ) -> Dict:
    """Computes Fe-N/S/O Streck- and Biegebewegungen for alle ligands.

    v3.5: Biege-Komponente is in INP- and OOP-Fractione bezueglich der
    cluster plane zerplaces. ``n_hat`` (cluster normal) is Pflichtparameter
    for the Aufteilung; without ihn is only the alte Skalar-bend-Groesse
    zurueckgiven (Abwaertskompatibilitaet).

    Parameters
    ----------
    evg : ndarray of shape (n_atoms, 3)
        Eigenvector the aktuellen Mode (thermisch scaled with u_rms).
    c2l : dict of {int: int}
        Mapping Gaussian-Center → Zeilenindex in ``evg``.
    coord_info : CoordInfo
        Koordinations-Information (enthaelt ``LigandInfo``-Objekte).
    atoms : list of dict
        Gaussian-atom list.
    idx_map : dict of {int: int}
        Center → Atom-Index.
    cfg : Config
        Benecessaryt ``sigma_eigvec``,
        ``significance_threshold_low`` (default 1.0),
        ``significance_threshold_high`` (default 3.0).
    n_hat : ndarray of shape (3,), optional
        cluster normaln-Einheitsvektor. If None, is only die
        skalare bend-Groesse zurueckgiven (Abwaertskompatibilitaet).
    u_rms : float
        Thermal scaling (to Sigma-Berechnung; evg is bereits
        scaled, sigma_ev muss correspond tod scaled werden).

    Returns
    -------
    dict of {str: dict}
        Pro ligands-Label:
            stretch, bend (Skalar), bend_inp, bend_oop, bend_inp_pct,
            bend_oop_pct, sigmas dazu, bend_significance, element, fe_idx.

    Notes
    -----
    Sigmas: jede Komponente (stretch, bend_inp, bend_oop) is eine
    Linearkombination zweier unabhaengiger Eigenvektor entries
    (e_lig and e_fe). Mit standard uncertainty ``sigma_ev = sigma_eigvec
    * u_rms`` per Komponente folgt sigma Komponente:
        sigma(component) = sigma_ev * sqrt(2)
    (ein Faktor sqrt(2) per Differenz, no additionallyer Faktor
    from the Projektion on einen *deterministischen* Einheitsvektor).
    """
    sigma_ev = float(cfg.sigma_eigvec) * float(u_rms)
    s_sb     = sigma_ev * math.sqrt(2.0)
    thr_low  = float(getattr(cfg, "significance_threshold_low",  1.0))
    thr_high = float(getattr(cfg, "significance_threshold_high", 3.0))

    result: Dict[str, Dict] = {}

    for lig in coord_info.ligands:
        label = lig.res_label

        # Fe-Eigenvektor
        fe_row = c2l.get(lig.fe_center)
        if fe_row is None or fe_row >= evg.shape[0]:
            # v1.0.4: previously silent — same bug class as the v1.0.2
            # _build_group_map fix. If the Fe atom's Gaussian center is
            # not present in (or has an out-of-range index into) the
            # eigenvector for this mode, the ligand's Fe-X stretch/bend
            # all become 0 without any indication in the output. Emit
            # a deduplicated UserWarning so silent zero-rows in the
            # Fe_S_Cys/Fe_N_His sheets become visible.
            _key = f"{label}|fe_lookup"
            if _key not in _WARNED_EMPTY_LIG:
                _WARNED_EMPTY_LIG.add(_key)
                _reason = ("Fe center not in eigenvector index"
                            if fe_row is None
                            else f"Fe row {fe_row} >= evg shape {evg.shape[0]}")
                warnings.warn(
                    f"Ligand '{label}': {_reason}. "
                    f"Fe-X stretch/bend will be 0 for this and any "
                    f"subsequent mode with the same problem. "
                    f"Check PDB-Gaussian atom matching.",
                    UserWarning, stacklevel=2)
            result[label] = _zero_lig(s_sb)
            continue

        # ligands-Eigenvektor
        lc_row = c2l.get(lig.lig_center)
        if lc_row is None or lc_row >= evg.shape[0]:
            # v1.0.4: same as above for the ligand donor atom.
            _key = f"{label}|lig_lookup"
            if _key not in _WARNED_EMPTY_LIG:
                _WARNED_EMPTY_LIG.add(_key)
                _reason = ("ligand center not in eigenvector index"
                            if lc_row is None
                            else f"ligand row {lc_row} >= evg shape {evg.shape[0]}")
                warnings.warn(
                    f"Ligand '{label}': {_reason}. "
                    f"Fe-X stretch/bend will be 0 for this and any "
                    f"subsequent mode with the same problem. "
                    f"Check PDB-Gaussian atom matching.",
                    UserWarning, stacklevel=2)
            result[label] = _zero_lig(s_sb)
            continue

        d_fe  = np.asarray(evg[fe_row],  dtype=float)
        d_lig = np.asarray(evg[lc_row],  dtype=float)
        bhat  = np.asarray(lig.bond_vec, dtype=float)

        if n_hat is not None:
            # v3.5: bend in INP/OOP-Fractione bzgl. cluster plane zerlegen
            split = _bend_split(d_lig, d_fe, bhat,
                                np.asarray(n_hat, dtype=float),
                                sigma_ev)
            stretch  = split["stretch"]
            bend_inp = split["bend_inp"]
            bend_oop = split["bend_oop"]
            bend_tot = math.sqrt(bend_inp * bend_inp + bend_oop * bend_oop)

            # Differenz-Signifikanz: is a the beiden Komponenten
            # signifikant groesser als the andere?
            bend_signif = classify_difference_significance(
                bend_oop, s_sb, bend_inp, s_sb, thr_low, thr_high)
            # Richtungs-Markierung anhaengen
            if bend_signif != "trivial":
                if bend_oop > bend_inp:
                    bend_signif = "OOP-" + bend_signif
                else:
                    bend_signif = "INP-" + bend_signif

            result[label] = {
                "stretch":            stretch,
                "bend":               bend_tot,
                "bend_inp":           bend_inp,
                "bend_oop":           bend_oop,
                "bend_inp_pct":       split["bend_inp_pct"],
                "bend_oop_pct":       split["bend_oop_pct"],
                "s_stretch":          s_sb,
                "s_bend":             s_sb,
                "s_bend_inp":         s_sb,
                "s_bend_oop":         s_sb,
                "bend_significance":  bend_signif,
                "element":            lig.lig_element,
                "fe_idx":             lig.fe_idx,
            }
        else:
            # Abwaertskompatibilitaet: only Skalar-bend
            rel     = d_lig - d_fe
            stretch = abs(float(rel @ bhat))
            perp    = rel - (rel @ bhat) * bhat
            bend    = float(np.linalg.norm(perp))
            result[label] = {
                "stretch":   stretch,
                "bend":      bend,
                "s_stretch": s_sb,
                "s_bend":    s_sb,
                "element":   lig.lig_element,
                "fe_idx":    lig.fe_idx,
            }

    return result


def analyze_his_hn(evg:      np.ndarray,
                   c2l:      Dict[int,int],
                   coord_info: CoordInfo,
                   cfg:      Config,
                   u_rms:    float = 1.0,
                   ) -> Dict:
    """Computes H-N bondaenderung for protonierte His-ligands.

    Parameters
    ----------
    evg : ndarray of shape (n_atoms, 3)
        Eigenvector MIT hydrogen-Atomen (thermisch scaled with u_rms by
        the caller).
    c2l : dict of {int: int}
        Mapping Gaussian-Center → Zeilenindex in ``evg``.
    coord_info : CoordInfo
        Koordinations-Information.
    cfg : Config
        Benecessaryt ``sigma_eigvec``.
    u_rms : float, optional
        Thermal scaling factor. Same value the caller used to scale the
        input ``evg``. Required for self-consistent sigma propagation:
        without this, ``s_hn_stretch`` would report the uncertainty of
        the *unscaled* eigenvector while ``hn_stretch`` reports the
        *scaled* displacement, making the two quantities incomparable
        (off by a factor of ~1/u_rms, typically ~20).

        v1.0.4: this argument is new. Pre-v1.0.4 callers (without
        ``u_rms``) get ``u_rms=1.0`` (no scaling) -- this matches the
        old behaviour bit-for-bit when ``u_rms`` was implicitly 1, but
        is still wrong if ``evg`` was scaled. Always pass the same
        ``u_rms`` you used to build ``evg``.

    Returns
    -------
    dict of {str: dict}
        Mapping Residuen-Label → ``{hn_stretch, s_hn_stretch, lig_element}``.
        Leeres Dict if no protonierten His present.
    """
    # v1.0.4 bugfix: sigma was missing the u_rms factor, while the
    # eigenvector rows ``d_n``/``d_h`` were thermally scaled. This made
    # s_hn_stretch ~20x too large for typical u_rms (~0.05 A), which is
    # the bug class identified in the v1.0.4 audit ("NH-Sigma fehlt
    # u_rms"). The corrected formula matches analyze_fe_ligand exactly:
    sigma_ev = float(cfg.sigma_eigvec) * float(u_rms)
    s_sb     = sigma_ev * math.sqrt(2.0)
    result: Dict[str, Dict] = {}

    for lig in coord_info.ligands:
        if not lig.his_protonated or lig.h_center is None:
            # Legitimate skip: this ligand is not a protonated His, so
            # there is no H-N bond to analyse. No warning.
            continue
        label = lig.res_label

        # Nutze his_hn_center (das N the H contributes, if applicable NE2 statt ND1)
        n_center = lig.his_hn_center if lig.his_hn_center is not None else lig.lig_center
        n_row = c2l.get(n_center)
        h_row = c2l.get(lig.h_center)
        if n_row is None or h_row is None:
            # v1.0.4: protonated His whose H or N is missing from c2l is a
            # genuine error (vs. legitimate deprot above, which has
            # h_center == None). Emit a deduplicated warning.
            _key = f"{label}|hn_lookup"
            if _key not in _WARNED_HN_SKIP:
                _WARNED_HN_SKIP.add(_key)
                _missing = []
                if n_row is None: _missing.append(f"N center {n_center}")
                if h_row is None: _missing.append(f"H center {lig.h_center}")
                warnings.warn(
                    f"Protonated His '{label}': "
                    f"{' and '.join(_missing)} missing from eigenvector "
                    f"index. PCET H-N stretching cannot be computed for "
                    f"this ligand. Check PDB-Gaussian atom matching.",
                    UserWarning, stacklevel=2)
            continue
        if n_row >= evg.shape[0] or h_row >= evg.shape[0]:
            # Same as above for out-of-range indices.
            _key = f"{label}|hn_oob"
            if _key not in _WARNED_HN_SKIP:
                _WARNED_HN_SKIP.add(_key)
                warnings.warn(
                    f"Protonated His '{label}': N row {n_row} or H row "
                    f"{h_row} out of range (evg shape {evg.shape[0]}). "
                    f"PCET H-N stretching cannot be computed for "
                    f"this ligand.",
                    UserWarning, stacklevel=2)
            continue

        d_n  = evg[n_row]
        d_h  = evg[h_row]
        rel  = d_h - d_n
        hn_hat = lig.hn_vec if lig.hn_vec is not None else rel / (np.linalg.norm(rel) + 1e-30)

        hn_stretch = abs(float(rel @ hn_hat))
        result[label] = {
            "hn_stretch":   hn_stretch,
            "s_hn_stretch": s_sb,
            "lig_element":  lig.lig_element,
        }

    return result


def _zero_lig(s: float) -> Dict:
    """Returns a Null-ligands-Dict (no Fe-Ligand found).

    Parameter ``s`` is Unsicherheit for stretch and bend gesetzt.
    """
    return {"stretch": 0., "bend": 0., "s_stretch": s, "s_bend": s,
            "element": "?", "fe_idx": 0}


# ===========================================================================
# Secondary structure-Analyse  (Bugfix B2: korrekte Indizierung)
# ===========================================================================

def analyze_ss_element(evg:        np.ndarray,
                        c2l:        Dict[int,int],
                        ss_centers: List[int],
                        atoms:      List[Dict],
                        idx_map:    Dict[int,int],
                        ss_type:    str,
                        u_rms:      float = 1.0,
                        sigma_eigvec: float = 5e-4,
                        ) -> Dict:
    """Analysiert Amplitude and bending a Sekundaerstruktur-elements.

    Parameters
    ----------
    evg : ndarray of shape (n_atoms, 3)
        Eigenvector the aktuellen Mode.
    c2l : dict of {int: int}
        Mapping Gaussian-Center → Zeilenindex in ``evg``.
    ss_centers : list of int
        Gaussian-Center-Nummern aller Atome dieses SS-elements.
    atoms : list of dict
        Gaussian-atom list.
    idx_map : dict of {int: int}
        Center → Atom-Index.
    ss_type : str
        ``"helix"`` or ``"sheet"``.
    u_rms : float, optional
        Thermal scaling for error propagation. Standard: ``1.0``.

    Returns
    -------
    dict
        Kennzahlen: ``amplitude_mean``, ``amplitude_max``,
        ``com_amplitude``, ``lateral_std``, ``lateral_amplitude``,
        ``stretching``, ``axial_amplitude``, ``tilting_angle``,
        ``internal_amplitude`` and zugehoerige ``s_``-Sigmas.

    Notes
    -----
    Bugfix B2: ``ss_centers`` are Gaussian-Center-Nummern (nicht
    PDB listnindizes). ``c2l[center]`` liefert the ``evg``-Zeilenindex.
    """
    _zero = {k: 0. for k in (
        "amplitude_mean","amplitude_max","com_amplitude",
        "lateral_std","lateral_amplitude","stretching",
        "axial_amplitude","tilting_angle","internal_amplitude",
        "s_amplitude_mean","s_com_amplitude","s_lateral_std",
        "s_lateral_amplitude","s_stretching","s_axial_amplitude",
        "s_tilting_angle","s_internal_amplitude")}

    if not ss_centers or evg.shape[0] == 0:
        return _zero

    # B2-Fix: Nutze c2l (Center → evg-Zeile) statt PDB listnpositionen
    valid_rows:   List[int]       = []
    valid_coords: List[List[float]] = []
    for ctr in ss_centers:
        row = c2l.get(ctr)
        if row is None or row >= evg.shape[0]:
            continue
        if ctr not in idx_map:
            continue
        a = atoms[idx_map[ctr]]
        valid_rows.append(row)
        valid_coords.append([a["x"], a["y"], a["z"]])

    if not valid_rows:
        return _zero

    e_sub = evg[valid_rows]
    c_sub = np.array(valid_coords)

    # Achse of the SS-elements (erster SVD-Vektor)
    if c_sub.shape[0] >= 2:
        ctr_c = c_sub.mean(0)
        _, _, Vt = np.linalg.svd(c_sub - ctr_c)
        axis = Vt[0]
    else:
        axis = np.array([0., 0., 1.])

    norms   = np.linalg.norm(e_sub, axis=1)
    ax_proj = e_sub @ axis
    perp    = e_sub - np.outer(ax_proj, axis)
    perp_n  = np.linalg.norm(perp, axis=1)

    amp_mean    = float(np.mean(norms))
    amp_max     = float(np.max(norms))
    com_amp     = float(np.linalg.norm(e_sub.mean(0)))
    ax_amp      = float(np.mean(np.abs(ax_proj)))
    lat_amp     = float(np.mean(perp_n))          # mittlere Querauslenkung (Å)
    lat_std     = float(np.std(perp_n))           # Std the Querauslenkung (Å)
    stretch     = float(np.std(ax_proj))           # differentielle stretching (Å)
    # Interne Deformation: Querauslenkung relativ to Schwerpunktbewegung
    e_rel       = e_sub - e_sub.mean(0)
    perp_rel    = e_rel - np.outer(e_rel @ axis, axis)
    int_amp     = float(np.mean(np.linalg.norm(perp_rel, axis=1)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        com_v  = e_sub.mean(0)
        tilt   = float(np.degrees(np.arctan2(
            float(np.linalg.norm(com_v - (com_v@axis)*axis)),
            abs(float(com_v@axis)) + 1e-30)))

    # error propagation: sigma for Meane
    # v1.0.4 bugfix: previously hardcoded as 1e-4, which silently overrode
    # the configured cfg.sigma_eigvec for SS-element sigmas only (all
    # other sigmas correctly use cfg.sigma_eigvec). Default 5e-4 here
    # matches the default Config value, so out-of-the-box behaviour is
    # 5x larger than pre-v1.0.4 -- but it is now self-consistent and
    # responds to cfg.sigma_eigvec changes in the TOML.
    n_v = len(valid_rows)
    s_amp = float(u_rms * sigma_eigvec * np.sqrt(n_v))

    return {
        "amplitude_mean":     amp_mean,  "amplitude_max":    amp_max,
        "com_amplitude":      com_amp,   "lateral_std":      lat_std,
        "lateral_amplitude":  lat_amp,   "stretching":       stretch,
        "axial_amplitude":    ax_amp,    "tilting_angle":    tilt,
        "internal_amplitude": int_amp,
        "s_amplitude_mean":   s_amp,    "s_com_amplitude":  s_amp,
        "s_lateral_std":      s_amp*1.5,"s_lateral_amplitude": s_amp,
        "s_stretching":       s_amp*1.4,"s_axial_amplitude":s_amp,
        "s_tilting_angle":    s_amp,    "s_internal_amplitude": s_amp,
    }


def analyze_all_ss(evg:          np.ndarray,
                    c2l:          Dict[int,int],
                    ss_center_map: Dict[str, List[int]],
                    atoms:        List[Dict],
                    idx_map:      Dict[int,int],
                    ss_elements:  List[Dict],
                    u_rms:        float = 1.0,
                    sigma_eigvec: float = 5e-4,
                    ) -> Dict[str, Dict]:
    """Analysiert alle SS-elements for a Mode.

    Parameters
    ----------
    evg : ndarray of shape (n_atoms, 3)
        Eigenvector the aktuellen Mode.
    c2l : dict of {int: int}
        Mapping Gaussian-Center → Zeilenindex in ``evg``.
    ss_center_map : dict of {str: list of int}
        SS-element-Name → Gaussian-Center-Nummern.
    atoms : list of dict
        Gaussian-atom list.
    idx_map : dict of {int: int}
        Center → Atom-Index.
    ss_elements : list of dict
        SS-element-Records from ``parse_pdb``.
    u_rms : float, optional
        Thermal scaling. Standard: ``1.0``.

    Returns
    -------
    dict of {str: dict}
        Mapping SS-element-Name → Analyse-Kennzahlen.
    """
    out: Dict[str, Dict] = {}
    for elem in ss_elements:
        name    = elem["name"]
        centers = ss_center_map.get(name, [])
        out[name] = analyze_ss_element(
            evg, c2l, centers, atoms, idx_map, elem["type"], u_rms,
            sigma_eigvec=sigma_eigvec)
    return out


# ===========================================================================
# OOP / INP-Klassifikation (binaer + fein, Hardening v3.0 #6)
# ===========================================================================

def classify_oop_inp(oop_pct: float,
                     binary_threshold: float = 60.0,
                     detail_thresholds: Tuple[float, float, float] = (60.0, 75.0, 90.0),
                     ) -> Tuple[str, str]:
    r"""Classifies a Mode anhand ihres OOP-Fractions in zwei Stufen.

    Aus the Out-of-plane-Quadratanteil ``oop_pct`` (in %) are zwei
    Labels abgeleitet:

      * **Broad-Label** (3 Stufen, Default-Schwelle 60 %):
        ``"Out-of-plane"`` / ``"In-plane"`` / ``"Torsional/Mixed"``.
        Identisch to frueheren Logik before Hardening v3.0; is fuer
        Excel-Faerbung and Plot-Legenden weiter used.

      * **Detail-Label** (7 Stufen, drei symmetrische Schwellen
        ``(low, mid, high)``, default ``(60, 75, 90)``):

        +-------------------+-------------------+
        | OOP%              | Label             |
        +===================+===================+
        | OOP% >= high      | "Pure OOP"         |
        +-------------------+-------------------+
        | mid <= OOP% < high| "Strong OOP"       |
        +-------------------+-------------------+
        | low <= OOP% < mid | "Majority OOP"|
        +-------------------+-------------------+
        | (100-low) <       | "Mixed"        |
        | OOP% < low        |                   |
        +-------------------+-------------------+
        | (100-mid) <       | "Majority INP"|
        | OOP% <= (100-low) |                   |
        +-------------------+-------------------+
        | (100-high) <      | "Strong INP"       |
        | OOP% <= (100-mid) |                   |
        +-------------------+-------------------+
        | OOP% <= (100-high)| "Pure INP"         |
        +-------------------+-------------------+

    Symmetry: the Detail-Label is symmetrisch um 50 %. For OOP% = 50
    erhaelt man "Mixed"; for OOP% = 100 "Pure OOP"; for OOP% = 0
    "Pure INP". The three thresholds ``low/mid/high`` liegen alle
    above 50 %.

    Anwendung in the PCET-Kontext: the Detail-Label hilft, modes mit
    different starker OOP-Mischung to unterscheiden — z. B.
    "stark OOP" (75-90 %) is with hoher Wahrscheinlichkeit eine
    cluster-breathingmode with Out-of-plane-Komponente, "mehrheitlich
    OOP" (60-75 %) hingegen typically a gemischte
    motion with klarem Schwerpunkt.

    Parameters
    ----------
    oop_pct : float
        Out-of-plane-Fraction in Prozent, ``0 <= oop_pct <= 100``.
    binary_threshold : float, optional
        Threshold for the Broad-Label in Prozent. default 60.
    detail_thresholds : (float, float, float), optional
        Drei aufsteigende Schwellen for the Detail-Label in Prozent.
        default (60, 75, 90).

    Returns
    -------
    broad : str
        ``"Out-of-plane"``, ``"In-plane"`` or ``"Torsional/Mixed"``.
    detail : str
        One of 7 labels: ``"Pure OOP"``, ``"Strong OOP"``,
        ``"Majority OOP"``, ``"Mixed"``, ``"Majority INP"``,
        ``"Strong INP"``, ``"Pure INP"``.
    """
    inp_pct = 100.0 - oop_pct

    # Broad-Label (3 Stufen)
    if oop_pct >= binary_threshold:
        broad = "Out-of-plane"
    elif inp_pct >= binary_threshold:
        broad = "In-plane"
    else:
        broad = "Torsional/Mixed"

    # Detail-Label (7 Stufen)
    low, mid, high = detail_thresholds
    if oop_pct >= high:
        detail = "Pure OOP"
    elif oop_pct >= mid:
        detail = "Strong OOP"
    elif oop_pct >= low:
        detail = "Majority OOP"
    elif oop_pct > 100.0 - low:
        detail = "Mixed"
    elif oop_pct > 100.0 - mid:
        detail = "Majority INP"
    elif oop_pct > 100.0 - high:
        detail = "Strong INP"
    else:
        detail = "Pure INP"

    return broad, detail


# ===========================================================================
# Kern-Klassifikation (D2h-Symmetriekoordinaten)
# ===========================================================================

def classify_kernel_mode_from_evg(evg_4x3:    np.ndarray,
                                    atoms:      List[Dict],
                                    idx_map:    Dict[int,int],
                                    fe_c:       List[int],
                                    s_c:        List[int],
                                    normal:     np.ndarray,
                                    min_total_disp: float = 0.003,
                                    ) -> Tuple[str, str, Dict, float]:
    r"""Heuristische geometrische Klassifikation the Kern-motion.

    .. warning::
       **Diese Funktion is HEURISTISCH.** Die zurueckgivenen Scores
       are geometrische Projektionen on qualitative motionsmuster
       and correspond to *nicht* the orthogonalen D2h-Irrep-Zerlegung.
       Mehrere Scores koennen to selben Irrep gehoeren (z. B. sind
       Breathing, Fe-stretching and S-stretching alle Ag-Fractione), und
       einige Scores mischen mehrere Irreps (z. B. Hinge-Fe enthaelt
       sowohl B2g-Rotation als also Au-Antisym-Twist).

       **Fuer a rigorose Symmetriezerlegung** (orthogonal,
       Summenregel, vergleichbar between Strukturen) verwende die
       SCSD-Methode in :func:`compute_scsd_for_mode_full` nach
       Kingsbury & Senge, *Chem. Sci.* **15**, 13638 (2024). Siehe
       also :func:`extract_dominant_scsd_irreps` for the SCSD-basierte
       Klassifikation.

    Achsenkonvention
    ----------------
    Die Heuristik nimmt an, dass the Cluster-Koordinatensystem so
    orientiert ist:

    - **x** entlang the Fe-Fe-Achse (Fe1 for -d_Fe, Fe2 for +d_Fe)
    - **y** entlang the S-S-Achse (S1 for +d_S, S2 for -d_S)
    - **z** = ``normal`` (cluster normal, out-of-plane)

    Diese Konvention stimmt with the kanonischen SCSD-Referenzgeometrie
    ueberein (see ``_SCSD_MODEL_COORDS``).

    D2h-Irrep-Zuordnung the Heuristik-Scores
    ----------------------------------------
    Eingangs-Verifikation per Symmetrieanalyse (auf reine D2h-Moden
    angewandt; see Tests):

    +------------------+--------------------------------------------------+
    | Heuristik-Score  | D2h-Irrep(s) (verifiziert)                       |
    +==================+==================================================+
    | Translation      | B1u + B2u + B3u (Schwerpunktsbewegung; for       |
    |                  | sauberen Eckart-orthogonalen modes ~ 0)          |
    +------------------+--------------------------------------------------+
    | Umbrella         | B1u (z-Komponente the Translation; ~ 0 for      |
    |                  | reine Vibrationen)                               |
    +------------------+--------------------------------------------------+
    | Hinge-Fe         | B2g (Rotation um y) + Au-Fractione                 |
    |                  | (sollte for Vibrationen ~ 0 sein)               |
    +------------------+--------------------------------------------------+
    | Hinge-S          | B3g (Rotation um x) + Au-Fractione                 |
    |                  | (sollte for Vibrationen ~ 0 sein)               |
    +------------------+--------------------------------------------------+
    | OOP-Twist        | Mischung from B2g + B3g (Mittel over             |
    |                  | Hinge-Fe + Hinge-S)                              |
    +------------------+--------------------------------------------------+
    | Fe-stretching     | Ag (Fe-Fe symmetrische stretching)                |
    +------------------+--------------------------------------------------+
    | S-stretching      | Ag (S-S symmetrische stretching)                  |
    +------------------+--------------------------------------------------+
    | Breathing        | Ag (radiale Mischung from Fe-Fe + S-S stretching)  |
    +------------------+--------------------------------------------------+
    | Rhombus-shear | naeherungsweise B1g (Rhombus-shear in the Ebene); |
    |                  | numerisch not unique isoliert               |
    +------------------+--------------------------------------------------+
    | Rotation-ip      | Mischung from Cluster-Rotation um z (B1g) and     |
    |                  | rhombischer In-plane-Scherung (B1g)              |
    +------------------+--------------------------------------------------+

    Diagnostische Bedeutung
    -----------------------
    Translation, Umbrella, Hinge-Fe, Hinge-S sollten for saubere
    Eckart-orthogonale DFT-Vibrationen *nahe Null* sein. Nichtnull-values
    deuten an, dass the Mode (a) translatorisch/rotatorisch verunreinigt
    ist, or (b) the cluster als Teil a groesseren
    Proteinbewegung mitbewegt — beides is diagnostisch wertvoll fuer
    PCET-relevante Niederfrequenzmoden.

    Fe-stretching, S-stretching, Breathing, Rhombus-shear are die
    *physikalisch* interessanten Fractione (echte interne Vibrationen).

    Parameters
    ----------
    evg_4x3 : ndarray of shape (4, 3)
        displacement vectors the vier Clusteratome (Fe1, Fe2, S1, S2).
    atoms : list of dict
        Gaussian-atom list.
    idx_map : dict of {int: int}
        Center → Atom-Index.
    fe_c : list of int
        Center-Nummern the Fe atoms.
    s_c : list of int
        Center-Nummern the S-Atome.
    normal : ndarray of shape (3,)
        normal vector the clusterebene.
    min_total_disp : float, optional
        Mindest-Gesamtauslenkung for sinnvolle Klassifikation.
        Standard: ``0.003``.

    Returns
    -------
    primary : str
        Heuristisches Primaer-Label. Mit the hoechsten Score.
    secondary : str
        Heuristisches Sekundaer-Label.
    scores : dict of {str: float}
        Heuristische Scores aller motionsmuster (values in [0, 1],
        not orthogonal, summieren i. A. not to 1).
    total_disp : float
        Gesamtauslenkung aller vier Clusteratome.

    See Also
    --------
    compute_scsd_for_mode_full : Rigorose D2h-Symmetriezerlegung
        after Kingsbury & Senge.
    extract_dominant_scsd_irreps : Determines primaere/sekundaere
        Irrep from SCSD-Output.
    """
    _norm_n = float(np.linalg.norm(normal))
    if _norm_n < 1e-12:
        # Cluster-Atome are kollinear; Normale is not defined.
        # SCSD kann not laufen without valide cluster plane.
        return "No core", "-", {}, 0.0
    n_hat = normal / _norm_n
    cl4   = fe_c + s_c

    if evg_4x3.shape[0] < 4 or not all(c in idx_map for c in cl4):
        return "No core", "-", {}, 0.0

    r = np.array([[atoms[idx_map[c]]["x"], atoms[idx_map[c]]["y"],
                   atoms[idx_map[c]]["z"]] for c in cl4])
    d = evg_4x3

    total_disp = float(np.linalg.norm(d))
    if total_disp < min_total_disp:
        return "No core contribution", "-", {}, total_disp

    fe1, fe2 = r[0], r[1]
    s1,  s2  = r[2], r[3]
    d_fe1, d_fe2 = d[0], d[1]
    d_s1,  d_s2  = d[2], d[3]

    fe_ax   = (fe2 - fe1); fe_ax_l = np.linalg.norm(fe_ax)
    if fe_ax_l < 1e-12: return "No core contribution", "-", {}, total_disp
    fe_hat  = fe_ax / fe_ax_l

    oop_ax  = n_hat
    ip_ax   = np.cross(oop_ax, fe_hat)
    ip_l    = np.linalg.norm(ip_ax)
    if ip_l < 1e-12: return "No core contribution", "-", {}, total_disp
    ip_hat  = ip_ax / ip_l

    ctr     = r.mean(0)
    def proj(v, ax):
        """Skalare Projektion von v on Einheitsvektor ax."""
        return float(v @ ax)
    def safe_hat(v):
        """Einheitsvektor von v; Nullvektor if |v| < 1e-12."""
        l = np.linalg.norm(v)
        return v/l if l > 1e-12 else np.zeros(3)

    scores: Dict[str, float] = {}
    com = (d_fe1+d_fe2+d_s1+d_s2)/4
    scores["Translation"]     = float(np.linalg.norm(com)) / total_disp

    oop_vals = [proj(d_fe1,oop_ax), proj(d_fe2,oop_ax),
                proj(d_s1,oop_ax), proj(d_s2,oop_ax)]
    scores["Umbrella"]        = abs(float(np.mean(oop_vals))) / total_disp
    scores["Hinge-Fe"]        = abs(proj(d_fe1,oop_ax)-proj(d_fe2,oop_ax))/2/total_disp
    scores["Hinge-S"]         = abs(proj(d_s1,oop_ax)-proj(d_s2,oop_ax))/2/total_disp
    scores["OOP-Twist"]       = (abs(proj(d_fe1,oop_ax)-proj(d_fe2,oop_ax))+
                                  abs(proj(d_s1,oop_ax)-proj(d_s2,oop_ax)))/4/total_disp
    scores["Fe-stretching"]    = abs(proj(d_fe1,fe_hat)-proj(d_fe2,fe_hat))/2/total_disp

    s_ax = s2-s1; s_l = np.linalg.norm(s_ax)
    if s_l > 1e-12:
        s_hat = s_ax/s_l
        scores["S-stretching"] = abs(proj(d_s1,s_hat)-proj(d_s2,s_hat))/2/total_disp
    else:
        scores["S-stretching"] = 0.

    breathing = (proj(d_fe1,safe_hat(fe1-ctr))+proj(d_fe2,safe_hat(fe2-ctr))+
                 proj(d_s1, safe_hat(s1-ctr)) +proj(d_s2, safe_hat(s2-ctr)))/4
    scores["Breathing"]       = abs(breathing) / total_disp

    d1=fe2-s1; d2=fe1-s2; d1l=np.linalg.norm(d1); d2l=np.linalg.norm(d2)
    if d1l > 1e-12 and d2l > 1e-12:
        sc = (proj(d_fe2-d_s1,d1/d1l)-proj(d_fe1-d_s2,d2/d2l))/2
        scores["Rhombus-shear"] = abs(sc)/total_disp
    else:
        scores["Rhombus-shear"] = 0.

    # Rotation-ip: Drehmoment on cluster normal n_hat projiziert
    # (Bug-Fix: frueherer Code nutzte globale z-Achse [2], unabhaengig von Orientierung)
    rot_ip = [
        float(np.dot(np.cross(fe1 - ctr, d_fe1), n_hat)),
        float(np.dot(np.cross(fe2 - ctr, d_fe2), n_hat)),
        float(np.dot(np.cross(s1  - ctr, d_s1 ), n_hat)),
        float(np.dot(np.cross(s2  - ctr, d_s2 ), n_hat)),
    ]
    scores["Rotation-ip"] = abs(float(np.mean(rot_ip))) / total_disp

    for k in scores:
        scores[k] = min(1.0, max(0.0, scores[k]))

    sorted_s  = sorted(scores.items(), key=lambda x: -x[1])
    primary   = sorted_s[0][0] if sorted_s else "n/a"
    secondary = sorted_s[1][0] if len(sorted_s) > 1 else "-"
    return primary, secondary, scores, total_disp


# ===========================================================================
# Haupt-Modusanalyse
# ===========================================================================

# SCORE_KEYS defined near top of module (after MAX_REASONABLE_KERNEL_D)


def evg_sub_extern(evg: np.ndarray, c2l: Dict[int, int],
                   ctrs: List[int]) -> np.ndarray:
    """Returns Eigenvektor-Zeilen for givene Center-Nummern .

    Sicherer Zugriff: fehlende Center or Index outside von evg
    are skipped. Ersetzt the frueheren Monkey-Patch.
    """
    idxs = [c2l[c] for c in ctrs if c in c2l and c2l[c] < evg.shape[0]]
    return evg[idxs] if idxs else np.zeros((0, 3))


def analyze_mode(bi:          BlockInfo,
                  col:         int,
                  filepath:    str,
                  atoms:       List[Dict],
                  idx_map:     Dict[int,int],
                  normal:      np.ndarray,
                  coord_info:  CoordInfo,
                  fe_c:        List[int],
                  s_c:         List[int],
                  cfg:         Config,
                  ) -> Optional[Dict]:
    """Performs the completee Analyse a einzelnen normal mode durch.

    Parameters
    ----------
    bi : BlockInfo
        Metadaten of the zugehoerigen Frequenz-Blocks.
    col : int
        Spalten-Index (0-basiert) the mode in the Block.
    filepath : str
        Path to the Gaussian ``.log``-file.
    atoms : list of dict
        Gaussian-atom list (ohne hydrogen).
    idx_map : dict of {int: int}
        Center → Atom-Index.
    normal : ndarray of shape (3,)
        normal vector the clusterebene.
    coord_info : CoordInfo
        Koordinations-Information.
    fe_c : list of int
        Center-Nummern the Fe atoms.
    s_c : list of int
        Center-Nummern the S-Atome.
    cfg : Config
        Konfiguration.

    Returns
    -------
    dict or None
        Ergebnis-Dict with allen Kennzahlen, or ``None`` for Error.
    """
    # Eigenvector laden (ohne H)
    centers, evg = _get_eigvec_smart(filepath, bi, col,
                                       include_hydrogen=False, cfg=cfg)
    if evg.shape[0] == 0:
        return None

    freq      = bi.freqs[col]      if col < len(bi.freqs)      else 0.0
    red_mass  = bi.red_masses[col] if col < len(bi.red_masses) else 1.0

    # Thermal scaling (Bugfix B3)
    u_rms = compute_thermal_amplitude(freq, red_mass, cfg.temp_k, cfg.amplitude)
    evg   = evg * u_rms

    c2l    = {c: i for i, c in enumerate(centers)}
    _norm_n2 = float(np.linalg.norm(normal))
    if _norm_n2 < 1e-12:
        # Defensiver Schutz: for degenerierter Cluster-Geometrie
        # (sollte through find_cluster bereits abgefangen sein, aber
        # wir vermeiden hier a Division through 0 zwingend).
        return None
    n_hat  = normal / _norm_n2
    _s     = cfg.sigma_eigvec * u_rms

    def _evg_sub(ctrs: List[int]) -> np.ndarray:
        """Extracts Eigenvektorzeileb for the givenen Gaussian-Center-Nummern."""
        idxs = [c2l[c] for c in ctrs if c in c2l]
        return evg[idxs] if idxs else np.zeros((0, 3))

    # ===== Ring 1: cluster core (Fe + S) ======================================
    evg_cl = _evg_sub(fe_c + s_c)
    if evg_cl.shape[0] > 0:
        kern_oop = _oop(evg_cl, n_hat) * 100.0
        kern_d   = _kd(evg_cl, None)
        s_ko     = _sigma_oop(evg_cl, n_hat, _s) * 100.0
        s_kd     = _sigma_kd(evg_cl, _s)
        if kern_d > MAX_REASONABLE_KERNEL_D:
            kern_d = kern_oop = s_ko = s_kd = float("nan")
    else:
        kern_oop = kern_d = s_ko = s_kd = 0.0

    # ===== Ring 2: Cluster-ligands-bond Atome ============================
    # (Cys-SG, His-ND1/NE2, Asp/Glu-O, ...). Ueber alle ligand.lig_center.
    # These Atome modulieren direkt the Fe-ligands-bonden and sind
    # spektroskopisch in the Fe-pDOS sichtbar (NRVS-relevante Banden 200-450 cm-1).
    lig_centers_uniq: List[int] = []
    _seen_lc: set = set()
    for _lig in coord_info.ligands:
        if _lig.lig_center not in _seen_lc:
            _seen_lc.add(_lig.lig_center)
            lig_centers_uniq.append(_lig.lig_center)
    evg_lig = _evg_sub(lig_centers_uniq)
    ring2 = _oop_ring_metrics(evg_lig, n_hat, _s)
    lig_oop_pct  = ring2["oop_pct"]
    lig_inp_pct  = ring2["inp_pct"]
    lig_d        = ring2["d"]
    s_lig_oop    = ring2["sigma_oop_pct"]
    s_lig_d      = ring2["sigma_d"]

    # ===== Ring 3: secondary sphere ==========================================
    # Vereinigung aller voll-ligands-Reste (group_map) + alle PCET-H-Bond-
    # acceptoren. Beschreibt the elektrostatische and H-Bond-Umfeld der
    # ligands, the for PCET-modulation and Resonanz-Raman-Banden relevant
    # ist. Atome are dedupliziert; Ring 3 enthaelt typically auch
    # Ring 2 als Untermenge (deliberately, da the ligands-bond Atome Teil
    # ihrer amino acid-Reste sind).
    second_centers: List[int] = []
    _seen_sc: set = set()
    for _gctr in coord_info.group_map.values():
        for _c in _gctr:
            if _c not in _seen_sc:
                _seen_sc.add(_c); second_centers.append(_c)
    if coord_info.pcet_info is not None:
        _accs = getattr(coord_info.pcet_info, "acceptor_centers_per_h", []) or []
        for _aclist in _accs:
            for _c in _aclist:
                if _c not in _seen_sc:
                    _seen_sc.add(_c); second_centers.append(_c)
    evg_2nd = _evg_sub(second_centers)
    ring3 = _oop_ring_metrics(evg_2nd, n_hat, _s)
    second_oop_pct  = ring3["oop_pct"]
    second_inp_pct  = ring3["inp_pct"]
    second_d        = ring3["d"]
    s_second_oop    = ring3["sigma_oop_pct"]
    s_second_d      = ring3["sigma_d"]

    # ===== Mode type: jetzt basierend on ligands-Sphaere (Ring 2) ============
    # Klassifikation on Basis von lig_oop_pct is universell also fuer
    # Cluster without PCET-Muster (reine Cys4-Cluster wie Ferredoxin) sinnvoll,
    # da alle [2Fe-2S]-Cluster ligands besitzen. Die globale OOP-Statistik
    # over alle Atome is not mehr computed (sie war for delokalisierten
    # Proteinmoden not aussagekraeftig).
    mode_type, mode_type_detail = classify_oop_inp(
        lig_oop_pct,
        binary_threshold=cfg.mode_type_threshold,
        detail_thresholds=cfg.mode_type_detail_thresholds,
    )

    # Cluster COM / Expansion / Rotation
    cl_com = cl_exp = cl_rot = 0.0
    if evg_cl.shape[0] >= 2 and math.isfinite(kern_d):
        cl_com = float(np.linalg.norm(evg_cl.mean(0)))
        coords_cl = np.array(
            [[atoms[idx_map[c]]["x"], atoms[idx_map[c]]["y"], atoms[idx_map[c]]["z"]]
             for c in fe_c + s_c if c in idx_map])
        if len(coords_cl) >= 2:
            ctr_cl = coords_cl.mean(0)
            rad    = coords_cl - ctr_cl
            rlen   = np.linalg.norm(rad, axis=1, keepdims=True)
            rlen[rlen < 1e-15] = 1.0
            r_hat  = rad / rlen
            nr     = min(evg_cl.shape[0], r_hat.shape[0])
            cl_exp = float(np.mean(np.sum(evg_cl[:nr] * r_hat[:nr], axis=1)))
            if (len(fe_c) >= 2 and
                    all(c in idx_map for c in fe_c[:2])):
                fe1c = np.array([atoms[idx_map[fe_c[0]]]["x"],
                                  atoms[idx_map[fe_c[0]]]["y"],
                                  atoms[idx_map[fe_c[0]]]["z"]])
                fe2c = np.array([atoms[idx_map[fe_c[1]]]["x"],
                                  atoms[idx_map[fe_c[1]]]["y"],
                                  atoms[idx_map[fe_c[1]]]["z"]])
                ax  = fe2c - fe1c
                al  = np.linalg.norm(ax)
                if al > 1e-12:
                    ax  /= al
                    rc   = evg_cl[:nr] - np.outer(evg_cl[:nr] @ ax, ax)
                    cl_rot = float(np.mean(np.linalg.norm(rc, axis=1)))

    # Gruppen
    group_res: Dict[str, Dict] = {}
    for gname, gctr in coord_info.group_map.items():
        evg_g = _evg_sub(gctr)
        if evg_g.shape[0] == 0:
            # Diagnostic warning (v1.0.2): the residue is in the group_map
            # but none of its Gaussian centers are present in the current
            # eigenvector mapping (c2l). With the v1.0.2 _build_group_map
            # fix this should be impossible for canonical Cys/His residues,
            # but we keep the branch as a safety net and emit a warning so
            # the same class of silent zero-row bug can be caught early in
            # future systems (e.g. Asp/Glu/Ser/Thr ligation, unusual PDBs).
            if gname not in _WARNED_EMPTY_GROUP:
                _WARNED_EMPTY_GROUP.add(gname)
                n_in_c2l = sum(1 for c in gctr if c in c2l)
                warnings.warn(
                    f"Group '{gname}': {len(gctr)} centers in group_map, "
                    f"but {n_in_c2l} of them are in the eigenvector index. "
                    f"All Groups_OOP/INP/Winkel/Tors values for this group "
                    f"will be 0. Check PDB-Gaussian atom matching.",
                    UserWarning, stacklevel=2)
            group_res[gname] = {k: 0. for k in
                ("oop","inp","angle","torsion","total","s_oop","s_inp","s_angle","s_tors")}
            continue
        g_oop  = _oop(evg_g, n_hat) * 100.0
        g_ang  = _ang(evg_g, n_hat, cfg.amplitude_threshold)
        g_tot  = float(np.sqrt(np.sum(evg_g**2)))
        s_oop  = _sigma_oop(evg_g, n_hat, _s) * 100.0
        s_ang  = _sigma_ang(evg_g, n_hat, cfg.amplitude_threshold, _s)

        # Torsion
        cg    = [np.array([atoms[idx_map[c]]["x"], atoms[idx_map[c]]["y"],
                            atoms[idx_map[c]]["z"]])
                 for c in gctr if c in idx_map]
        ctr_g = np.mean(cg, 0) if cg else np.zeros(3)
        tors  = []
        # v1.0.4: index drift bugfix.
        # evg_g is built via _evg_sub(gctr), which filters out gctr
        # entries not in c2l. So evg_g has shape (n_in_c2l, 3) and its
        # i-th row corresponds to the i-th gctr center FOR WHICH
        # c IN c2l. The original loop used `ai` (the gctr index) to
        # address both gctr AND evg_g, which silently mis-aligns rows
        # whenever c2l is missing any center earlier in gctr. We now
        # maintain an explicit ai_local counter that advances only for
        # centers we actually consumed from evg_g, guaranteeing that
        # rv (computed from atom coords) and evg_g[ai_local] always
        # refer to the SAME atom. In normal operation c2l covers all
        # heavy-block centers and the fix is a no-op; if c2l ever has
        # gaps it prevents cross-atom contamination of the torsion sum.
        ai_local = 0
        for c in gctr:
            if c not in c2l or c not in idx_map:
                continue
            if ai_local >= evg_g.shape[0]:
                break
            rv = np.array([atoms[idx_map[c]]["x"], atoms[idx_map[c]]["y"],
                            atoms[idx_map[c]]["z"]]) - ctr_g
            rl = np.linalg.norm(rv)
            if rl < 1e-10:
                ai_local += 1
                continue
            tv = np.cross(n_hat, rv/rl); tl = np.linalg.norm(tv)
            if tl > 1e-10:
                tors.append(abs(float(evg_g[ai_local] @ (tv/tl))))
            ai_local += 1
        g_tors = float(np.mean(tors)) if tors else 0.0

        group_res[gname] = {
            "oop": g_oop, "inp": 100.-g_oop,
            "angle": g_ang, "torsion": g_tors, "total": g_tot,
            "s_oop": s_oop, "s_inp": s_oop,
            "s_angle": s_ang, "s_tors": 0.,
        }

    # Fe-ligands-Analyse (v3.5: bend in INP/OOP gesplittet)
    fe_lig_res = analyze_fe_ligand(evg, c2l, coord_info, atoms, idx_map, cfg,
                                    n_hat=n_hat, u_rms=u_rms)

    # H-N Analyse (nur if H-Atome present and His protoniert)
    #
    # Wichtig for PCET-Modenklassifikation: for protonierten His-ligands
    # is the N-H-Streckkoordinate erfasst and in the Mode-Composition-
    # Output uebernommen. A stille Fehlbehandlung an dieser Stelle
    # wuerde PCET-relevante Information verlieren — therefore loggt der
    # except-Zweig the Frequency and Modennummer the betroffenen Mode,
    # so that the REPORT nachvollziehbar bleibt.
    his_hn_res: Dict[str,Dict] = {}
    if cfg.include_hn_vibration and coord_info.his_ligand_labels:
        # Eigenvector MIT H nachladen
        try:
            centers_h, evg_h = _get_eigvec_smart(
                filepath, bi, col, include_hydrogen=True, cfg=cfg)
            c2l_h = {c: i for i, c in enumerate(centers_h)}
            # v1.0.4: pass u_rms so analyze_his_hn can scale its sigma
            # consistently with the (u_rms-scaled) eigenvector input.
            his_hn_res = analyze_his_hn(
                evg_h * u_rms, c2l_h, coord_info, cfg, u_rms=u_rms)
        except Exception as _exc_hn:
            # Bug v3.0 (vormals undefined: ``runlog`` and ``r``):
            # warnings.warn is the in the Modul uebliche Konvention; das
            # main script faengt UserWarnings ab and routet sie ueber
            # ``runlog`` weiter. Modennummer and Frequency from ``bi``/``col``
            # extrahieren — beide are hier definitiv in the Scope, ``r`` noch
            # nicht.
            import warnings as _w
            _mn = bi.mode_nums[col] if col < len(bi.mode_nums) else "?"
            _w.warn(
                f"H-N-Analyse fehlgeschlagen "
                f"(Mode {_mn} @ {freq:.2f} cm-1): "
                f"{type(_exc_hn).__name__}: {_exc_hn}. "
                f"PCET-relevante His-N-H-stretching fehlt for diese Mode.",
                UserWarning, stacklevel=2)

    # Kern-Klassifikation
    kern_primary = kern_secondary = "n/a"
    kern_scores: Dict[str, float] = {}
    pts_ref_cl = pts_dist_cl = None
    cl4 = fe_c + s_c
    if all(c in idx_map and c in c2l for c in cl4):
        pts_ref_cl  = np.array([[atoms[idx_map[c]]["x"],
                                   atoms[idx_map[c]]["y"],
                                   atoms[idx_map[c]]["z"]] for c in cl4])
        pts_dist_cl = pts_ref_cl.copy()
        for ci_i, c in enumerate(cl4):
            li = c2l[c]
            if li < evg.shape[0]:
                pts_dist_cl[ci_i] += evg[li]
        evg_4x3 = pts_dist_cl - pts_ref_cl
        kern_primary, kern_secondary, kern_scores, _ = \
            classify_kernel_mode_from_evg(
                evg_4x3, atoms, idx_map, fe_c, s_c, normal)

    # Kern-Lokalisierungsgrad: Fraction the Quadratsumme in the [2Fe-2S]-Kern
    # (= Fraction the kinetischen energy for mass-weighted Cartesian Eigenvektoren)
    # value 0.0 = no Kern-Lokalisierung, 1.0 = vollständig in the Kern
    kern_loc = 0.0
    for _kc in (fe_c + s_c):
        _kli = c2l.get(_kc)
        if _kli is not None and _kli < evg.shape[0]:
            kern_loc += float(np.sum(evg[_kli] ** 2))

    # Normalisiere kern_loc on [0, 1] relativ to Gesamt-Quadratsumme
    _evg_total = float(np.sum(evg ** 2))
    kern_loc_frac = (kern_loc / _evg_total) if _evg_total > 1e-30 else 0.0
    kern_loc_frac = min(1.0, max(0.0, kern_loc_frac))

    # ===================================================================
    # v3.7: Marcus-Hush-Reorganisations-Modulationen per Kanal
    # ===================================================================
    # We compute per Mode the geometrische modulation dr_X and den
    # Marcus-Hush-reorganization contribution lambda_X for alle bond 
    # channels X in {FeFe, FeN, FeS, NH, HA}. NO classification,
    # NO thresholds, NO filter - only the raw modulation data.
    # Die system aggregation is done later in the runner.
    reorg_per_mode: Dict[str, Dict[str, float]] = {}
    reorg_subchannels: list = []

    if getattr(cfg, "pcet_enabled", False):
        from . import reorganization as reorg_mod

        # We need the H-Eigenvector (Cluster-Modes with H-Atomen).
        if (coord_info.et_info is not None
                and coord_info.atoms_h is not None):
            try:
                centers_h_pcet, evg_h_pcet = _get_eigvec_smart(
                    filepath, bi, col, include_hydrogen=True, cfg=cfg)
                c2l_h_pcet = {c: i for i, c in enumerate(centers_h_pcet)}

                # Channels einmalig cachen via coord_info-Attribut
                cached_channels = getattr(coord_info, "_reorg_channels_v37", None)
                if cached_channels is None:
                    r0 = float(getattr(cfg, "pcet_acceptor_r0_a", 2.8))
                    sig = float(getattr(cfg, "pcet_acceptor_sigma_a", 0.4))
                    cached_channels = reorg_mod.build_channels_from_coord_info(
                        coord_info=coord_info,
                        atoms_h=coord_info.atoms_h,
                        idx_map_h=coord_info.idx_map_h,
                        acceptor_r0_a=r0, acceptor_sigma_a=sig,
                    )
                    coord_info._reorg_channels_v37 = cached_channels

                # Eigenvector thermisch skalieren -> displacementen in A
                e_atoms_full = evg_h_pcet * u_rms

                # Pro-Mode-Berechnung
                ch_results = reorg_mod.compute_mode_modulations(
                    e_atoms=e_atoms_full,
                    omega_cm1=freq,
                    mode_red_mass_amu=red_mass,
                    channels=cached_channels,
                )

                # Aggregation per parent_channel (kompakte Diagnose-Sicht)
                reorg_per_mode = reorg_mod.aggregate_by_parent(ch_results)
                
                # Sub-Channel-contributions als kompakte Liste mitschleppen
                # (nur Name + lambda_pair + lambda_mode, the spart Speicher
                # gegenvia the vollen ChannelResult). Is in the runner zu
                # system total per Sub-Channel aggregiert.
                reorg_subchannels = [
                    (r.name, r.parent_channel, r.weight,
                     r.lambda_pair_cm1, r.lambda_mode_cm1, r.dr_signed_a)
                    for r in ch_results
                    if np.isfinite(r.lambda_pair_cm1)
                ]

            except Exception as _exc_reorg:
                # For jeder Errorquelle: leere Aggregate, weiterlaufen
                from . import reorganization as _rm
                reorg_per_mode = {ch: {
                    "dr_rss_a": float("nan"), "dr_sum_signed_a": float("nan"),
                    "lambda_pair_cm1": float("nan"),
                    "lambda_mode_cm1": float("nan"),
                    "n_subchannels": 0,
                } for ch in _rm.CHANNELS}
                reorg_subchannels = []
                import warnings as _w_reorg
                _mn_p = bi.mode_nums[col] if col < len(bi.mode_nums) else "?"
                _w_reorg.warn(
                    f"Reorg computation failed "
                    f"(Mode {_mn_p} @ {freq:.2f} cm-1): "
                    f"{type(_exc_reorg).__name__}: {_exc_reorg}.",
                    UserWarning, stacklevel=2)

    mn = bi.mode_nums[col] if col < len(bi.mode_nums) else -1

    return {
        "number":    mn,
        "freq":      freq,
        "red_mass":  red_mass,
        "frc_const": bi.frc_consts[col] if col < len(bi.frc_consts) else 0.,
        "sym":       bi.syms[col]       if col < len(bi.syms)       else "A",
        "precision": "high" if bi.is_hp else "standard",
        "mode_type": mode_type,
        "mode_type_detail": mode_type_detail,
        # Ring 2: Cluster-ligands-bond Atome (v3.5)
        "lig_oop_pct":   lig_oop_pct,   "lig_inp_pct":   lig_inp_pct,
        "lig_d":         lig_d,
        "s_lig_oop":     s_lig_oop,     "s_lig_d":       s_lig_d,
        # Ring 3: secondary sphere (v3.5)
        "second_oop_pct": second_oop_pct, "second_inp_pct": second_inp_pct,
        "second_d":       second_d,
        "s_second_oop":   s_second_oop,   "s_second_d":     s_second_d,
        # Ring 1: cluster core (Fe + S)
        "kern_oop":  kern_oop,  "kern_inp":  100.-kern_oop,
        "kern_d":    kern_d,
        "s_kern_oop": s_ko,    "s_kern_d":  s_kd,
        "cl_com":    cl_com,   "cl_exp":    cl_exp,  "cl_rot": cl_rot,
        "groups":    group_res,
        "fe_lig":    fe_lig_res,
        "his_hn":    his_hn_res,
        "u_rms":     u_rms,
        "kern_primary":   kern_primary,
        "kern_secondary": kern_secondary,
        "kern_scores":    kern_scores,
        "kern_loc":       kern_loc,
        # v3.7: Marcus-Hush-Reorganisations-Modulationen per bond Kanal.
        # Pro Mode ein Dict {parent_channel: {dr_rss_a, lambda_pair_cm1,
        # lambda_mode_cm1, n_subchannels, dr_sum_signed_a}}.
        # Aggregation over Modes occurs in the runner (compute_total_reorg
        # and compute_modulation_spectra).
        "reorg_per_mode": reorg_per_mode,
        "reorg_subchannels": reorg_subchannels,
        "pts_ref":   pts_ref_cl,
        "pts_dist":  pts_dist_cl,
        # Interne Arrays for SS and Calpha (werden after SS-Analyse entfernt)
        "_centers":  centers,
        "_evg":      evg,
        "_c2l":      c2l,
    }


def analyze_mode_with_fallback(candidates:  List[BlockInfo],
                                 col:         int,
                                 filepath:    str,
                                 atoms:       List[Dict],
                                 idx_map:     Dict[int,int],
                                 normal:      np.ndarray,
                                 coord_info:  CoordInfo,
                                 fe_c:        List[int],
                                 s_c:         List[int],
                                 cfg:         Config,
                                 mode_num:    int = None,
                                 ) -> Tuple[Optional[Dict], str]:
    """Analysiert a Mode with Fallback on alle Kandidaten-blocks.

    Parameters
    ----------
    candidates : list of BlockInfo
        Alle verfuegbaren blocks for diese Modenummer.
    col : int
        Spalten-Index in ``candidates[0]`` (dem besten Block).
    filepath : str
        Path to the Gaussian ``.log``-file.
    atoms, idx_map, normal, coord_info, fe_c, s_c, cfg :
        Wie in ``analyze_mode``.
    mode_num : int, optional
        Globale Modenummer for korrekte Spaltenwahl in jedem Kandidaten.

    Returns
    -------
    result : dict or None
        Analyseergebnis, or ``None`` for Error.
    reason : str
        Errorgrund (leer for Erfolg).

    Notes
    -----
    Als einzige validation is geprueft ob alle Eigenvektor-elements
    finit are (kein NaN, no Inf). A Amplitudenbegrenzung gibt es
    nicht: if Gaussian a Rechnung zurueckgibt, is the Ergebnis
    ausgiven werden.
    """
    last_reason = "keine Kandidaten"

    # HP-blocks haben immer Vorrang.
    # Standard-blocks are only used if for diese Mode
    # ueberhaupt no HP-Block present is (reiner Fallback).
    hp_cands  = [bi for bi in candidates if bi.is_hp]
    std_cands = [bi for bi in candidates if not bi.is_hp]
    ordered   = hp_cands if hp_cands else std_cands

    for bi in ordered:
        try:
            local_col = col
            if mode_num is not None and mode_num in bi.mode_nums:
                local_col = bi.mode_nums.index(mode_num)

            r = analyze_mode(bi, local_col, filepath, atoms, idx_map,
                              normal, coord_info, fe_c, s_c, cfg)

            if r is None:
                last_reason = (f"{'HP' if bi.is_hp else 'Std'}-Block: "
                               f"evg leer (data_offset={bi.data_offset})")
                continue

            evg = r.get("_evg")
            if evg is None or evg.shape[0] == 0:
                last_reason = "evg leer after analyze_mode"
                continue
            if not np.all(np.isfinite(evg)):
                last_reason = "evg nicht-finit"
                continue

            return r, ""

        except Exception as exc:
            last_reason = f"{'HP' if bi.is_hp else 'Std'}-Block Exception: {exc}"

    return None, last_reason


# ===========================================================================
# SCSD — Symmetry-Coordinate Structural Decomposition (rigoros, orthogonal)
# ===========================================================================
#
# **Methode**: Symmetry-Coordinate Structural Decomposition nach
# Kingsbury & Senge, *Chem. Sci.* **15**, 13638 (2024).
# Web-Implementierung: https://www.kingsbury.id.au/scsd
# Python-Bibliothek: ``scsdpy`` (https://pypi.org/project/scsdpy/)
#
# **Prinzip**: Die Verschiebung (dist - ref) the vier Clusteratome wird
# in a *orthogonale* Basis from D2h-Symmetriekoordinaten zerplaces:
#   Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u
# Jeder Irrep-contribution ``SCSD_d<Irr>`` is the Projektion on eine
# einzige Symmetriekoordinate. Die contributions are orthogonal, summieren
# quadratisch to Gesamtdeformation and are direkt zwischen
# differenten Strukturen vergleichbar — *im Gegensatz* to den
# heuristischen Scores from :func:`classify_kernel_mode_from_evg`.
#
# **Achsenkonvention** (identical with the Heuristik-Klassifikation):
#   x = Fe-Fe-Achse, y = S-S-Achse, z = cluster normal
#
# **Kanonische Referenzgeometrie**: Damit Ergebnisse differenter
# Rechnungen direkt vergleichbar sind, is IMMER dieselbe feste
# Referenzgeometrie used — *unabhaengig* vom tatsaechlich
# computethe cluster. Dies is the zentrale Konvention the Kingsbury-
# Methode and macht beispielsweise the Vergleich WT-protoniert vs.
# WT-deprotoniert vs. Mut-protoniert unmittelbar interpretierbar.
#
# Kanonische [2Fe-2S]-D2h-Referenzwerte (Meane from Rieske/Ferredoxin-
# Kristallstrukturen, Literatur):
#   Fe-Fe = 2.73 A  →  d_fe = 1.365 A
#   Fe-S  = 2.20 A  →  d_s  = sqrt(2.20^2 - 1.365^2) ~ 1.726 A
#
# Modell-Koordinaten (D2h-symmetrisch, Rhombus in the xy-Ebene):
#   Fe1 = (-d_fe, 0, 0)     Fe2 = (+d_fe, 0, 0)
#   S1  = (0, +d_s, 0)      S2  = (0, -d_s, 0)
#
# **Erweiterung on andere Referenzgeometrien**: If spaeter eine
# zweite SCSD-Analyse with a system-specificen Referenz (z. B. der
# tatsaechlich optimierten Geometrie a bestimmten Modells)
# wuenschenswert ist, kann ``_get_scsd_model`` with anderen
# Modell-Koordinaten parametrisiert werden. Solange beide Analysen
# parallel laufen, bleibt the kanonische Kingsbury-Referenz die
# Grundlage for the System-uebergreifenden Vergleich.

#: Kanonische Fe-Fe-Halbachse [A] for the SCSD-Referenzmodell.
_SCSD_D_FE: float = 1.365   # Fe-Fe = 2.73 A

#: Kanonische S-displacement [A] for the SCSD-Referenzmodell.
#: Abgeleitet aus: d_s = sqrt(Fe-S^2 − d_fe^2) = sqrt(2.20^2 − 1.365^2)
_SCSD_D_S:  float = float(np.sqrt(max(2.20**2 - _SCSD_D_FE**2, 0.1)))

#: Kanonische Modell-Koordinaten (D2h), Form (4, 3).
_SCSD_MODEL_COORDS: np.ndarray = np.array([
    [-_SCSD_D_FE,  0.,       0.],   # Fe1
    [ _SCSD_D_FE,  0.,       0.],   # Fe2
    [  0.,        _SCSD_D_S, 0.],   # S1
    [  0.,       -_SCSD_D_S, 0.],   # S2
], dtype=float)


def _get_scsd_model(_dist_ref_unused: Dict = None) -> Optional[object]:
    """Creates the D2h-Referenzmodell for the SCSD-Zerlegung.

    Parameters
    ----------
    _dist_ref_unused : dict, optional
        Wird ignoriert.  Der Parameter existiert only for API-Kompatibilitaet
        with aelteren Aufrufen.  Das Modell used stets the kanonische
        Referenzgeometrie (see Modulkopf).

    Returns
    -------
    object or None
        scsdpy-Modell-Objekt, or ``None`` if scsdpy not installiert
        or the Modell-Erstellung fehlgeschlagen ist.

    Notes
    -----
    Bugfix B1 + Hardening: Die Referenzgeometrie is jetzt kanonisch und
    fest (Fe-Fe = 2.73 A, Fe-S = 2.20 A, D2h-Rhombus).  Ergebnisse
    differenter Strukturen are so that direkt vergleichbar, analog zum
    Vorgehen the SCSD-Website (Kingsbury & Senge, Chem. Sci. 2024, 15, 13638).
    """
    try:
        from scsd.scsd import scsd_model as _scsd_model_cls
    except ImportError:
        return None
    except Exception as e:
        import warnings; warnings.warn(f"SCSD-Import: {e}", UserWarning)
        return None

    try:
        return _scsd_model_cls("2Fe2S_canonical", _SCSD_MODEL_COORDS, "D2h")
    except Exception as e:
        # Fallback: Modell without Namen versuchen
        try:
            return _scsd_model_cls(_SCSD_MODEL_COORDS, "D2h")
        except Exception:
            import warnings; warnings.warn(f"SCSD model: {e}", UserWarning)
            return None


def _parse_scsd_result(mat: List) -> Dict[str, float]:
    """
    Bugfix B1: Unterstuetzt numpy structured array UND Listen-Format.
    """
    result: Dict[str, float] = {}
    if mat is None:
        return result

    rows = mat[:-1] if hasattr(mat, "__len__") and len(mat) > 0 else []

    for row in rows:
        try:
            # Numpy structured array (neue scsdpy)
            if hasattr(row, "dtype") and hasattr(row.dtype, "names") \
                    and row.dtype.names:
                cols = row.dtype.names
                k = str(row[cols[0]])
                v = float(row[cols[1]])
            else:
                # Liste / Tuple (aeltere scsdpy)
                k = str(row[0])
                v = float(row[1])
            result[k] = v
        except (IndexError, KeyError, TypeError, ValueError):
            continue
    return result


def _run_scsd(coords_4x3: np.ndarray, model: Optional[object]) -> Dict[str, float]:
    """
    Performs SCSD-Analyse durch.
    Versucht mehrere API-Formate for Kompatibilitaet with differenten scsdpy-Versionen.
    Bein the first Error is a einmalige Diagnose ausgiven.
    """
    try:
        from scsd.scsd import scsd_matrix

        ats4 = [("Fe", tuple(coords_4x3[0])),
                ("Fe", tuple(coords_4x3[1])),
                ("S",  tuple(coords_4x3[2])),
                ("S",  tuple(coords_4x3[3]))]

        # Verschiedene Initialisierungsformate versuchen
        sm = None
        last_err = None
        for args, kwargs in [
            ([ats4],                    {"model": model}),   # Standard
            ([ats4, model],             {}),                  # Ohne Keyword
            ([np.array(coords_4x3)],   {"model": model}),   # Nur Koordinaten
        ]:
            try:
                sm = scsd_matrix(*args, **kwargs); break
            except Exception as e:
                last_err = e; continue

        if sm is None:
            raise RuntimeError(f"scsd_matrix initialization failed: {last_err}")

        # calc_scsd with and without arguments versuchen
        calc_ok = False
        for kwargs in [{"bhopping": False}, {}, {"by_graph": True}]:
            try:
                sm.calc_scsd(**kwargs); calc_ok = True; break
            except Exception as e:
                last_err = e; continue

        if not calc_ok:
            raise RuntimeError(f"calc_scsd failed: {last_err}")

        # Ergebnis from scsd_matrix or simple_scsd lesen
        for attr in ["scsd_matrix", "simple_scsd", "result"]:
            mat = getattr(sm, attr, None)
            if mat is not None:
                result = _parse_scsd_result(mat)
                if result:
                    return result

        # Einmalige Diagnostik if nichts funktioniert hat
        if not getattr(_run_scsd, "_diag_done", False):
            _run_scsd._diag_done = True
            attrs = [a for a in dir(sm) if not a.startswith("_")]
            import warnings; warnings.warn(f"SCSD-Diagnose: verfuegbare Attribute: {attrs}", UserWarning)
        return {}

    except Exception as e:
        err = str(e)
        if not getattr(_run_scsd, "_shown", False):
            _run_scsd._shown = True
            import warnings; warnings.warn(f"SCSD Error: {err[:200]}", UserWarning)
            import warnings; warnings.warn("SCSD Note: pip install --upgrade scsdpy", UserWarning)
        return {}


def _cluster_geometry(pts: np.ndarray) -> Dict[str, float]:
    """Computes alle sechs Clusterabstaende (Fe-Fe, Fe-S, S-S) in Angstrom.

    Erwartet pts als (4,3)-Array: [Fe1, Fe2, S1, S2].
    Returns leeres Dict zurueck if weniger als 4 Atome present.
    """
    if pts.shape[0] < 4: return {}
    fe1, fe2, s1, s2 = pts[0], pts[1], pts[2], pts[3]
    return {
        "fe_fe":   float(np.linalg.norm(fe2-fe1)),
        "fe1_s1":  float(np.linalg.norm(s1-fe1)),
        "fe2_s1":  float(np.linalg.norm(s1-fe2)),
        "fe1_s2":  float(np.linalg.norm(s2-fe1)),
        "fe2_s2":  float(np.linalg.norm(s2-fe2)),
        "s_s":     float(np.linalg.norm(s2-s1)),
    }


def compute_scsd_for_mode_full(pts_ref:  np.ndarray,
                                 pts_dist: np.ndarray,
                                 model,
                                 u_rms:   float = 1.0,
                                 sigma_coord:  float = 1e-3,
                                 sigma_eigvec: float = 5e-4,
                                 ) -> Dict:
    r"""Vollstaendige SCSD-Symmetriezerlegung a Mode.

    Computes the orthogonale D2h-Symmetriekoordinaten-Zerlegung
    sowohl the Gleichgewichtsgeometrie (``pts_ref``) als also der
    ausgelenkten Geometrie (``pts_dist = pts_ref + Eigenvektor·u_rms``)
    against the kanonische Kingsbury-Referenz (see Modulkopf).
    Returns per Irrep drei values:

      * ``SCSD_<Irr>_ref``  — Projektion the Gleichgewichtsgeometrie
      * ``SCSD_<Irr>_dist`` — Projektion the ausgelenkten Geometrie
      * ``SCSD_d<Irr>``     — Differenz (= Symmetrieanteil the mode)

    Zusaetzlich are zurueckgiven:

      * ``scsd_primary`` / ``scsd_secondary`` — the zwei dominanten
        Irreps the *Verschiebung* (sortiert after |SCSD_d<Irr>|).
        Ergaenzt the Heuristik-Klassifikation aus
        :func:`classify_kernel_mode_from_evg`.
      * ``total_ref`` / ``total_dist`` / ``total_d`` —
        Summenregel over alle Irreps (rms).
      * bond lengthn-changeen (``d_fe_fe``, ``d_fe1_s1``, ...)
        and atomare displacementen (``disp_Fe1`` ...).

    Parameters
    ----------
    pts_ref : ndarray of shape (4, 3)
        Gleichgewichtskoordinaten the vier Clusteratome (Fe1, Fe2, S1, S2)
        in Angstrom. Reihenfolge muss to Modell-Reihenfolge passen
        (see ``_SCSD_MODEL_COORDS``: Fe1, Fe2, S1, S2).
    pts_dist : ndarray of shape (4, 3)
        Ausgelenkte Koordinaten = pts_ref + Eigenvector * u_rms.
    model : object
        scsdpy-Modell-Objekt from :func:`_get_scsd_model`. Returns the
        kanonische D2h-Referenz and the Symmetrieoperatoren.
    u_rms : float, optional
        Thermal scaling for error propagation. Standard: ``1.0``.

    Returns
    -------
    dict
        Vollstaendiger SCSD-Output (see oben).

    Notes
    -----
    **Methode**: Kingsbury & Senge, *Chem. Sci.* **15**, 13638 (2024);
    Implementierung via ``scsdpy``. **Achsenkonvention**: x = Fe-Fe,
    y = S-S, z = cluster normal.

    Die values ``SCSD_d<Irr>`` are the *physikalisch interessante*
    Groesse — sie quantifizieren, in welche Symmetrieanteile the cluster
    for dieser Mode auseinanderbricht. ``SCSD_<Irr>_ref`` quantifiziert
    in welchem Mass the Gleichgewichtsstruktur selbst bereits von der
    kanonischen D2h-Referenz abweicht.

    Im Gegensatz to Heuristik in :func:`classify_kernel_mode_from_evg`
    are the SCSD-values ORTHOGONAL — saubere D2h-Vibrationen erscheinen
    in genau a Irrep.

    See Also
    --------
    extract_dominant_scsd_irreps : Reads dominante Irreps from the Output.
    classify_kernel_mode_from_evg : Heuristische Alternative.
    """
    out: Dict = {}
    # v1.0.4 bugfix: previously hardcoded as 5e-7 (geometry) and 5e-6
    # (eigenvector), which silently overrode the configured cfg.sigma_coord
    # and cfg.sigma_eigvec for SCSD sigmas only. The hardcoded values are
    # about 1000x smaller than the cfg defaults, so SCSD sigmas were
    # systematically under-reported. Defaults here match Config (1e-3 and
    # 5e-4) so the corrected behaviour is consistent with the rest of the
    # code. To recover pre-v1.0.4 behaviour exactly, callers can still
    # pass sigma_coord=5e-7, sigma_eigvec=5e-6.
    s_geo_ref  = float(np.sqrt(2) * sigma_coord)
    s_geo_dist = float(np.sqrt(2) *
                        np.sqrt(sigma_coord**2 + (sigma_eigvec * u_rms)**2))

    modes_ref  = _run_scsd(pts_ref,  model)
    modes_dist = _run_scsd(pts_dist, model)
    irreps     = sorted(modes_ref.keys())
    out["scsd_irreps"] = irreps

    for irr in irreps:
        vr = modes_ref.get(irr, 0.)
        vd = modes_dist.get(irr, 0.)
        out[f"SCSD_{irr}_ref"]    = vr
        out[f"s_SCSD_{irr}_ref"]  = s_geo_ref
        out[f"SCSD_{irr}_dist"]   = vd
        out[f"s_SCSD_{irr}_dist"] = s_geo_dist
        out[f"SCSD_d{irr}"]       = vd - vr
        out[f"s_SCSD_d{irr}"]     = float(np.sqrt(s_geo_ref**2 + s_geo_dist**2))

    geo_ref  = _cluster_geometry(pts_ref)
    geo_dist = _cluster_geometry(pts_dist)
    for k, vr in geo_ref.items():
        vd = geo_dist.get(k, vr)
        out[f"{k}_ref"]  = vr
        out[f"{k}_dist"] = vd
        out[f"d_{k}"]    = vd - vr
        out[f"s_d_{k}"]  = s_geo_dist

    out["total_ref"]  = float(np.sqrt(sum(v**2 for v in modes_ref.values())))
    out["total_dist"] = float(np.sqrt(sum(v**2 for v in modes_dist.values())))
    out["total_d"]    = out["total_dist"] - out["total_ref"]

    for i, lbl in enumerate(["Fe1","Fe2","S1","S2"]):
        if i < pts_ref.shape[0] and i < pts_dist.shape[0]:
            out[f"disp_{lbl}"] = float(np.linalg.norm(pts_dist[i]-pts_ref[i]))
    out["com_disp_cluster"] = float(np.linalg.norm(
        (pts_dist - pts_ref).mean(0)))

    # Rigorose Irrep-Klassifikation from SCSD: dominante D2h-Irreps der
    # Verschiebung (dist - ref). Is in the SCSD-Sheet als 'Kern-Modus
    # (SCSD)' angezeigt; ergaenzt the Heuristik-Klassifikation aus
    # ``classify_kernel_mode_from_evg``.
    primary_scsd, secondary_scsd, _ = extract_dominant_scsd_irreps(out)
    out["scsd_primary"]   = primary_scsd
    out["scsd_secondary"] = secondary_scsd

    return out


def extract_dominant_scsd_irreps(scsd_dict: Dict) -> Tuple[str, str, Dict[str, float]]:
    r"""Determines the dominanten D2h-Irreps from a SCSD-Output.

    Sortiert the :math:`|d\mathrm{SCSD}_\mathrm{Irr}|`-values (Verschiebung
    in Irrep-Koordinaten against the kanonische Kingsbury-Referenz) und
    gibt the zwei groessten contributions . Im Gegensatz zur
    Heuristik-Klassifikation in :func:`classify_kernel_mode_from_evg`
    are diese Labels rigoros — sie kommen from the orthogonalen
    Symmetriezerlegung the ``scsdpy``-Routine (Kingsbury & Senge,
    *Chem. Sci.* **15**, 13638 (2024)).

    Parameters
    ----------
    scsd_dict : dict
        Output von :func:`compute_scsd_for_mode_full`. Erwartet wird
        ``scsd_irreps`` (List of gerechneten Irrep-Labels) und
        ``SCSD_d<irr>``-Felder per Irrep.

    Returns
    -------
    primary : str
        Dominantes Irrep-Label (z. B. ``"Ag"``, ``"B1g"``).
        ``"n/a"`` if no Irreps or alle |d| < 1e-12.
    secondary : str
        Zweitgroesstes Irrep-Label, or ``"-"`` if not defined.
    contributions : dict
        ``{irrep: |d|}`` aller computeden Irreps, sortiert absteigend.
    """
    irreps = scsd_dict.get("scsd_irreps", [])
    if not irreps:
        return "n/a", "-", {}

    abs_d = {}
    for irr in irreps:
        v = scsd_dict.get(f"SCSD_d{irr}", 0.0)
        try:
            abs_d[irr] = abs(float(v))
        except (TypeError, ValueError):
            abs_d[irr] = 0.0

    # Nach |d| absteigend sortieren
    sorted_d = sorted(abs_d.items(), key=lambda kv: -kv[1])
    contributions = {k: v for k, v in sorted_d}

    if not sorted_d or sorted_d[0][1] < 1e-12:
        return "n/a", "-", contributions
    primary = sorted_d[0][0]
    secondary = sorted_d[1][0] if len(sorted_d) > 1 and sorted_d[1][1] >= 1e-12 else "-"
    return primary, secondary, contributions


__version__ = "1.4"  # modenanalyse v1.4
