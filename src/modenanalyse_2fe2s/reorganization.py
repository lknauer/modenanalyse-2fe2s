# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

# -*- coding: utf-8 -*-
"""
reorganization.py
=================

Mode-aufgeloeste vibratorische bond reorganization energyn fuer
[2Fe-2S]-Cluster (v3.7).

Wissenschaftlicher Hintergrund: Was wir berechnen
-------------------------------------------------
We compute per vibrational mode i and per bond Kanal X den
contribution

    lambda_X(i) = (1/2) mu_i omega_i^2 (Delta r_X(i))^2

mit X in {FeFe, FeN, FeS, NH, HA}. Summed over alle Modes:

    Lambda_X^total = sum_i lambda_X(i)

ist the **thermal bond-length fluctuation energy** der
givenen bond X for the givenen Temperature (bzw. die
zero-point-energy contribution for T -> 0).

Diese Groesse hat dieselbe Form wie the klassische Marcus-Hush-
reorganization energy

    lambda_Marcus = (1/2) sum_i mu_i omega_i^2 (Delta Q_i)^2

aber a andere physikalische Bedeutung. In Marcus-Hush is Delta Q_i
die *statische* Verschiebung the mode-Koordinate between Reaktant-
und Produkt-Zustand. Hier is Delta r_X(i) the *thermisch* scalede
bond modulation from the Normal-Modes (e_atoms = e_normalized * u_rms,
mit u_rms from the QHO-Erwartungswert).

Was unser Lambda_X NICHT ist
-----------------------------
* Kein direkter Marcus-Theorie-value for a spezifische Reaktion --
  for this braeuchte man Reaktant- and Produkt-Geometrien getrennt
  (zwei Frequenzrechnungen, Mode-Mapping).
* Keine Aktivierungs-energy for einen ET- or PT-Prozess.
* Keine kinetische Marcus-Rate-Vorhersage.

Was unser Lambda_X IST
----------------------
* A wohldefinede spektroskopische Observable: the energy der
  thermischen (bzw. ZPE-) Fluktuation the bond X, dekonstruiert in
  contributions the einzelnen vibrational moden.
* Geeignet for NRVS-Bandenidentifikation: frequency rangee with grossen
  Lambda_X-Beicontribute are the Banden, in denen bond X stark
  moduliert wird.
* Geeignet for System-Vergleiche for gleicher Temperature (WT vs. Mut,
  prot vs. deprot): the values skalieren konsistent.
* Verwandt with the Huang-Rhys-Faktor and the Franck-Condon-Theorie
  vibratorischer Spektren.

Pro-Mode-contribution
----------------
    Delta r_X(i) = (e_a(i) - e_b(i)) . axis(r_a, r_b)         [klassisch]
    Delta r_HA(i) = (e_H(i) - e_N(i)) . axis(r_N, r_Akz)      [Reaktionskoord.]

Konvention: e_atoms are thermisch scaled (multipliziert with u_rms
aus core.compute_thermal_amplitude). So is Delta r_X eine
modulation in Angstrom with the QHO-Wahrscheinlichkeitsverteilung als
Skala.

Reduzierte Masse: zwei Konventionen
-----------------------------------
* mu_mode (aus Gaussian/ORCA): mode-specific, garantiert orthogonale
  Eigenmodes - the thermodynamisch konsistente Wahl for the *Summe*
  Lambda_X_total. For T -> 0 is Lambda_X(i) ~ (hbar*omega/4) * alpha^2,
  wobei alpha the mode-Fraction on the bond ist.

* mu_pair (bond-specific, mu = m_a*m_b/(m_a+m_b)): isolierter-
  Oszillator-Naeherung - vergleichbar between Modes derselben bond,
  but not without Doppelzaehlung addierbar.

Im Output are beide present (lambda_pair_cm1 and lambda_mode_cm1).
Fuer system aggregation is mu_mode used.

Modulations-Spektren M_X(omega)
-------------------------------
Frequenz-aufgeloestes Modulations-Spektrum per Kanal X:

    M_X(omega) = sum_i |dr_X(i)| * G(omega - omega_i, sigma)

mit Gauss-Verbreiterung G(x, sigma) = exp(-x^2 / (2 sigma^2)).
sigma is by default 5 cm^-1 and per TOML einstellbar.

Kumulative Reorg-Spektren Lambda_X(omega)
-----------------------------------------
    Lambda_X(omega_cut) = sum_i [omega_i <= omega_cut] * lambda_X(i)

Zeigt, wo entlang the Frequenzachse the bond fluctuation energy
akkumuliert wird. Direkter Vergleich with NRVS-Spektren possible.

Was dieses Modul NICHT mehr macht (im Unterschied to v3.6)
----------------------------------------------------------
* Keine P_CPET/P_PT/P_ET-Scores
* Keine Klassifikations-Schwellen (stark/moderat)
* Keine primary_class, no Kategorisierung
* Keine Lokalisations-Filter (kern_loc, lig_loc)
* Keine tanh-Saettigung with lambda0
* Keine TOML-Parameter for Schwellen or Filter

Oeffentliche functions
-----------------------
compute_mode_modulations
    Pro Mode: dr_X and lambda_X for alle channels.
compute_total_reorganization
    Aggregates lambda_X over alle Modes.
compute_modulation_spectra
    Computes M_X(omega) als gauss-verbreiterte Spektren.
compute_cumulative_reorganization
    Computes Lambda_X(omega) als kumulative Reorg-Kurven.
"""

from __future__ import annotations
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ===========================================================================
# Physikalische Konstanten
# ===========================================================================

_AMU_KG: float = 1.66053906660e-27
_C_CMS: float  = 2.99792458e10
_H_JS: float   = 6.62607015e-34
_HC_JCM: float = _H_JS * _C_CMS

_ATOMIC_MASSES_AMU: Dict[str, float] = {
    "H":  1.00794,
    "C":  12.0107,
    "N":  14.0067,
    "O":  15.9994,
    "F":  18.9984,
    "S":  32.065,
    "Cl": 35.453,
    "Fe": 55.845,
}

_ATOMIC_NUM_TO_ELEM: Dict[int, str] = {
    1: "H",  6: "C",  7: "N",  8: "O",  9: "F",
    16: "S",  17: "Cl",  26: "Fe",
}


# ===========================================================================
# bond channels
# ===========================================================================

CHANNELS: Tuple[str, ...] = ("FeFe", "FeN", "FeS", "NH", "HA")
"""Die fuenf bond channels for reorganization energyn.

* FeFe : Fe-Fe-Cluster-Atmung (ET-Reorganisation)
* FeN  : Fe-N(His)-bond (PT-relevant in Rieske-Systemen)
* FeS  : Fe-S(Cys)-bond (allgemein in [2Fe-2S]-Clustern)
* NH   : N-H-Stretch the protonierten His
* HA   : H...acceptor-Distanz (hydrogenbruecke)

NH and HA are for deprotoniertem System on Null gesetzt, because
keine N-H-bond existiert.
"""


# ===========================================================================
# Helper: bond Modulationen entlang a Achse
# ===========================================================================

def signed_dr_along_axis(
    e_a: np.ndarray, e_b: np.ndarray,
    r_a: np.ndarray, r_b: np.ndarray,
) -> float:
    """Computes the signed bond-length modulation.

        Delta r = (e_a - e_b) . r_hat

    Konvention: r_hat zeigt von b after a. Positives Delta r heisst
    "bond dehnt sich"; negatives "bond staucht sich".

    Parameters
    ----------
    e_a, e_b : np.ndarray, shape (3,)
        displacement vectors the beiden Atome (in A, thermisch scaled).
    r_a, r_b : np.ndarray, shape (3,)
        Gleichgewichts-Positionen (in A).

    Returns
    -------
    dr_signed : float
        bond lengthnaenderung in A.
    """
    bond = r_a - r_b
    bond_len = float(np.linalg.norm(bond))
    if bond_len < 1e-10:
        return 0.0
    axis = bond / bond_len
    return float(np.dot(e_a - e_b, axis))


# ===========================================================================
# Helper: reduzierte Masse for ein Atompaar
# ===========================================================================

def reduced_mass_pair(elem_a: str, elem_b: str,
                       fallback_amu: float = 12.0107) -> float:
    """Reduced mass mu = m_a * m_b / (m_a + m_b) in amu.

    For unbekanntem element is on C-Masse zurueckgegriffen, mit
    Warnung. Das vermeidet the Score-Killer from v3.6.0.

    Parameters
    ----------
    elem_a, elem_b : str
        element-Symbole.
    fallback_amu : float
        Fallback-Masse for unbekannte elements (Default: C = 12.011).

    Returns
    -------
    mu_amu : float
    """
    ma = _ATOMIC_MASSES_AMU.get(elem_a)
    mb = _ATOMIC_MASSES_AMU.get(elem_b)
    if ma is None or mb is None:
        warnings.warn(
            f"reduced_mass_pair: unbekanntes element ({elem_a!r} or "
            f"{elem_b!r}). Fallback on {fallback_amu} amu.",
            UserWarning, stacklevel=2)
        ma = ma if ma is not None else fallback_amu
        mb = mb if mb is not None else fallback_amu
    return float(ma * mb / (ma + mb))


# ===========================================================================
# Helper: Marcus-Hush lambda per Mode
# ===========================================================================

def lambda_pair_cm1(dr_a: float, omega_cm1: float, mu_amu: float) -> float:
    """Computes the Reorg-contribution a Mode entlang a bond in cm^-1.

    Formel:
        lambda = (1/2) mu omega^2 (dr)^2

    For thermisch scaleden dr-valuesn (e_atoms = e_normalized * u_rms,
    default convention in the Aufrufer compute_mode_modulations) is das
    the QHO-Erwartungswert the bond length-fluctuation energy
    per Mode. Aufsummiert summed over all modes gives the thermal
    bond reorganization energy the bond -- eine
    spektroskopische Observable, the with the klassischen Marcus-Hush-
    reorganization energy verwandt but not identical ist
    (see module docstring).

    Parameters
    ----------
    dr_a : float
        bond modulation in A (Vorzeichen egal, only Quadratbetrag zaehlt).
    omega_cm1 : float
        Mode-Frequency in cm^-1.
    mu_amu : float
        Reduced mass in amu (mode-specific or bond-specific).

    Returns
    -------
    lambda_cm1 : float
        Reorg-contribution dieser Mode entlang the bond in cm^-1.
    """
    if not (np.isfinite(dr_a) and np.isfinite(omega_cm1)
            and np.isfinite(mu_amu)):
        return float("nan")
    if mu_amu <= 0 or omega_cm1 <= 0:
        return 0.0
    omega_si = 2.0 * np.pi * _C_CMS * omega_cm1
    mu_si = mu_amu * _AMU_KG
    dq_si = dr_a * 1.0e-10
    e_joule = 0.5 * mu_si * omega_si**2 * dq_si**2
    return float(e_joule / _HC_JCM)


# ===========================================================================
# Pro-Mode-Berechnung: dr_X and lambda_X for alle channels
# ===========================================================================

@dataclass
class ChannelGeometry:
    """Geometrie-Metadaten for einen bond Kanal.

    Saves the Atom-Indizes and element-Symbole, so that per Mode
    the dr_X and lambda_X effizient computed are koennen.

    Attributes
    ----------
    name : str
        Kanal-Name (e.g. "FeFe", "FeN_His255", "FeS_Cys207", "NH_His255",
        "HA_His255_O6265").
    idx_a, idx_b : int
        Indices the beiden Atome in the atoms_h-Array (h=mit hydrogen,
        konsistent with the H-Eigenvektor).
    r_a, r_b : np.ndarray, shape (3,)
        Gleichgewichts-Positionen the beiden Atome (in A).
    elem_a, elem_b : str
        element-Symbole.
    mu_pair_amu : float
        bondspezifische reduzierte Masse (vorcomputed).
    parent_channel : str
        Aggregat-Kategorie ("FeFe", "FeN", "FeS", "NH", "HA").
        Mehrere ChannelGeometry koennen dasselbe parent_channel haben
        (e.g. FeN_His255 and FeN_His259 -> beide below "FeN").
    weight : float
        Optionaler Gewichtungsfaktor (e.g. acceptor-Gauss-Gewichtung
        over the Gleichgewichtsabstand). default 1.0.
    idx_donor : int or None
        OPTIONAL: donor-Atom-Index for "reaction coordinaten-Modus".
        If gesetzt, is dr not als (e_a - e_b) . axis(r_a, r_b)
        computed, but rather als (e_a - e_donor) . axis(r_b, r_donor).
        Das is the echte PT-reaction coordinate for HA-channels:
        "motion of the H relativ to donor-N, projiziert auf
        N->acceptor-Achse". Misst spezifisch the PT-Komponente, nicht
        the relative acceptor-motion. None: klassische bond 
        Modulation.
    r_donor : np.ndarray, shape (3,) or None
        Position of the donor-Atoms (nur relevant if idx_donor != None).
    """
    name: str
    idx_a: int
    idx_b: int
    r_a: np.ndarray
    r_b: np.ndarray
    elem_a: str
    elem_b: str
    mu_pair_amu: float
    parent_channel: str
    weight: float = 1.0
    idx_donor: Optional[int] = None
    r_donor: Optional[np.ndarray] = None


@dataclass
class ChannelResult:
    """Pro Mode per Kanal: dr and lambda."""
    name: str
    parent_channel: str
    dr_signed_a: float
    lambda_pair_cm1: float
    lambda_mode_cm1: float
    weight: float = 1.0


def compute_mode_modulations(
    e_atoms: np.ndarray,
    omega_cm1: float,
    mode_red_mass_amu: Optional[float],
    channels: Sequence[ChannelGeometry],
) -> List[ChannelResult]:
    """Computes dr_X and lambda_X for alle bond channels a Mode.

    Zwei Berechnungs-Modi per after ChannelGeometry.idx_donor:
    
    1. Klassische bond modulation (idx_donor None):
       dr = (e[idx_a] - e[idx_b]) . axis(r_a, r_b)
       Misst the modulation the bond length a-b.
    
    2. reaction coordinaten-Modus (idx_donor gesetzt):
       dr = (e[idx_a] - e[idx_donor]) . axis(r_b, r_donor)
       Misst the motion von a relativ to donor, projiziert auf
       the donor->b-Achse. Fuer HA-channels: H-motion relativ
       to donor-N, in N->acceptor-Richtung. Das is the echte
       PT-reaction coordinate.

    Parameters
    ----------
    e_atoms : np.ndarray, shape (N_atoms, 3)
        Thermisch scalede displacement vectors aller Atome (in A).
    omega_cm1 : float
        Mode-Frequency in cm^-1.
    mode_red_mass_amu : float or None
        Mode-eigene reduzierte Masse from Gaussian/ORCA. If None,
        is lambda_mode on NaN gesetzt; lambda_pair bleibt defined.
    channels : Sequence of ChannelGeometry
        List of bond channels.

    Returns
    -------
    results : list of ChannelResult
        Pro Kanal ein Result.
    """
    results: List[ChannelResult] = []
    n_atoms = len(e_atoms)
    for ch in channels:
        if ch.idx_a < 0 or ch.idx_b < 0 \
                or ch.idx_a >= n_atoms or ch.idx_b >= n_atoms:
            results.append(ChannelResult(
                name=ch.name, parent_channel=ch.parent_channel,
                dr_signed_a=float("nan"),
                lambda_pair_cm1=float("nan"),
                lambda_mode_cm1=float("nan"),
                weight=ch.weight))
            continue
        
        # Auswahl: reaction coordinaten-Modus or klassische Modulation
        if (ch.idx_donor is not None and ch.r_donor is not None
                and 0 <= ch.idx_donor < n_atoms):
            # reaction coordinaten-Modus:
            # dr = motion von a relativ to donor, projiziert on donor->b-Achse
            # axis_unit zeigt von donor (r_donor) after b (r_b).
            bond = ch.r_b - ch.r_donor
            bond_len = float(np.linalg.norm(bond))
            if bond_len < 1e-10:
                dr = 0.0
            else:
                axis = bond / bond_len
                dr = float(np.dot(e_atoms[ch.idx_a] - e_atoms[ch.idx_donor], axis))
        else:
            # Klassische bond modulation a-b
            dr = signed_dr_along_axis(
                e_atoms[ch.idx_a], e_atoms[ch.idx_b], ch.r_a, ch.r_b)
        
        lam_pair = lambda_pair_cm1(dr, omega_cm1, ch.mu_pair_amu)
        if mode_red_mass_amu is not None and mode_red_mass_amu > 0:
            lam_mode = lambda_pair_cm1(dr, omega_cm1, mode_red_mass_amu)
        else:
            lam_mode = float("nan")
        results.append(ChannelResult(
            name=ch.name, parent_channel=ch.parent_channel,
            dr_signed_a=dr, lambda_pair_cm1=lam_pair,
            lambda_mode_cm1=lam_mode, weight=ch.weight))
    return results


# ===========================================================================
# Aggregation: per Mode per Parent-Channel
# ===========================================================================

def aggregate_by_subchannel(
    results_per_mode: List[List[ChannelResult]],
) -> Dict[str, Dict[str, float]]:
    """Aggregates per SUB-Kanal (e.g. FeS_Cys207, FeN_His255 etc.) ueber
    alle Modes.
    
    Im Gegensatz to aggregate_by_parent, the per Mode the Sub-channels zu
    a Parent-Kanal zusammenfasst, liefert diese Funktion per Sub-
    Kanal seine eigene system total-Reorg. So kann man sehen, welche
    *einzelne* bond at the meisten to Reorganisation beicontributes.

    Parameters
    ----------
    results_per_mode : list of list of ChannelResult
        Pro Mode the List of ChannelResult from compute_mode_modulations.

    Returns
    -------
    sub_totals : dict of {channel_name: {key: float}}
        Mit keyn "lambda_total_pair_cm1", "lambda_total_mode_cm1",
        "n_modes_contributing", "parent_channel", "weight".
    """
    sub_totals: Dict[str, Dict[str, float]] = {}
    for mode_results in results_per_mode:
        for r in mode_results:
            key = r.name
            if key not in sub_totals:
                sub_totals[key] = {
                    "parent_channel":         r.parent_channel,
                    "weight":                 r.weight,
                    "lambda_total_pair_cm1":  0.0,
                    "lambda_total_mode_cm1":  0.0,
                    "n_modes_contributing":   0,
                }
            d = sub_totals[key]
            # v3.7.4: konsistente Count -- jede Mode, the zu
            # IRGENDEINEM the beiden Lambdas (pair or mode) beicontributes,
            # is gezaehlt.
            contributed = False
            if np.isfinite(r.lambda_pair_cm1) and r.lambda_pair_cm1 > 0:
                d["lambda_total_pair_cm1"] += r.weight * r.lambda_pair_cm1
                contributed = True
            if np.isfinite(r.lambda_mode_cm1) and r.lambda_mode_cm1 > 0:
                d["lambda_total_mode_cm1"] += r.weight * r.lambda_mode_cm1
                contributed = True
            if contributed:
                d["n_modes_contributing"] += 1
    return sub_totals


def aggregate_by_parent(results: List[ChannelResult]) -> Dict[str, Dict[str, float]]:
    """Aggregates Channel-Results after parent_channel.

    For mehreren Sub-channelsn (e.g. FeN_His255 + FeN_His259) wird:

    * dr_rss_a  : RSS (Root Sum of Squares), sqrt(sum_i w_i * dr_i^2),
                  vorzeichenfrei -- so that positive and negative contributions
                  not aufheben. Streng genommen a gewichtete
                  L2-Norm the sub-channel-dr-values; in the spektroskopischen
                  Sprachgebrauch oft als RSS bezeichnet.
    * lambda_X  : einfache Summe (reorg contributions addieren linear).

    Parameters
    ----------
    results : list of ChannelResult
        Pro-Kanal-Ergebnisse from compute_mode_modulations.

    Returns
    -------
    agg : dict of {parent_channel: dict}
        Mit keyn "dr_rss_a", "lambda_pair_cm1", "lambda_mode_cm1".
    """
    agg: Dict[str, Dict[str, float]] = {ch: {
        "dr_rss_a": 0.0, "dr_sum_signed_a": 0.0,
        "lambda_pair_cm1": 0.0, "lambda_mode_cm1": 0.0,
        "n_subchannels": 0,
    } for ch in CHANNELS}
    for r in results:
        pc = r.parent_channel
        if pc not in agg:
            continue
        if not (np.isfinite(r.dr_signed_a) and np.isfinite(r.lambda_pair_cm1)):
            continue
        w = r.weight
        agg[pc]["dr_rss_a"]        += w * r.dr_signed_a**2
        agg[pc]["dr_sum_signed_a"] += w * r.dr_signed_a
        agg[pc]["lambda_pair_cm1"] += w * r.lambda_pair_cm1
        if np.isfinite(r.lambda_mode_cm1):
            agg[pc]["lambda_mode_cm1"] += w * r.lambda_mode_cm1
        agg[pc]["n_subchannels"]   += 1
    for pc in agg:
        # RSS (Wurzel-aus-Quadratsumme) statt Quadratsumme
        agg[pc]["dr_rss_a"] = float(np.sqrt(agg[pc]["dr_rss_a"]))
    return agg


# ===========================================================================
# Bau the Channel list from coord_info + atoms
# ===========================================================================

def build_channels(
    coord_info,
    atoms: List[Dict],
    idx_map: Dict[int, int],
    *,
    fe1_center: int,
    fe2_center: int,
    s_centers: Sequence[int],
    his_n_centers: Sequence[int],
    his_n_labels: Sequence[str],
    his_h_centers: Sequence[int],
    cys_s_centers: Sequence[int],
    cys_s_labels: Sequence[str],
    fe_centers_per_his_n: Sequence[int],
    fe_centers_per_cys_s: Sequence[int],
    acceptors_per_h: Sequence[Sequence[int]],
    acceptor_elem_per_h: Sequence[Sequence[str]],
    eq_distances_per_h: Sequence[Sequence[float]],
    acceptor_r0_a: float = 2.8,
    acceptor_sigma_a: float = 0.4,
) -> List[ChannelGeometry]:
    """Baut the Liste aller bond channels for ein System.

    Parameters
    ----------
    coord_info, atoms, idx_map :
        Wie in core.py: List of Atom-Dicts and Center-zu listn-Mapping.
    fe1_center, fe2_center : int
        Gauss-Center-IDs the zwei Fe atoms.
    s_centers : list of int
        Gauss-Center-IDs the zwei Cluster-S-Atome.
    his_n_centers, his_n_labels :
        Pro His-ligands: N-Center-ID and Label "His 255" etc.
    his_h_centers :
        Pro His-ligands: H-Center-ID (1:1 with his_n_centers).
    cys_s_centers, cys_s_labels :
        Pro Cys-ligands: S-Center-ID and Label "Cys 207" etc.
    fe_centers_per_his_n, fe_centers_per_cys_s :
        Pro His/Cys-ligands: an welches Fe gebunden (Center-ID).
    acceptors_per_h :
        Pro His-H: List of acceptor-Center-IDs.
    acceptor_elem_per_h :
        Pro His-H: List of acceptor-elements.
    eq_distances_per_h :
        Pro His-H: List of Gleichgewichts-H...acceptor-Abstaende.
    acceptor_r0_a, acceptor_sigma_a :
        Gauss-Gewichtungs-Parameter for the H...acceptor-Distanz.
        Beiliegende acceptoren bekommen mehr Gewicht. r0 is optimaler
        bondabstand, sigma the tolerance.

    Returns
    -------
    channels : list of ChannelGeometry
    """
    channels: List[ChannelGeometry] = []
    
    def _atom_pos(center_id: int) -> Optional[np.ndarray]:
        if center_id not in idx_map:
            return None
        a = atoms[idx_map[center_id]]
        return np.array([a["x"], a["y"], a["z"]], dtype=float)
    
    def _atom_elem(center_id: int) -> str:
        if center_id not in idx_map:
            return "X"
        a = atoms[idx_map[center_id]]
        an = int(a.get("atomic_num", a.get("an", 0)))
        return _ATOMIC_NUM_TO_ELEM.get(an, "X")
    
    def _atom_listidx(center_id: int) -> int:
        return idx_map.get(center_id, -1)

    # === FeFe ============================================================
    pa = _atom_pos(fe1_center); pb = _atom_pos(fe2_center)
    if pa is not None and pb is not None:
        channels.append(ChannelGeometry(
            name="FeFe", idx_a=_atom_listidx(fe1_center),
            idx_b=_atom_listidx(fe2_center),
            r_a=pa, r_b=pb, elem_a="Fe", elem_b="Fe",
            mu_pair_amu=reduced_mass_pair("Fe", "Fe"),
            parent_channel="FeFe"))

    # === FeN (pro His) ===================================================
    for i, (n_c, label, fe_c) in enumerate(zip(
            his_n_centers, his_n_labels, fe_centers_per_his_n)):
        pa = _atom_pos(fe_c); pb = _atom_pos(n_c)
        if pa is None or pb is None:
            continue
        channels.append(ChannelGeometry(
            name=f"FeN_{label.replace(' ', '')}",
            idx_a=_atom_listidx(fe_c), idx_b=_atom_listidx(n_c),
            r_a=pa, r_b=pb, elem_a="Fe", elem_b="N",
            mu_pair_amu=reduced_mass_pair("Fe", "N"),
            parent_channel="FeN"))

    # === FeS (pro Cys) ===================================================
    for i, (s_c, label, fe_c) in enumerate(zip(
            cys_s_centers, cys_s_labels, fe_centers_per_cys_s)):
        pa = _atom_pos(fe_c); pb = _atom_pos(s_c)
        if pa is None or pb is None:
            continue
        channels.append(ChannelGeometry(
            name=f"FeS_{label.replace(' ', '')}",
            idx_a=_atom_listidx(fe_c), idx_b=_atom_listidx(s_c),
            r_a=pa, r_b=pb, elem_a="Fe", elem_b="S",
            mu_pair_amu=reduced_mass_pair("Fe", "S"),
            parent_channel="FeS"))
    
    # FeS also for the Cluster-S-Atome (Fe-S-bridging)
    for fe_c in (fe1_center, fe2_center):
        for s_c in s_centers:
            pa = _atom_pos(fe_c); pb = _atom_pos(s_c)
            if pa is None or pb is None:
                continue
            # Nur, if the beiden Atome tatsaechlich gebunden are (< 3.0 A)
            if np.linalg.norm(pa - pb) > 3.0:
                continue
            channels.append(ChannelGeometry(
                name=f"FeS_Cluster_{fe_c}_{s_c}",
                idx_a=_atom_listidx(fe_c), idx_b=_atom_listidx(s_c),
                r_a=pa, r_b=pb, elem_a="Fe", elem_b="S",
                mu_pair_amu=reduced_mass_pair("Fe", "S"),
                parent_channel="FeS"))

    # === NH and HA (pro His-H) ===========================================
    # HA-Konvention v3.7.1: per His-H is only the HAUPT-acceptor (mit
    # hoechstem Gauss-Gewicht over the Gleichgewichtsabstand) used,
    # and dr_HA is in the reaction coordinaten-Modus computed:
    #   dr_HA = (e_H - e_N_donor) . axis_N_to_Akz
    # Das misst spezifisch the H-motion relativ to donor in Richtung
    # acceptor (= echte PT-reaction coordinate), not the unspezifische
    # H...acceptor-Distanz-Modulation.
    for i, (n_c, h_c, label) in enumerate(zip(
            his_n_centers, his_h_centers, his_n_labels)):
        pa = _atom_pos(n_c); pb = _atom_pos(h_c)
        if pa is not None and pb is not None and h_c >= 0:
            # NH-Stretch: klassische bond modulation N-H
            channels.append(ChannelGeometry(
                name=f"NH_{label.replace(' ', '')}",
                idx_a=_atom_listidx(n_c), idx_b=_atom_listidx(h_c),
                r_a=pa, r_b=pb, elem_a="N", elem_b="H",
                mu_pair_amu=reduced_mass_pair("N", "H"),
                parent_channel="NH"))
            
            # HA: per H only HAUPT-acceptor wählen (höchstes Gauss-Gewicht)
            if i < len(acceptors_per_h):
                accs = acceptors_per_h[i]
                elems = acceptor_elem_per_h[i] if i < len(acceptor_elem_per_h) else []
                eqs = eq_distances_per_h[i] if i < len(eq_distances_per_h) else []
                
                # Beste Wahl: acceptor with max Gauss-Gewicht
                best_idx = -1; best_w = -1.0
                for j, a_c in enumerate(accs):
                    if _atom_pos(a_c) is None:
                        continue
                    eq = eqs[j] if j < len(eqs) else float(
                        np.linalg.norm(_atom_pos(a_c) - pb))
                    w = float(np.exp(-((eq - acceptor_r0_a) ** 2) /
                                      (2.0 * acceptor_sigma_a ** 2)))
                    if w > best_w:
                        best_w = w; best_idx = j
                
                if best_idx >= 0:
                    a_c = accs[best_idx]
                    pa_acc = _atom_pos(a_c)
                    elem = elems[best_idx] if best_idx < len(elems) else _atom_elem(a_c)
                    if elem == "X":
                        elem = _atom_elem(a_c)
                    if elem != "X":
                        # reaction coordinaten-Modus:
                        #   idx_a = H, idx_donor = N, axis N->acceptor
                        # r_a is hier only for Konsistenz gefuellt; relevant
                        # is (r_b, r_donor) wegen idx_donor != None.
                        #
                        # v3.7.4: weight=1.0. Das Gauss-Gewicht (best_w) wird
                        # AUSSCHLIESSLICH for the Hauptakzeptor-Auswahl
                        # used (binaere Logik: is a H-Bridge-acceptor
                        # in pcet_hbond_cutoff_a Reichweite?). Eine
                        # additionallye Daempfung in the Lambda-Aggregation
                        # waere doppelte Strafe -- the geometrische
                        # dr_HA-modulation is for weiten acceptoren ohnehin
                        # klein, without dass man sie nochmal with dem
                        # Gauss-Gewicht reduzieren muesste.
                        channels.append(ChannelGeometry(
                            name=f"HA_{label.replace(' ', '')}_{a_c}",
                            idx_a=_atom_listidx(h_c),
                            idx_b=_atom_listidx(a_c),
                            r_a=pb, r_b=pa_acc,
                            elem_a="H", elem_b=elem,
                            mu_pair_amu=reduced_mass_pair("H", elem),
                            parent_channel="HA",
                            weight=1.0,
                            idx_donor=_atom_listidx(n_c),
                            r_donor=pa,  # donor-N-Position
                        ))

    return channels


# ===========================================================================
# system aggregation: Total-Reorg per Kanal
# ===========================================================================

def compute_total_reorganization(
    per_mode_aggregates: List[Dict[str, Dict[str, float]]],
) -> Dict[str, Dict[str, float]]:
    """Aggregates pro-Mode-contributions to system total-Reorgs per Kanal.

    Marcus-Hush-Konvention: Lambda_total = sum_i lambda_i.
    We use lambda_mode (mode-specifice reduzierte Masse) als
    primaere Quantitaet, because diese contributions without Doppelzaehlung
    addieren. lambda_pair is ebenfalls summiert, but only als
    bond-specifice Vergleichsgroesse.

    Parameters
    ----------
    per_mode_aggregates : list of dict
        Pro Mode the aggregate-Dict from aggregate_by_parent.

    Returns
    -------
    totals : dict of {parent_channel: {key: float}}
        Mit keyn "lambda_total_pair_cm1", "lambda_total_mode_cm1",
        "n_modes_contributing" per Kanal.
    """
    totals: Dict[str, Dict[str, float]] = {ch: {
        "lambda_total_pair_cm1": 0.0,
        "lambda_total_mode_cm1": 0.0,
        "n_modes_contributing":  0,
    } for ch in CHANNELS}
    for agg in per_mode_aggregates:
        for ch in CHANNELS:
            if ch not in agg:
                continue
            d = agg[ch]
            lp = d.get("lambda_pair_cm1", 0.0)
            lm = d.get("lambda_mode_cm1", 0.0)
            contributed = False
            if np.isfinite(lp) and lp > 0:
                totals[ch]["lambda_total_pair_cm1"] += lp
                contributed = True
            if np.isfinite(lm) and lm > 0:
                totals[ch]["lambda_total_mode_cm1"] += lm
                contributed = True
            # v3.7.4: zaehle the Mode, sobald sie to IRGENDEINEM the beiden
            # Lambdas (pair or mode) beicontributes. Vorher was only bei
            # lambda_mode hochgezaehlt, was n_modes_contributing for den
            # pair-value systematisch unterzaehlte, if for einzelne
            # Modes mode_red_mass_amu not verfuegbar war.
            if contributed:
                totals[ch]["n_modes_contributing"] += 1
    return totals


# ===========================================================================
# Frequenz-aufgeloestes Modulations-Spektrum M_X(omega)
# ===========================================================================

def compute_modulation_spectra(
    freqs: np.ndarray,
    per_mode_aggregates: List[Dict[str, Dict[str, float]]],
    grid_cm1: np.ndarray,
    sigma_cm1: float = 5.0,
) -> Dict[str, np.ndarray]:
    """Computes M_X(omega) als gauss-verbreiterte Modulations-Spektren.

    Pro Kanal X:
        M_X(omega) = sum_i |dr_X(i)| * G(omega - omega_i, sigma)

    with G(x, sigma) = exp(-x^2 / (2 sigma^2)) / (sigma * sqrt(2 pi)).

    Die Gewichtung occurs over the RSS-Aggregation dr_rss (vorzeichen-
    frei), so that positive and negative Modulationen nicht
    gegenseitig ausloeschen.

    Parameters
    ----------
    freqs : np.ndarray, shape (N_modes,)
        Frequenzen aller Modes in cm^-1.
    per_mode_aggregates : list of dict
        Pro Mode the aggregate-Dict from aggregate_by_parent.
    grid_cm1 : np.ndarray, shape (N_grid,)
        Frequenz-Raster in cm^-1.
    sigma_cm1 : float
        Gauss-Verbreiterung in cm^-1.

    Returns
    -------
    spectra : dict of {parent_channel: np.ndarray, shape (N_grid,)}
        Pro Kanal the Modulations-Spektrum.
        Einheit: A * cm  (dr * Verbreitungs-Normierung)
    """
    spectra: Dict[str, np.ndarray] = {
        ch: np.zeros_like(grid_cm1, dtype=float) for ch in CHANNELS}
    if sigma_cm1 <= 0:
        return spectra
    norm = 1.0 / (sigma_cm1 * np.sqrt(2.0 * np.pi))
    inv_2sig2 = 1.0 / (2.0 * sigma_cm1 * sigma_cm1)
    
    for f, agg in zip(freqs, per_mode_aggregates):
        if not np.isfinite(f):
            continue
        # Gauss-Faktor for alle Gridpunkte
        gauss = norm * np.exp(-((grid_cm1 - f) ** 2) * inv_2sig2)
        for ch in CHANNELS:
            if ch not in agg:
                continue
            dr = agg[ch].get("dr_rss_a", 0.0)
            if np.isfinite(dr) and dr > 0:
                spectra[ch] += dr * gauss
    return spectra


def compute_co_modulation_spectra(
    spectra: Dict[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    """Computes Co-Modulations-Spektren als geometrische Mittel zweier
    M_X-Spektren.

    Diese are hilfreich, um Banden to identifizieren, in denen *gleich-
    zeitig* mehrere bonden moduliert are (PCET-Indikator).

    Definitionen:
      C_PCET     = sqrt(M_HA  * M_FeFe)  (H-A + Cluster-Atmung)
      C_PT_FeN   = sqrt(M_HA  * M_FeN)   (H-A + Fe-N(His))
      C_ET_FeS   = sqrt(M_FeFe * M_FeS)  (Cluster-Atmung + Cys-Beteiligung)

    Parameters
    ----------
    spectra : dict
        Output von compute_modulation_spectra.

    Returns
    -------
    co_spectra : dict
        Mit keyn "C_PCET", "C_PT_FeN", "C_ET_FeS".
    """
    co: Dict[str, np.ndarray] = {}
    if "HA" in spectra and "FeFe" in spectra:
        co["C_PCET"]   = np.sqrt(np.maximum(0, spectra["HA"] * spectra["FeFe"]))
    if "HA" in spectra and "FeN" in spectra:
        co["C_PT_FeN"] = np.sqrt(np.maximum(0, spectra["HA"] * spectra["FeN"]))
    if "FeFe" in spectra and "FeS" in spectra:
        co["C_ET_FeS"] = np.sqrt(np.maximum(0, spectra["FeFe"] * spectra["FeS"]))
    return co


# ===========================================================================
# Kumulative Reorg-Kurven Lambda_X(omega)
# ===========================================================================

def compute_cumulative_reorganization(
    freqs: np.ndarray,
    per_mode_aggregates: List[Dict[str, Dict[str, float]]],
    grid_cm1: np.ndarray,
    use_mode_mass: bool = True,
) -> Dict[str, np.ndarray]:
    """Computes Lambda_X(omega) = sum_{i: omega_i <= omega} lambda_X(i).

    By default used the mode-specifice reduzierte Masse
    (lambda_mode), because diese contributions without Doppelzaehlung
    addieren (Marcus-Hush-Konvention).

    Parameters
    ----------
    freqs : np.ndarray, shape (N_modes,)
        Frequenzen aller Modes.
    per_mode_aggregates : list of dict
    grid_cm1 : np.ndarray, shape (N_grid,)
        Frequenz-Raster in cm^-1.
    use_mode_mass : bool
        If True: lambda_mode (Mode-mu) - Marcus-Hush-korrekt for Summen.
        If False: lambda_pair (bond mu) - bindungs-spezifischer.

    Returns
    -------
    cum : dict of {parent_channel: np.ndarray, shape (N_grid,)}
        Kumulative Reorg-Kurve per Kanal in cm^-1.
    """
    key = "lambda_mode_cm1" if use_mode_mass else "lambda_pair_cm1"
    
    # Sortiere Modes after Frequency for effiziente Akkumulation
    order = np.argsort(freqs)
    sorted_freqs = freqs[order]
    
    cum: Dict[str, np.ndarray] = {ch: np.zeros_like(grid_cm1, dtype=float)
                                    for ch in CHANNELS}
    
    # Akkumuliere Lambda per Kanal in derselben Reihenfolge
    for ch in CHANNELS:
        running = 0.0
        # Gehe for each grid-Punkt through sortierte Modes
        # Effizient: ein Marsch through beide Listen
        mi = 0  # mode-index in sortierter Liste
        for gi, gf in enumerate(grid_cm1):
            # Akkumuliere alle Modes with freq <= gf
            while mi < len(sorted_freqs) and sorted_freqs[mi] <= gf:
                orig_idx = order[mi]
                if orig_idx < len(per_mode_aggregates):
                    val = per_mode_aggregates[orig_idx].get(ch, {}).get(key, 0.0)
                    if np.isfinite(val):
                        running += val
                mi += 1
            cum[ch][gi] = running
    return cum


__all__ = [
    "CHANNELS",
    "ChannelGeometry",
    "ChannelResult",
    "signed_dr_along_axis",
    "reduced_mass_pair",
    "lambda_pair_cm1",
    "compute_mode_modulations",
    "aggregate_by_parent",
    "aggregate_by_subchannel",
    "build_channels",
    "build_channels_from_coord_info",
    "compute_total_reorganization",
    "compute_modulation_spectra",
    "compute_co_modulation_spectra",
    "compute_cumulative_reorganization",
]


# ===========================================================================
# High-Level-Adapter: from coord_info direkt Channels bauen
# ===========================================================================

def build_channels_from_coord_info(
    coord_info,
    atoms_h: List[Dict],
    idx_map_h: Dict[int, int],
    *,
    acceptor_r0_a: float = 2.8,
    acceptor_sigma_a: float = 0.4,
) -> List[ChannelGeometry]:
    """High-Level-Adapter: extrahiert alles, was wir brauchen, direkt aus
    coord_info and seinen Sub-Strukturen (pcet_info, et_info, ligands).

    Parameters
    ----------
    coord_info : CoordInfo
        Aus geometry.find_coordinating_residues. Muss .ligands haben,
        sowie .pcet_info and .et_info (aus pcet_et.build_pcet_info /
        build_et_info).
    atoms_h : list of dict
        Voll-atom list with H-Atomen (gleiche Konvention wie der
        H-Eigenvektor).
    idx_map_h : dict {center_id: idx in atoms_h}
    acceptor_r0_a, acceptor_sigma_a :
        Gauss-Gewicht for the H...acceptor-Distanz.

    Returns
    -------
    channels : list of ChannelGeometry
    """
    pcet_info = getattr(coord_info, "pcet_info", None)
    et_info   = getattr(coord_info, "et_info", None)
    
    if et_info is None:
        return []

    fe1_c, fe2_c = et_info.fe_centers
    
    # Cluster-S-Centers: cluster_centers without Fe1 and Fe2
    s_centers = [c for c in et_info.cluster_centers
                 if c not in (fe1_c, fe2_c)]
    
    # His-ligands: from coord_info.ligands, His-element=N
    # His-N and His-H mappen from pcet_info (wenn present):
    his_n_centers: List[int] = []
    his_n_labels:  List[str] = []
    his_h_centers: List[int] = []
    fe_centers_per_his_n: List[int] = []
    
    # Konstruiere Mapping center -> ligand Info for Fe-N-bond
    n_to_lig = {}
    for lig in coord_info.ligands:
        if lig.lig_element == "N":
            n_to_lig[lig.lig_center] = lig
    
    if pcet_info is not None and pcet_info.enabled:
        # His with H: from pcet_info
        for n_c, h_c in zip(pcet_info.his_n_centers, pcet_info.his_h_centers):
            # finde Label and Fe-bond
            lig = n_to_lig.get(n_c)
            if lig is not None:
                his_n_centers.append(n_c)
                his_n_labels.append(lig.res_label)
                his_h_centers.append(h_c)
                fe_centers_per_his_n.append(lig.fe_center)
            else:
                # H-tragendes N koennte the andere Imidazol-N sein (Nε vs Nδ)
                # Suche the zugehoerige Fe-koordinierende N over his_hn_center
                for L in coord_info.ligands:
                    if L.lig_element == "N" and L.his_hn_center == n_c:
                        his_n_centers.append(n_c)
                        his_n_labels.append(L.res_label)
                        his_h_centers.append(h_c)
                        fe_centers_per_his_n.append(L.fe_center)
                        break
    else:
        # Auch for deprotonierten His: Fe-N-bond is relevant for FeN-Reorg
        for lig in coord_info.ligands:
            if lig.lig_element == "N":
                his_n_centers.append(lig.lig_center)
                his_n_labels.append(lig.res_label)
                his_h_centers.append(-1)  # no H
                fe_centers_per_his_n.append(lig.fe_center)
    
    # Cys-ligands: from coord_info.ligands, Cys-element=S
    cys_s_centers: List[int] = []
    cys_s_labels:  List[str] = []
    fe_centers_per_cys_s: List[int] = []
    for lig in coord_info.ligands:
        if lig.lig_element == "S":
            cys_s_centers.append(lig.lig_center)
            cys_s_labels.append(lig.res_label)
            fe_centers_per_cys_s.append(lig.fe_center)
    
    # acceptor listn (nur if pcet_info enabled)
    if pcet_info is not None and pcet_info.enabled:
        acceptors_per_h = list(pcet_info.acceptor_centers_per_h)
        eq_distances_per_h = list(pcet_info.eq_distances_per_h)
        # element-Symbole per acceptor: from atom_dict
        acceptor_elem_per_h: List[List[str]] = []
        for accs in acceptors_per_h:
            elems = []
            for a_c in accs:
                if a_c in idx_map_h:
                    a = atoms_h[idx_map_h[a_c]]
                    an = int(a.get("atomic_num", a.get("an", 0)))
                    elems.append(_ATOMIC_NUM_TO_ELEM.get(an, "X"))
                else:
                    elems.append("X")
            acceptor_elem_per_h.append(elems)
    else:
        acceptors_per_h = [[] for _ in his_h_centers]
        acceptor_elem_per_h = [[] for _ in his_h_centers]
        eq_distances_per_h = [[] for _ in his_h_centers]
    
    return build_channels(
        coord_info=coord_info, atoms=atoms_h, idx_map=idx_map_h,
        fe1_center=fe1_c, fe2_center=fe2_c, s_centers=s_centers,
        his_n_centers=his_n_centers, his_n_labels=his_n_labels,
        his_h_centers=his_h_centers,
        cys_s_centers=cys_s_centers, cys_s_labels=cys_s_labels,
        fe_centers_per_his_n=fe_centers_per_his_n,
        fe_centers_per_cys_s=fe_centers_per_cys_s,
        acceptors_per_h=acceptors_per_h,
        acceptor_elem_per_h=acceptor_elem_per_h,
        eq_distances_per_h=eq_distances_per_h,
        acceptor_r0_a=acceptor_r0_a, acceptor_sigma_a=acceptor_sigma_a,
    )
