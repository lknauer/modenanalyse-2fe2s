# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

# -*- coding: utf-8 -*-
"""
embedding.py
=========================
Feature-Matrix and Dimensionsreduktions-Einbettungen.

Oeffentliche functions
-----------------------
build_feature_matrix
    Creates Basis- and erweiterte Feature-Matrix (ohne H-Atome).
compute_embeddings
    Computes PCA, t-SNE and UMAP; fuehrt HDBSCAN-Clustering durch.
characterize_clusters
    Charakterisiert HDBSCAN-Cluster with Z-Scores and Fisher-F.
compute_ss_umap_cluster
    UMAP clustering of secondary-structure amplitudes.

Bugfixes (gegenvia Vorversion)
---------------------------------
B4  UMAP_Cluster-Sheet:  Feature-Spalten are jetzt befuellt.
    H-Atome and H-N-motionen fliessen NICHT in the Feature-Matrix ein.
B19 HDBSCAN:             Wird on ALLEN computeden Embeddings ausgefuehrt.

scipy-Fix
---------
``ConstantInputWarning`` von ``scipy.stats.f_oneway`` is in
``characterize_clusters`` automatically unterdrueckt (tritt auf, wenn
ein SS-Feature over alle modes konstant ist).
"""
from __future__ import annotations
import warnings
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from .config  import Config
from .logio      import RunLog
from .geometry import CoordInfo


# ===========================================================================
# Feature-Matrix
# ===========================================================================

from .core import SCORE_KEYS as _SCORE_KEYS



def build_feature_matrix(
        results:    List[Dict],
        coord_info: CoordInfo,
) -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
    """Creates Basis- and erweiterte Feature-Matrix for Embeddings.

    H-Atome and H-N-motionen (``his_hn``) are not enthalten, damit
    hydrogen not the UMAP-Topologie beeinflusst.

    Parameters
    ----------
    results : list of dict
        Modenanalyse-Ergebnisse.
    coord_info : CoordInfo
        Koordinations-Information (liefert Gruppen- and ligands-Labels).

    Returns
    -------
    X_basis : ndarray of shape (n_modes, 3)
        Basis-Features: ``[OOP%, kern_OOP%, kern_|d|]``.
    X_extended : ndarray of shape (n_modes, n_features)
        Erweiterte Features: Basis + Kern-Scores + Gruppen-OOP/INP
        + Fe-ligands-Stretch/Bend.
    feat_basis : list of str
        Feature-Namen for the Basis-Matrix.
    feat_ext : list of str
        Feature-Namen for the erweiterte Matrix.
    """
    group_names = list(coord_info.group_map.keys())
    # ligands-elements for Fe-Lig Feature-Namen
    lig_labels = list({lig.res_label for lig in coord_info.ligands})
    lig_labels.sort()
    
    # v3.7.2: Welche Gruppen are Cys? Cys-OOP-Features are geometrisch
    # immer ~0 (Cys-Sgamma liegt in the cluster plane), bringen also keine
    # Information. Wir entfernen sie from the Feature-Matrix. Identifikation
    # over the res_label-Praefix ("Cys ...") in the ligands.
    cys_labels: Set[str] = {lig.res_label for lig in coord_info.ligands
                              if lig.lig_element == "S"}

    feat_b: List[str] = ["lig_OOP%", "kern_OOP%", "kern_|d|"]
    # Per-ligands: stretch, bend_inp, bend_oop (in dieser Reihenfolge,
    # konsistent with the Schleife unten).
    _per_lig: List[str] = []
    for l in lig_labels:
        _per_lig += [f"{l}_stretch", f"{l}_bend_inp", f"{l}_bend_oop"]
    
    # OOP-Features: only for Nicht-Cys-Gruppen (His, if applicable Asp/Glu),
    # because Cys-OOP per Konstruktion ~0 ist.
    _oop_groups = [g for g in group_names if g not in cys_labels]
    _inp_groups = list(group_names)  # INP for alle Gruppen (Cys hat sinnvolle INP)
    
    feat_e: List[str] = (feat_b + _SCORE_KEYS
                         + [f"{g}_OOP" for g in _oop_groups]
                         + [f"{g}_INP" for g in _inp_groups]
                         + _per_lig)

    rows_b: List[List[float]] = []
    rows_e: List[List[float]] = []

    for r in results:
        ks   = r.get("kern_scores", {})
        grp  = r.get("groups", {})
        flig = r.get("fe_lig", {})

        # v3.5: Basis-Achsen are jetzt (lig_oop_pct, kern_oop, kern_d).
        # Die globale oop_pct was entfernt; lig_oop_pct (ligands-Sphaere)
        # is physikalisch aussagekraeftiger for the Mode-Klassifikation.
        b = [r.get("lig_oop_pct", 0.) / 100.,
             r["kern_oop"] / 100.,
             r["kern_d"]]
        b = [v if np.isfinite(v) else 0. for v in b]

        e = b[:]
        for sk in _SCORE_KEYS:
            v = ks.get(sk, 0.)
            e.append(float(v) if np.isfinite(float(v)) else 0.)
        # OOP only for Nicht-Cys-Gruppen
        for g in _oop_groups:
            gv = grp.get(g, {})
            e.append(gv.get("oop", 0.) / 100.)
        # INP for alle Gruppen
        for g in _inp_groups:
            gv = grp.get(g, {})
            e.append(gv.get("inp", 0.) / 100.)
        for l in lig_labels:
            lv = flig.get(l, {})
            e.append(lv.get("stretch", 0.))
            # v3.5: bend in INP- and OOP-Fractione getrennt (verbessert
            # Trennkraft of the Embeddings: OOP- and INP-Biegemoden landen
            # in differenten Cluster-Regionen).
            e.append(lv.get("bend_inp", lv.get("bend", 0.)))
            e.append(lv.get("bend_oop", 0.))

        rows_b.append(b)
        rows_e.append(e)

    X_b = np.nan_to_num(np.array(rows_b, dtype=float))
    X_e = np.nan_to_num(np.array(rows_e, dtype=float))
    return X_b, X_e, feat_b, feat_e


# ===========================================================================
# Embeddings + HDBSCAN
# ===========================================================================

def _auto_mcs(n: int) -> int:
    """Determines HDBSCAN ``min_cluster_size`` automatically.

    Heuristik (v3.1): ``max(5, min(30, round(sqrt(n))))``. Die alte
    3-Prozent-Regel war for DFT-Modes to permissiv; for N=1600 lieferte
    sie mcs=48, was for 18-dim Feature-Raeumen to 0 Cluster fuehrt
    (alles noise). Mit sqrt-Skalierung statt 3% bekommen wir mcs=40
    for N=1600, with Obergrenze 30 correspond tod 30; the findet bei
    typicalen DFT-Frequenz-Datensaetzen kleine Mode-Gruppen besser.
    """
    return max(5, min(30, int(round(n ** 0.5))))


def _auto_embed_params(n: int,
                       n_neighbors_override: Optional[int] = None,
                       perplexity_override:  Optional[float] = None,
                       ) -> Tuple[int, float]:
    """Computes UMAP- and t-SNE-Parameter after Empfehlung the Original-Paper.

    Beide Parameter are on ca. 1 % the Datenpunkte gesetzt,
    begrenzt on [5, 15].

    * **UMAP** ``n_neighbors``: McInnes et al. (2018) empfehlen
      values between 5 and 50; the Obergrenze 15 haelt alle
      frequency window (0-100, 0-500 cm⁻^1) on demselben value
      for direkte Vergleichbarkeit.  Formel: ``max(5, min(15, N // 100))``

    * **t-SNE** ``perplexity``: Van the Maaten & Hinton (2008) empfehlen
      values between 5 and 50. Gleiche Formel: ``max(5, min(15, N // 100))``

    For N = 1 539 (typicales [2Fe-2S]-System) gives sich
    ``n_neighbors = perplexity = 15`` - exakt the Paper-Standard.

    Parameters
    ----------
    n : int
        Number of Datenpunkte.
    n_neighbors_override : int, optional
        Manueller Override from ``cfg.umap_n_neighbors``.
    perplexity_override : float, optional
        Manueller Override from ``cfg.tsne_perplexity``.

    Returns
    -------
    n_neighbors : int
    perplexity : float

    References
    ----------
    McInnes, L., Healy, J. & Melville, J. (2018).
    UMAP: Uniform Manifold Approximation and Projection for Dimension
    Reduction. *arXiv*:1802.03426.

    Van the Maaten, L. & Hinton, G. (2008).
    Visualizing data using t-SNE.
    *Journal of Machine Learning Research*, 9, 2579-2605.
    """
    # Paper-empfohlene Formel: ~1% von N, Bereich [5, 15]
    # Obergrenze 15 = Paper-default (McInnes et al. 2018, Van the Maaten & Hinton 2008)
    # Hält alle frequency window (0-100, 100-300, 300-500, 0-500 cm⁻^1)
    # on demselben value → direkte Vergleichbarkeit the UMAP-Strukturen
    _auto = max(5, min(15, n // 100))
    n_nb  = int(n_neighbors_override) if n_neighbors_override else _auto
    perp  = float(perplexity_override) if perplexity_override  else float(_auto)
    # Sicherheitsgrenze: Parameter dürfen n not überschreiten
    n_nb  = min(n_nb, max(2, n - 1))
    perp  = min(perp, max(2., n - 1.))
    return n_nb, perp


def _hdbscan_on(Z2d: np.ndarray, min_size: int) -> np.ndarray:
    """Performs HDBSCAN-Clustering on a 2-D-Einbettung durch.

    Parameters
    ----------
    Z2d : ndarray of shape (n, 2)
        2-D-Koordinaten (e.g. UMAP- or t-SNE-Projektion).
    min_size : int
        Minimale Clustergroesse for HDBSCAN.

    Returns
    -------
    ndarray of int, shape (n,)
        Cluster-Labels; -1 bedeutet noise/kein Cluster.
        Returns a Array with -1 zurueck if hdbscan not installiert ist.
    """
    try:
        import hdbscan
        cl = hdbscan.HDBSCAN(min_cluster_size=min_size,
                              min_samples=1,
                              cluster_selection_method="eom",
                              prediction_data=True)
        return cl.fit_predict(Z2d)
    except ImportError:
        return np.full(len(Z2d), -1, dtype=int)


def compute_embeddings(X_b:     np.ndarray,
                        X_e:     np.ndarray,
                        feat_e:  List[str],
                        results: List[Dict],
                        runlog:  RunLog,
                        ) -> Tuple[Dict[str, np.ndarray], Dict]:
    """Computes UMAP on the erweiterten Feature-Matrix and fuehrt
    HDBSCAN-Clustering durch.

    v3.7.2: PCA, t-SNE and UMAP_Basis were entfernt. UMAP on der
    vollen Feature-Matrix bewahrt globale and lokale Struktur at the besten
    for DFT-Mode-Daten. justification in the Handbuch (Embeddings-Kapitel).

    Parameters
    ----------
    X_b : ndarray of shape (n_modes, 3)
        Basis-Feature-Matrix (NICHT mehr used, only fuer
        Backward-Kompatibilitaet in the Aufruf-Signatur).
    X_e : ndarray of shape (n_modes, n_features)
        Erweiterte Feature-Matrix.
    feat_e : list of str
        Feature-Namen the erweiterten Matrix.
    results : list of dict
        Modenanalyse-Ergebnisse (fuer Cluster-Charakterisierung).
    runlog : RunLog
        Fuer Warnmeldungen.

    Returns
    -------
    coords : dict of {str: ndarray of shape (n_modes, 2)}
        2D-Koordinaten ``{"UMAP": ...}``.
    cl_data : dict of {str: (labels, chars, cids)}
        HDBSCAN-Ergebnisse (wenn >= 2 Cluster found).
    """
    coords:  Dict[str, np.ndarray] = {}
    cl_data: Dict                  = {}

    n_modes = len(results)
    if n_modes < 3:
        runlog.warn(
            f"Zu wenige modes for Embedding ({n_modes} < 3) -- "
            f"Embedding analysis skipped.")
        return coords, cl_data

    # Standardisierung (pro Feature: Mean 0, Var 1) - necessary damit
    # UMAP-Distanzen not von Feature-Skalen dominiert werden.
    try:
        from sklearn.preprocessing import StandardScaler
        Xs = StandardScaler().fit_transform(X_e)
    except ImportError:
        Xs = X_e.copy()

    # UMAP on erweiterter Feature-Matrix
    try:
        import umap as umap_lib
        n_nb, _ = _auto_embed_params(
            len(results),
            n_neighbors_override=getattr(runlog.cfg, "umap_n_neighbors", None))
        coords["UMAP"] = umap_lib.UMAP(
            n_components=2, n_neighbors=n_nb,
            random_state=42).fit_transform(Xs)
        _umap_msg = (f"UMAP: n_neighbors={n_nb}"
                     f"  (auto, N={len(results)}, "
                     f"{X_e.shape[1]} Features)")
        runlog.info(_umap_msg)
        print(f"    {_umap_msg}")
    except ImportError:
        runlog.warn(
            "umap-learn not installed -- UMAP skipped. "
            "Installation: pip install umap-learn")

    # HDBSCAN on UMAP
    try:
        import hdbscan as _hdbscan_mod  # noqa: F401
        _hdbscan_avail = True
    except ImportError:
        _hdbscan_avail = False

    if not _hdbscan_avail:
        runlog.warn("hdbscan-Paket not installiert -- Cluster-Sheets "
                    "skipped (embedding coordinates remain "
                    "preserved in the sheets). Installation: "
                    "pip install hdbscan")
    else:
        mcs = _auto_mcs(len(results))
        for method, Z2d in coords.items():
            try:
                labels  = _hdbscan_on(Z2d, mcs)
                n_cl    = len(set(labels) - {-1})
                n_noise = int((labels == -1).sum())
                _hdb_msg = (f"HDBSCAN on {method}: "
                            f"{n_cl} clusters, {n_noise} noise (mcs={mcs})")
                runlog.info(_hdb_msg)
                print(f"    {_hdb_msg}")
                if n_cl >= 2:
                    chars, cids = characterize_clusters(
                        Z2d, labels, X_e, feat_e, results)
                    cl_data[method] = (labels, chars, cids)
                else:
                    runlog.info(
                        f"  ({method}: no sinnvollen Cluster found, "
                        f"sheet skipped)")
            except Exception as e:
                runlog.warn(f"HDBSCAN {method}: {e}")

    return coords, cl_data

    return coords, cl_data


# ===========================================================================
# Cluster-Charakterisierung
# ===========================================================================

def characterize_clusters(
        Z2d:        np.ndarray,
        labels:     np.ndarray,
        feat_matrix: np.ndarray,
        feat_names: List[str],
        results:    List[Dict],
        top_n:      int = 10,
) -> Tuple[Dict, List[int]]:
    """Charakterisiert HDBSCAN-Cluster statistisch.

    Parameters
    ----------
    Z2d : ndarray of shape (n_modes, 2)
        2D-Einbettungskoordinaten.
    labels : ndarray of shape (n_modes,)
        HDBSCAN-Cluster-Labels (``-1`` = noise).
    feat_matrix : ndarray of shape (n_modes, n_features)
        Feature-Matrix.
    feat_names : list of str
        Feature-Namen.
    results : list of dict
        Modenanalyse-Ergebnisse.
    top_n : int, optional
        Number of diskriminantester Features and repraesen-tativster modes.
        Standard: ``10``.

    Returns
    -------
    chars : dict
        Cluster-Charakterisierungen (Z-Scores, Fisher-F, Top-Moden).
    cluster_ids : list of int
        Sortierte Cluster-IDs (ohne ``-1``).

    Notes
    -----
    scipy-Fix: ``ConstantInputWarning`` von ``f_oneway`` is unterdrueckt.
    Sie tritt on if ein Feature over alle modes a Clusters
    konstant is (e.g. Null-Amplituden for bestimmte SS-elements).
    """
    try:
        from scipy.stats import f_oneway
    except ImportError:
        def f_oneway(*args):
            """Fallback if scipy not installiert: gibt (0, 1) ."""
            return (0., 1.)

    cluster_ids = sorted(k for k in set(labels) if k >= 0)
    n_feats     = feat_matrix.shape[1]
    g_mean      = feat_matrix.mean(0)
    g_std       = feat_matrix.std(0)
    g_std[g_std < 1e-12] = 1.0

    fisher = np.zeros(n_feats)
    for fi in range(n_feats):
        groups = [feat_matrix[labels == k, fi]
                  for k in cluster_ids if (labels == k).sum() >= 2]
        if len(groups) >= 2:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        message="Each of the input arrays is constant",
                        category=RuntimeWarning)
                    F, _ = f_oneway(*groups)
                fisher[fi] = float(F) if np.isfinite(F) else 0.
            except Exception:
                fisher[fi] = 0.  # numerisch not defined

    centers2d = {k: Z2d[labels == k].mean(0)
                 for k in cluster_ids if (labels == k).any()}
    chars: Dict = {}

    for k in list(cluster_ids) + ([-1] if -1 in labels else []):
        mask  = labels == k
        n_k   = int(mask.sum())
        sub   = feat_matrix[mask]
        means = sub.mean(0) if n_k > 0 else np.zeros(n_feats)
        z_sc  = (means - g_mean) / g_std

        disc     = np.abs(z_sc) * np.sqrt(np.maximum(fisher, 0))
        top_fi   = np.argsort(-disc)[:top_n]
        top_disc = [(feat_names[fi], float(z_sc[fi]), float(disc[fi]))
                    for fi in top_fi]

        # mask hat Länge len(labels) - for SS-UMAP < len(results)
        # Schneide range on len(mask) um IndexError to vermeiden
        freqs_k = np.array([results[i]["freq"]
                             for i in range(min(len(results), len(mask)))
                             if mask[i]])
        freq_stats = {
            "min":  float(freqs_k.min())  if len(freqs_k) else 0.,
            "max":  float(freqs_k.max())  if len(freqs_k) else 0.,
            "mean": float(freqs_k.mean()) if len(freqs_k) else 0.,
        }

        typ_dist: Dict[str, int] = {}
        for i in range(min(len(results), len(mask))):
            if mask[i]:
                t = results[i]["mode_type"]
                typ_dist[t] = typ_dist.get(t, 0) + 1

        idxs = np.where(mask)[0]
        if k >= 0 and k in centers2d and len(idxs) > 0:
            dists   = np.linalg.norm(Z2d[idxs] - centers2d[k], axis=1)
            top_idx = idxs[np.argsort(dists)[:top_n]]
        else:
            top_idx = idxs[:top_n]

        top_modes = [{"freq": results[i]["freq"],
                      "mode_type": results[i]["mode_type"],
                      "number":    results[i]["number"]}
                     for i in top_idx]
        all_modes = sorted([{"freq": results[i]["freq"],
                              "mode_type": results[i]["mode_type"],
                              "number":    results[i]["number"]}
                             for i in idxs], key=lambda x: x["freq"])

        chars[k] = {
            "n": n_k, "z_scores": z_sc, "means": means,
            "top_disc": top_disc, "top_n_modes": top_modes,
            "all_modes": all_modes, "freq_stats": freq_stats,
            "typ_dist": typ_dist, "is_noise": k == -1,
        }

    chars["_fisher"] = fisher
    return chars, cluster_ids


# ===========================================================================
# SS-UMAP-Clustering
# ===========================================================================

def compute_ss_umap_cluster(results:     List[Dict],
                              ss_elements: List[Dict],
                              runlog=None,
                              ) -> tuple:
    """UMAP clustering of secondary-structure amplitudes.

    Parameters
    ----------
    results : list of dict
        Modenanalyse-Ergebnisse; required ``results[i]["ss"]``.
    ss_elements : list of dict
        SS-element-Records from ``parse_pdb``.

    Returns
    -------
    Z2d : ndarray of shape (n_valid, 2) or None
        2D-UMAP-Koordinaten the modes with SS-Daten.
    full_labels : ndarray of shape (n_modes,) or None
        HDBSCAN-Labels for alle modes (``-99`` = no SS-Datum).
    feat_names : list of str
        Feature-Namen the SS-Amplituden-Matrix.
    X_norm : ndarray or None
        Normierte Feature-Matrix.
    valid_idx : list of int
        Indizes the modes with SS-Daten in ``results``.
    cluster_chars : dict
        Cluster-Charakterisierung (Z-Scores, Fisher-F, Top-Moden)
        from ``characterize_clusters``.
    """
    if not ss_elements:
        return None, None, [], None, []

    ss_names = [e["name"] for e in ss_elements]
    metrics  = ["amplitude_mean","amplitude_max","com_amplitude",
                "bending_std","bending_mean","stretching",
                "axial_amplitude","tilting_angle","internal_amplitude"]

    rows:       List[List[float]] = []
    valid_idx:  List[int]          = []
    for i, r in enumerate(results):
        ss = r.get("ss", {})
        if not ss: continue
        row = [float(ss.get(sn, {}).get(m, 0.))
               for sn in ss_names for m in metrics]
        rows.append(row)
        valid_idx.append(i)

    if len(rows) < 5:
        msg = (f"SS-UMAP: only {len(rows)} modes with SS data "
               f"(min: 5). SS-UMAP skipped.")
        if runlog is not None:
            runlog.warn(msg)
        else:
            print(f"    {msg}")
        return None, None, [], None, []

    feat_names = [f"{sn}_{m}" for sn in ss_names for m in metrics]
    X = np.array(rows, dtype=float)
    mu = X.mean(0, keepdims=True)
    sg = X.std(0, keepdims=True); sg[sg < 1e-12] = 1.
    X_norm = (X - mu) / sg

    try:
        from umap import UMAP
        nn, _ = _auto_embed_params(len(rows))
        Z2d  = UMAP(n_components=2, n_neighbors=nn,
                    random_state=42, low_memory=False).fit_transform(X_norm)
    except Exception as e_u:
        msg = f"SS-UMAP: UMAP fehlgeschlagen, Fallback on PCA ({type(e_u).__name__}: {e_u})"
        if runlog is not None:
            runlog.warn(msg)
        else:
            print(f"    {msg}")
        from sklearn.decomposition import PCA
        Z2d = PCA(n_components=2, random_state=42).fit_transform(X_norm)

    mcs    = _auto_mcs(len(rows))
    labels = _hdbscan_on(Z2d, mcs)
    n_cl   = len(set(labels) - {-1})
    n_ns   = int((labels == -1).sum())
    _ss_umap_msg = (f"SS-UMAP: {n_cl} clusters, {n_ns} noise, "
                    f"mcs={mcs}, n_modes={len(rows)}")
    if runlog is not None:
        runlog.info(_ss_umap_msg)
    else:
        print(f"    {_ss_umap_msg}")

    full_labels = np.full(len(results), -99, dtype=int)
    for li, gi in enumerate(valid_idx):
        full_labels[gi] = int(labels[li])

    # Cluster-Charakterisierung (Z-Scores + F-values + Top-Moden)
    # characterize_clusters gibt (chars_dict, cluster_ids)-Tuple zurueck
    cluster_chars: dict = {}
    _cids_ss = sorted(k for k in set(labels) if k >= 0)
    if _cids_ss and X_norm is not None:
        cluster_chars, _ = characterize_clusters(
            Z2d, labels, X_norm, feat_names, results)

    return Z2d, full_labels, feat_names, X_norm, valid_idx, cluster_chars

__version__ = "1.4"  # modenanalyse v1.4
