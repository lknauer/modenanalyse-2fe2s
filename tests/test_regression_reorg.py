# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Headline-Regressionstest: pinnt die aggregierten Reorganisations-
Kennzahlen des realen Cys4-Fixture-Logfiles auf enge Toleranzen.

Motivation
----------
Die uebrige Testsuite prueft *Verhaeltnisse* (FeS > FeFe), *Grenzfaelle*
(ZPE vs. klassisch) und *Invarianten* (Additivitaet, OOP+INP=1), aber
KEINE absoluten Ende-zu-Ende-Zahlen im Default-Lauf. Ein globaler
Skalenfehler (z.B. ein Faktor-2 im 1/2 mu omega^2-Vorfaktor oder eine
Angstrom<->m-Verwechslung) wuerde alle lambda gleichmaessig skalieren und
saemtliche Verhaeltnis-/Grenzwerttests unbemerkt passieren lassen.

Dieser Test faengt genau das: er faehrt die volle Pipeline einmal auf dem
mitgelieferten Cys_2Fe-2S-Fixture und vergleicht Lambda_FeS / Lambda_FeFe
(pair und mode) gegen fixierte Referenzwerte mit +/-2 % Toleranz.

Referenzwerte
-------------
Ermittelt am Stand des Branches ``critical-review-fixes`` (nach dem
Massen-/Normalisierungs-Audit) auf dem Fixture
``Cys_2Fe-2S_red3_hpfrq_opt.log`` bei ``temp_k=5.0``, ``freq_max=800.0``.
Aendert sich hier etwas ausserhalb der Toleranz, ist das ein Signal fuer
eine echte (beabsichtigte oder versehentliche) Physik-Aenderung -- dann
die Werte bewusst neu pinnen und im CHANGELOG dokumentieren.

Anders als ``test_smoke.py`` ist dieser Test NICHT als ``slow`` markiert:
er soll im Default-Lauf mitlaufen, weil er die wertvollste Absicherung
gegen eine stille Physik-Regression ist. Die Pipeline laeuft dafuer
einmal (~30 s, modul-weite Fixture).
"""

from __future__ import annotations
import lzma
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

# src-Layout: Paketpfad fuer Standalone-Lauf verfuegbar machen.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "src"))


DATA_DIR = Path(__file__).parent / "data"
LOG_XZ = DATA_DIR / "Cys_2Fe-2S_red3_hpfrq_opt.log.xz"

# Fixierte Referenzwerte (cm^-1) und relative Toleranz.
REF = {
    "FeFe": {"pair": 54.1181, "mode": 21.6186, "n": 66},
    "FeS":  {"pair": 286.6039, "mode": 206.3856, "n": 66},
}
REL_TOL = 0.02  # +/-2 %


@pytest.fixture(scope="module")
def reorg_totals(tmp_path_factory):
    """Faehrt die Pipeline einmalig und liefert das Reorg_total-Dict."""
    import openpyxl
    from modenanalyse_2fe2s import Config, run_analysis

    assert LOG_XZ.exists(), f"Test-Datei fehlt: {LOG_XZ}"

    workdir = tmp_path_factory.mktemp("reorg_reg")
    log_path = workdir / "Cys_2Fe-2S_red3_hpfrq_opt.log"
    with lzma.open(LOG_XZ, "rb") as src, open(log_path, "wb") as dst:
        shutil.copyfileobj(src, dst)

    output_dir = workdir / "results"
    cfg = Config(
        log_file=str(log_path),
        output_dir=str(output_dir),
        pdb_file="",
        temp_k=5.0,
        freq_max=800.0,
        analyze_scsd=False,   # nicht noetig fuer Reorg-Zahlen; spart Zeit
        pcet_enabled=True,
        use_cache=False,
    )
    run_analysis(cfg)

    excel = list((output_dir / "0-800_cm-1").glob("*_analysis.xlsx"))
    assert excel, f"_analysis.xlsx nicht gefunden in {output_dir}"
    wb = openpyxl.load_workbook(excel[0], read_only=True, data_only=True)
    ws = wb["Reorg_total"]

    vals = {}
    for row in list(ws.values)[1:]:
        if not row or row[0] is None:
            continue
        ch = str(row[0]).strip()
        if ch in ("FeFe", "FeN", "FeS", "NH", "HA"):
            vals[ch] = {
                "pair": float(row[1]) if row[1] is not None else 0.0,
                "mode": float(row[2]) if row[2] is not None else 0.0,
                "n":    int(row[3]) if row[3] is not None else 0,
            }
    return vals


@pytest.mark.parametrize("ch", ["FeFe", "FeS"])
@pytest.mark.parametrize("quantity", ["pair", "mode"])
def test_reorg_headline_value_pinned(reorg_totals, ch, quantity):
    """Lambda_{ch}_{quantity} liegt innerhalb +/-2 % des Referenzwerts."""
    got = reorg_totals[ch][quantity]
    ref = REF[ch][quantity]
    assert np.isfinite(got), f"Lambda_{ch}_{quantity} nicht finite: {got}"
    rel = abs(got - ref) / ref
    assert rel <= REL_TOL, (
        f"Lambda_{ch}_{quantity} = {got:.4f} cm-1 weicht {rel*100:.2f} % vom "
        f"gepinnten Referenzwert {ref:.4f} cm-1 ab (Toleranz {REL_TOL*100:.0f} %). "
        f"Falls beabsichtigt: Referenz in REF neu setzen + CHANGELOG."
    )


def test_reorg_mode_counts_pinned(reorg_totals):
    """Anzahl beitragender Modes ist stabil (Fixture-spezifisch)."""
    for ch in ("FeFe", "FeS"):
        assert reorg_totals[ch]["n"] == REF[ch]["n"], (
            f"{ch}: n_modes = {reorg_totals[ch]['n']}, erwartet {REF[ch]['n']}"
        )


def test_his_channels_zero(reorg_totals):
    """Cys4-Cluster: His-abhaengige Kanaele (FeN, NH, HA) sind null."""
    for ch in ("FeN", "NH", "HA"):
        assert reorg_totals[ch]["mode"] == 0.0, (
            f"{ch} muss 0 sein (kein His), ist {reorg_totals[ch]['mode']}"
        )
