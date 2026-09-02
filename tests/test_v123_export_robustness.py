# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regression tests for v1.2.3: export robustness in multi-window mode.

Background
----------
A user run (multi-window, ``freq_windows = [[0, 100]]``) produced an
interpolated Excel that looked "un-interpolated". The REPORT showed why:

1. ``interpolation Excel failed: [Errno 13] Permission denied`` -- the
   file of a previous run was open in Excel, the new one was silently
   dropped, and the stale file was mistaken for the result.
2. ``97 context modes right`` followed by ``interp_boundary_mode='context'
   but no context modes available`` -- the per-window export drew its
   context candidates from ``results``, which in multi-window mode is
   already filtered to the window hull.
3. ``output_dir: ./out\\0-100_cm-1`` in the REPORT although REPORT and
   overall export live in ``./out``.

Run::

    python3 -m pytest tests/test_v123_export_robustness.py -v
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

openpyxl = pytest.importorskip("openpyxl")


def _cfg(**kw):
    from modenanalyse_2fe2s.config import Config
    cfg = Config(log_file="dummy.log")
    for k, v in kw.items():
        setattr(cfg, k, v)
    return cfg


def _runlog():
    from modenanalyse_2fe2s.logio import RunLog
    return RunLog(_cfg())


# =============================================================================
# 1) Locked target file -> fallback name, loud warning
# =============================================================================

class _LockedWorkbook:
    """Workbook stand-in whose save() raises PermissionError for ``locked``."""

    def __init__(self, locked):
        self.locked = set(locked)
        self.saved = []

    def save(self, path):
        if path in self.locked:
            raise PermissionError(13, "Permission denied", path)
        self.saved.append(path)


def test_save_workbook_plain(tmp_path):
    from modenanalyse_2fe2s.export import _save_workbook
    out = str(tmp_path / "x_analysis_interp0.05.xlsx")
    wb = _LockedWorkbook([])
    rl = _runlog()
    assert _save_workbook(wb, out, rl) == out
    assert wb.saved == [out]
    assert rl.warnings == []


def test_save_workbook_locked_falls_back_to_numbered_name(tmp_path, capsys):
    from modenanalyse_2fe2s.export import _save_workbook
    out = str(tmp_path / "x_analysis_interp0.05.xlsx")
    wb = _LockedWorkbook([out])
    rl = _runlog()
    written = _save_workbook(wb, out, rl)
    assert written == str(tmp_path / "x_analysis_interp0.05_1.xlsx")
    assert wb.saved == [written]
    # reported in the REPORT ...
    assert len(rl.warnings) == 1
    assert "locked" in rl.warnings[0] and "_1.xlsx" in rl.warnings[0]
    # ... and on the console
    assert "locked" in capsys.readouterr().out


def test_save_workbook_skips_locked_fallbacks(tmp_path):
    from modenanalyse_2fe2s.export import _save_workbook
    out = str(tmp_path / "x.xlsx")
    wb = _LockedWorkbook([out, str(tmp_path / "x_1.xlsx")])
    written = _save_workbook(wb, out, _runlog())
    assert written == str(tmp_path / "x_2.xlsx")


def test_save_workbook_real_openpyxl_roundtrip(tmp_path):
    """The real Workbook goes through the same helper (no monkeypatching)."""
    from modenanalyse_2fe2s.export import _save_workbook
    wb = openpyxl.Workbook()
    out = str(tmp_path / "real.xlsx")
    assert _save_workbook(wb, out, _runlog()) == out
    assert os.path.exists(out)


# =============================================================================
# 2) Context modes for a window come from the full pool
# =============================================================================

def _modes(*freqs):
    return [{"freq": float(f)} for f in freqs]


def test_window_context_uses_modes_outside_hull():
    """The situation of the user run: window 0-100, analysed modes end at
    99.94, context modes 100-105 analysed separately."""
    from modenanalyse_2fe2s.runner import _window_context
    results = _modes(4.75, 50.0, 99.94)          # already window-filtered
    ctx_above = _modes(100.3, 102.0, 104.9, 105.0, 105.1)
    right, left = _window_context(results + ctx_above, 0.0, 100.0, 5.0)
    assert [r["freq"] for r in right] == [100.3, 102.0, 104.9, 105.0]
    assert left == []


def test_window_context_neighbour_window_still_supplies_context():
    from modenanalyse_2fe2s.runner import _window_context
    pool = _modes(45.0, 49.0, 50.0, 52.0, 55.0, 56.0)
    right, left = _window_context(pool, 50.0, 100.0, 5.0)
    assert [r["freq"] for r in left] == [45.0, 49.0]
    assert right == []


def test_window_context_open_window_has_no_right_context():
    from modenanalyse_2fe2s.runner import _window_context
    pool = _modes(100.5, 700.0, 703.0)
    right, left = _window_context(pool, 500.0, None, 5.0)
    assert right == []
    assert left == []


def test_window_context_sorted_even_if_pool_is_not():
    from modenanalyse_2fe2s.runner import _window_context
    pool = _modes(104.0, 101.0, 99.0, 96.0, 97.5)
    right, left = _window_context(pool, 100.0, 100.0 + 0.0, 5.0)
    assert [r["freq"] for r in right] == [101.0, 104.0]
    assert [r["freq"] for r in left] == [96.0, 97.5, 99.0]


# =============================================================================
# 3) REPORT shows the base directory in multi-window mode
# =============================================================================

def test_report_output_dir_multi_window(tmp_path):
    from modenanalyse_2fe2s.logio import RunLog
    base = str(tmp_path / "out")
    cfg = _cfg(output_dir=base, freq_min=0.0, freq_max=100.0,
               freq_windows=[(0.0, 100.0)])
    rl = RunLog(cfg)
    report = tmp_path / "REPORT.txt"
    rl.write_befund(str(report))
    txt = report.read_text(encoding="utf-8")
    line = next(l for l in txt.splitlines() if "output_dir:" in l)
    assert base in line
    assert "0-100_cm-1" not in line
    assert "freq_windows: [0-100] cm-1" in txt
    assert "ignored" in txt          # freq_min/freq_max hint


def test_report_output_dir_single_window_unchanged(tmp_path):
    from modenanalyse_2fe2s.logio import RunLog
    base = str(tmp_path / "out")
    cfg = _cfg(output_dir=base, freq_min=0.0, freq_max=100.0)
    rl = RunLog(cfg)
    report = tmp_path / "REPORT.txt"
    rl.write_befund(str(report))
    txt = report.read_text(encoding="utf-8")
    line = next(l for l in txt.splitlines() if "output_dir:" in l)
    assert "0-100_cm-1" in line
    assert "freq_filter: 0.0 - 100.0 cm-1" in txt
