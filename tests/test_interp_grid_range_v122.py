# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regression tests for v1.2.2: the interpolation grid spans exactly the
requested frequency range.

Background
----------
Before v1.2.2 the interpolated Excel files (``_analysis_interp*.xlsx`` and
``_analysis_SSE_interp*.xlsx``) derived their frequency grid from
``cfg.freq_min`` / ``cfg.freq_max``, and fell back to
*first mode - interp_edge_extend* whenever those were ``None``. Two
situations made that surprising for users:

1. In multi-window mode the *overall* export deliberately clears
   ``freq_min``/``freq_max`` (so that ``Config.outdir()`` does not create a
   frequency subfolder). The top-level grid therefore started at the first
   real mode, e.g. 4.2179 cm^-1 for a 0-100 request.
2. With only an upper bound set (``freq_max = 800``, the shipped template
   default) the output folder was labelled ``0-800_cm-1`` while the grid
   still started at the first mode.

``export.interp_grid_bounds`` now implements a single rule: **a requested
range is reproduced exactly**; only a run without any requested range
follows the data.

Run::

    python3 -m pytest tests/test_interp_grid_range_v122.py -v
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

# Frequencies of a typical low-frequency run: no mode below 4.7 cm^-1,
# none above 99.9 -- exactly the situation that produced the truncated grid.
MODE_FREQS = [4.7179, 12.4, 33.9, 51.2, 78.6, 99.9]


def _cfg(**kw):
    from modenanalyse_2fe2s.config import Config
    cfg = Config(log_file="dummy.log")
    for k, v in kw.items():
        setattr(cfg, k, v)
    return cfg


def _grid(cfg, lo, hi):
    """Rebuild the grid exactly as the exporters do."""
    return np.arange(lo, hi + cfg.interp_step / 2, cfg.interp_step)


# =============================================================================
# 1) A requested range is reproduced exactly
# =============================================================================

def test_requested_window_is_honoured_exactly():
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_min=0.0, freq_max=100.0)
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS)
    assert lo == 0.0, f"grid must start at the requested 0.0, got {lo}"
    assert hi == 100.0, f"grid must end at the requested 100.0, got {hi}"


def test_grid_first_point_is_zero_not_first_mode():
    """The concrete user-visible symptom: first row of the Excel."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_min=0.0, freq_max=100.0, interp_step=0.05)
    grid = _grid(cfg, *interp_grid_bounds(cfg, MODE_FREQS))
    assert grid[0] == pytest.approx(0.0, abs=1e-12)
    assert grid[-1] == pytest.approx(100.0, abs=1e-9)
    # pre-v1.2.2 behaviour would have started at 4.7179 - 0.5 = 4.2179
    assert grid[0] < 4.2179


def test_points_below_first_mode_interpolate_to_zero():
    """Grid extension is only useful if the empty range reads as 0.0,
    not as an extrapolated value."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_min=0.0, freq_max=100.0, interp_step=0.05)
    grid = _grid(cfg, *interp_grid_bounds(cfg, MODE_FREQS))
    vals = np.interp(grid, MODE_FREQS, np.ones(len(MODE_FREQS)),
                     left=0.0, right=0.0)
    below = grid < MODE_FREQS[0]
    assert below.any(), "test setup: grid must reach below the first mode"
    assert np.all(vals[below] == 0.0)


def test_lower_bound_only_keeps_data_driven_upper_edge():
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_min=50.0, freq_max=None)
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS)
    assert lo == 50.0
    assert hi == pytest.approx(max(MODE_FREQS) + cfg.interp_edge_extend)


# =============================================================================
# 2) An upper bound alone implies the lower bound 0.0
# =============================================================================

def test_upper_bound_alone_implies_zero_lower_bound():
    """Matches Config.get_windows() (lo defaults to 0.0) and the
    ``0-800_cm-1`` subfolder label produced by Config.freq_label()."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_max=800.0)
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS)
    assert lo == 0.0
    assert hi == 800.0
    # the folder label the user sees must describe the same range
    assert cfg.freq_label() == "0-800_cm-1"


# =============================================================================
# 3) No requested range -> unchanged data-driven grid
# =============================================================================

def test_no_window_keeps_data_driven_grid():
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg()
    assert cfg.freq_min is None and cfg.freq_max is None
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS)
    assert lo == pytest.approx(min(MODE_FREQS) - cfg.interp_edge_extend)
    assert hi == pytest.approx(max(MODE_FREQS) + cfg.interp_edge_extend)


# =============================================================================
# 4) Explicit overrides (multi-window overall export)
# =============================================================================

def test_explicit_bounds_win_over_cleared_cfg():
    """runner.py clears freq_min/freq_max for the overall export so that
    no frequency subfolder is created; the window hull is passed
    explicitly instead."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg()                       # freq_min/freq_max both None
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS, grid_min=0.0, grid_max=100.0)
    assert (lo, hi) == (0.0, 100.0)


def test_explicit_bounds_override_cfg_values():
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_min=50.0, freq_max=100.0)
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS, grid_min=0.0, grid_max=300.0)
    assert (lo, hi) == (0.0, 300.0)


def test_window_hull_matches_config_windows():
    """The hull runner.py computes for the overall export must cover every
    window of the multi-window config."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_windows=[(0.0, 50.0), (50.0, 100.0)])
    windows = cfg.get_windows()
    hull_lo = min(lo for lo, _ in windows)
    hull_hi = max(hi for _, hi in windows)
    assert (hull_lo, hull_hi) == (0.0, 100.0)
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS,
                                grid_min=hull_lo, grid_max=hull_hi)
    assert (lo, hi) == (0.0, 100.0)


# =============================================================================
# 5) Robustness
# =============================================================================

def test_infinite_upper_bound_falls_back_to_data():
    """An open window (hi = inf) must not produce an infinite grid."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg(freq_min=50.0)
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS, grid_max=float("inf"))
    assert lo == 50.0
    assert np.isfinite(hi)
    assert hi == pytest.approx(max(MODE_FREQS) + cfg.interp_edge_extend)


def test_degenerate_range_falls_back_to_data():
    """hi <= lo would make np.arange return an empty grid (Excel with no
    frequency columns); fall back to the data range instead."""
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg()
    lo, hi = interp_grid_bounds(cfg, MODE_FREQS, grid_min=100.0, grid_max=10.0)
    assert lo < hi
    assert lo == pytest.approx(min(MODE_FREQS) - cfg.interp_edge_extend)
    assert hi == pytest.approx(max(MODE_FREQS) + cfg.interp_edge_extend)
    assert len(_grid(cfg, lo, hi)) > 0


def test_single_mode_run_produces_nonempty_grid():
    from modenanalyse_2fe2s.export import interp_grid_bounds
    cfg = _cfg()
    lo, hi = interp_grid_bounds(cfg, [42.0])
    assert len(_grid(cfg, lo, hi)) > 0


# =============================================================================
# 6) Both interpolation exporters accept the explicit bounds
# =============================================================================

def test_both_interp_exporters_accept_grid_bounds():
    """Guards against the SSE exporter drifting out of sync with the core
    exporter (they used to duplicate the bounds logic)."""
    import inspect
    from modenanalyse_2fe2s import export

    for fn in (export.export_interpolated_excel, export.export_sse_interp_excel):
        params = inspect.signature(fn).parameters
        assert "grid_min" in params, f"{fn.__name__} lacks grid_min"
        assert "grid_max" in params, f"{fn.__name__} lacks grid_max"
        assert params["grid_min"].default is None
        assert params["grid_max"].default is None


def test_export_payload_carries_grid_bounds():
    from modenanalyse_2fe2s.export import ExportPayload
    fields = ExportPayload.__dataclass_fields__
    assert "grid_min" in fields and "grid_max" in fields
    payload = ExportPayload(results=[], coord_info=None, dist_ref={},
                            logname="x", cfg=None, runlog=None)
    assert payload.grid_min is None and payload.grid_max is None
