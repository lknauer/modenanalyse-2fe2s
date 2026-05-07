# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit-Tests fuer modenanalyse_2fe2s.config.

Testet:
* Validierungs-Logik in Config.validate() fuer alle Edge-Cases
  (freq-Reihenfolge, Window-Format, negative Werte, ungueltige Modi)
* TOML-Roundtrip via Config.from_toml + Config.to_toml
* Backward-Compat: v3.6/v3.7.1-Schluessel werden mit UserWarning
  akzeptiert statt zu crashen (_legacy_silent_drop)
* Section-Flatten ([input], [freq], [pcet] etc.)

Aufruf::

    python -m pytest tests/test_config.py -v
    python tests/test_config.py
"""
import os
import sys
import tempfile
import warnings

# Pfad-Setup fuer direkten Aufruf via python tests/test_config.py
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "src"))

import dataclasses
import pytest

from modenanalyse_2fe2s.config import Config


# =============================================================================
# Default-Konstruktion und Validierung
# =============================================================================

def test_default_config_constructs():
    """Config() ohne Argumente baut eine Instanz mit Defaults."""
    cfg = Config()
    assert cfg.log_file == ""
    assert cfg.output_dir == ""
    assert cfg.cluster_index == 0
    assert cfg.pcet_enabled is True
    assert cfg.analyze_scsd is True


def test_default_config_validate_complains_about_logfile():
    """Default-Config: validate() meldet fehlendes log_file."""
    cfg = Config()
    errors = cfg.validate()
    assert any("log_file" in e for e in errors), (
        f"Erwartet log_file-Fehler, bekommen: {errors}")


def test_validate_freq_min_ge_max():
    """freq_min >= freq_max wird gemeldet."""
    cfg = Config(freq_min=500.0, freq_max=300.0)
    errors = cfg.validate()
    assert any("freq_min" in e and "freq_max" in e for e in errors), (
        f"Erwartet freq_min/max-Fehler, bekommen: {errors}")


def test_validate_freq_min_lt_max_ok():
    """freq_min < freq_max wird akzeptiert (kein freq-Fehler)."""
    cfg = Config(freq_min=100.0, freq_max=500.0)
    errors = cfg.validate()
    assert not any("freq_min" in e and "freq_max" in e for e in errors), (
        f"Unerwartet freq-Fehler: {errors}")


def test_validate_negative_temp_k():
    """Negative oder 0 K Temperatur wird abgefangen."""
    for bad_t in (-10.0, 0.0):
        cfg = Config(temp_k=bad_t)
        errors = cfg.validate()
        assert any("temp_k" in e for e in errors), (
            f"temp_k={bad_t}: Fehler erwartet, bekommen: {errors}")


def test_validate_freq_windows_invalid_format():
    """freq_windows mit kaputtem Format wird abgefangen."""
    cfg = Config(freq_windows=[(100, 300), "kaputt"])
    errors = cfg.validate()
    assert any("freq_windows" in e and "[1]" in e for e in errors), (
        f"Erwartet freq_windows[1]-Fehler, bekommen: {errors}")


def test_validate_freq_windows_lo_ge_hi():
    """freq_windows mit lo >= hi wird abgefangen."""
    cfg = Config(freq_windows=[(100, 300), (500, 400)])
    errors = cfg.validate()
    assert any("freq_windows" in e and "[1]" in e for e in errors), (
        f"Erwartet freq_windows[1]-Reihenfolge-Fehler, bekommen: {errors}")


def test_validate_freq_windows_valid_list():
    """Eine gueltige freq_windows-Liste erzeugt keinen freq_windows-Fehler."""
    cfg = Config(freq_windows=[(0, 100), (100, 300), (300, 500)])
    errors = cfg.validate()
    assert not any("freq_windows" in e for e in errors), (
        f"Unerwartet freq_windows-Fehler: {errors}")


def test_validate_negative_cluster_index():
    """cluster_index < 0 wird abgefangen."""
    cfg = Config(cluster_index=-1)
    errors = cfg.validate()
    assert any("cluster_index" in e for e in errors), (
        f"Erwartet cluster_index-Fehler, bekommen: {errors}")


def test_validate_invalid_interp_boundary_mode():
    """Nicht-erlaubter interp_boundary_mode wird abgefangen."""
    cfg = Config(interp_boundary_mode="foobar")
    errors = cfg.validate()
    assert any("interp_boundary_mode" in e for e in errors), (
        f"Erwartet interp_boundary_mode-Fehler, bekommen: {errors}")


def test_validate_valid_interp_modes():
    """Alle drei erlaubten Boundary-Modi gehen durch."""
    for mode in ("context", "zero", "nearest"):
        cfg = Config(interp_boundary_mode=mode)
        errors = cfg.validate()
        assert not any("interp_boundary_mode" in e for e in errors), (
            f"Mode '{mode}': unerwartet abgelehnt: {errors}")


def test_validate_negative_pcet_hbond_cutoff():
    """pcet_hbond_cutoff_a <= 0 wird abgefangen."""
    cfg = Config(pcet_hbond_cutoff_a=-1.0)
    errors = cfg.validate()
    assert any("pcet_hbond_cutoff_a" in e for e in errors), (
        f"Erwartet pcet_hbond_cutoff_a-Fehler, bekommen: {errors}")


def test_validate_negative_reorg_spectrum_params():
    """reorg_spectrum_sigma_cm1/step_cm1 muessen positiv sein."""
    cfg = Config(reorg_spectrum_sigma_cm1=-5.0)
    errors = cfg.validate()
    assert any("reorg_spectrum_sigma_cm1" in e for e in errors), errors

    cfg = Config(reorg_spectrum_step_cm1=0.0)
    errors = cfg.validate()
    assert any("reorg_spectrum_step_cm1" in e for e in errors), errors


def test_validate_negative_fe_coord_cutoffs():
    """fe_coord_cutoff_n/s/o muessen positiv sein."""
    cfg = Config(fe_coord_cutoff_s=-1.0)
    errors = cfg.validate()
    assert any("fe_coord_cutoff" in e for e in errors), (
        f"Erwartet fe_coord_cutoff-Fehler, bekommen: {errors}")


def test_validate_negative_sigma_eigvec_coord():
    """sigma_eigvec, sigma_coord duerfen nicht negativ sein."""
    cfg = Config(sigma_eigvec=-1.0)
    errors = cfg.validate()
    assert any("sigma_eigvec" in e for e in errors), errors

    cfg = Config(sigma_coord=-0.5)
    errors = cfg.validate()
    assert any("sigma_coord" in e for e in errors), errors


def test_validate_zero_sigma_eigvec_ok():
    """sigma_eigvec = 0 ist erlaubt (deaktiviert Fehlerfortpflanzung)."""
    cfg = Config(sigma_eigvec=0.0)
    errors = cfg.validate()
    assert not any("sigma_eigvec" in e for e in errors), errors


# =============================================================================
# TOML-Read: Section-Flatten
# =============================================================================

def _write_temp_toml(content: str) -> str:
    """Schreibt content in ein temp .toml-File und gibt den Pfad zurueck."""
    fd, path = tempfile.mkstemp(suffix=".toml")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(content)
        return path
    except Exception:
        os.close(fd)
        os.unlink(path)
        raise


def test_from_toml_minimal():
    """Minimales TOML mit nur log_file und output_dir wird gelesen."""
    toml = (
        '[input]\n'
        'log_file = "test.log"\n'
        'output_dir = "results"\n'
    )
    path = _write_temp_toml(toml)
    try:
        cfg = Config.from_toml(path)
        assert cfg.log_file == "test.log"
        assert cfg.output_dir == "results"
    finally:
        os.unlink(path)


def test_from_toml_section_flatten():
    """Sub-Sections [input], [freq], [pcet] werden flach in Config gemappt."""
    toml = (
        '[input]\n'
        'log_file = "x.log"\n'
        'output_dir = "out"\n'
        '\n'
        '[freq]\n'
        'freq_max = 800.0\n'
        '\n'
        '[pcet]\n'
        'pcet_hbond_cutoff_a = 4.5\n'
    )
    path = _write_temp_toml(toml)
    try:
        cfg = Config.from_toml(path)
        assert cfg.log_file == "x.log"
        assert cfg.output_dir == "out"
        assert cfg.freq_max == 800.0
        assert cfg.pcet_hbond_cutoff_a == 4.5
    finally:
        os.unlink(path)


def test_from_toml_freq_windows_list_to_tuple():
    """freq_windows als TOML [[100,300],[300,500]] wird zu List[Tuple]."""
    toml = (
        '[input]\n'
        'log_file = "x.log"\n'
        'output_dir = "out"\n'
        '[freq]\n'
        'freq_windows = [[0, 100], [100, 300], [300, 500]]\n'
    )
    path = _write_temp_toml(toml)
    try:
        cfg = Config.from_toml(path)
        assert cfg.freq_windows == [(0.0, 100.0), (100.0, 300.0), (300.0, 500.0)]
        assert all(isinstance(w, tuple) for w in cfg.freq_windows)
    finally:
        os.unlink(path)


def test_from_toml_unknown_key_raises():
    """Unbekannte Top-Level-Schluessel werfen ValueError."""
    toml = (
        '[input]\n'
        'log_file = "x.log"\n'
        'output_dir = "out"\n'
        'voellig_unbekannt = 42\n'
    )
    path = _write_temp_toml(toml)
    try:
        with pytest.raises(ValueError, match="voellig_unbekannt|unknown|unbekannt"):
            Config.from_toml(path)
    finally:
        os.unlink(path)


def test_from_toml_duplicate_key_raises():
    """Schluessel sowohl Top-Level als auch in Sektion -> ValueError."""
    toml = (
        'log_file = "x.log"\n'
        '[input]\n'
        'log_file = "y.log"\n'
    )
    path = _write_temp_toml(toml)
    try:
        with pytest.raises(ValueError, match="doppelt|duplicate"):
            Config.from_toml(path)
    finally:
        os.unlink(path)


# =============================================================================
# Backward-Compat: v3.6/v3.7.1-Legacy-Schluessel
# =============================================================================

def test_from_toml_legacy_v36_pcet_warns_not_crashes():
    """v3.6 PCET-Score-Parameter werden mit UserWarning ignoriert."""
    toml = (
        '[input]\n'
        'log_file = "x.log"\n'
        'output_dir = "out"\n'
        '[pcet]\n'
        'pcet_strong_threshold = 0.5\n'
        'pcet_moderate_threshold = 0.2\n'
    )
    path = _write_temp_toml(toml)
    try:
        with warnings.catch_warnings(record=True) as warned:
            warnings.simplefilter("always")
            cfg = Config.from_toml(path)
            user_warnings = [w for w in warned if issubclass(w.category, UserWarning)]
            assert len(user_warnings) >= 2, (
                f"Erwartet mindestens 2 UserWarnings (eine pro legacy key), "
                f"bekommen: {len(user_warnings)}")
            # Die Felder duerfen NICHT auf der Config gelandet sein
            assert not hasattr(cfg, "pcet_strong_threshold")
        assert cfg.log_file == "x.log"  # Rest ist normal verarbeitet
    finally:
        os.unlink(path)


def test_from_toml_legacy_v371_tsne_warns():
    """v3.7.1 t-SNE-Parameter wird mit UserWarning ignoriert."""
    toml = (
        '[input]\n'
        'log_file = "x.log"\n'
        'output_dir = "out"\n'
        '[embedding]\n'
        'tsne_perplexity = 30.0\n'
    )
    path = _write_temp_toml(toml)
    try:
        with warnings.catch_warnings(record=True) as warned:
            warnings.simplefilter("always")
            cfg = Config.from_toml(path)
            user_warnings = [w for w in warned if issubclass(w.category, UserWarning)]
            assert any("tsne" in str(w.message).lower() for w in user_warnings), (
                f"Erwartet UserWarning ueber tsne_perplexity, bekommen: "
                f"{[str(w.message) for w in user_warnings]}")
    finally:
        os.unlink(path)


def test_from_toml_legacy_lambda0_params_warn():
    """v3.6 cpet/pt/et_lambda0_cm1-Parameter werden gewarnt + ignoriert."""
    toml = (
        '[input]\n'
        'log_file = "x.log"\n'
        'output_dir = "out"\n'
        '[pcet]\n'
        'cpet_lambda0_cm1 = 100.0\n'
        'pt_lambda0_cm1 = 80.0\n'
        'et_lambda0_cm1 = 120.0\n'
    )
    path = _write_temp_toml(toml)
    try:
        with warnings.catch_warnings(record=True) as warned:
            warnings.simplefilter("always")
            cfg = Config.from_toml(path)
            assert len(warned) >= 3, (
                f"Erwartet mindestens 3 Warnings, bekommen: {len(warned)}")
    finally:
        os.unlink(path)


# =============================================================================
# TOML-Roundtrip: from_toml(to_toml(cfg)) ~ cfg
# =============================================================================

def test_to_toml_then_from_toml_roundtrip():
    """Roundtrip: schreiben + lesen produziert die gleichen Felder."""
    cfg_orig = Config(
        log_file="some.log",
        output_dir="some_results",
        pdb_file="some.pdb",
        temp_k=5.0,
        freq_max=800.0,
        freq_windows=[(0, 100), (100, 300), (300, 500)],
        cluster_index=0,
        pcet_hbond_cutoff_a=4.5,
        analyze_scsd=True,
    )
    fd, path = tempfile.mkstemp(suffix=".toml")
    os.close(fd)
    try:
        cfg_orig.to_toml(path)
        cfg_loaded = Config.from_toml(path)
        # Vergleich der relevanten Felder
        for fname in ("log_file", "output_dir", "pdb_file",
                      "temp_k", "freq_max", "cluster_index",
                      "pcet_hbond_cutoff_a", "analyze_scsd"):
            assert getattr(cfg_orig, fname) == getattr(cfg_loaded, fname), (
                f"Feld {fname}: Original {getattr(cfg_orig, fname)!r}, "
                f"geladen {getattr(cfg_loaded, fname)!r}")
        # freq_windows: List[Tuple] vs List[Tuple]
        assert cfg_orig.freq_windows == cfg_loaded.freq_windows
    finally:
        os.unlink(path)


def test_config_has_expected_n_fields():
    """Anzahl Config-Felder bleibt stabil (Regression-Test)."""
    fields = dataclasses.fields(Config)
    # 48 Felder in v1.0.0:
    #   - 47 nach NIS-Cleanup in v3.7.5 (vorher 58 mit 11 NIS-Feldern;
    #     NIS-Berechnungen erfolgen jetzt ausschliesslich in nisspec3).
    #   - +1 in v1.0.0: analyze_all_clusters (Multi-Cluster-Modus).
    assert len(fields) == 48, (
        f"Erwartet 48 Felder, bekommen {len(fields)}. Falls bewusst ein "
        f"Feld hinzugefuegt/entfernt wurde, bitte hier aktualisieren.")


# =============================================================================
# Standalone-Aufruf
# =============================================================================

if __name__ == "__main__":
    # Manueller Aufruf ohne pytest -- iteriert ueber alle Tests im Modul
    import inspect

    test_funcs = [
        (name, fn) for name, fn in inspect.getmembers(sys.modules[__name__])
        if name.startswith("test_") and callable(fn)
    ]

    print(f"Running {len(test_funcs)} tests in {os.path.basename(__file__)}\n")

    n_pass = 0
    n_fail = 0
    for name, fn in test_funcs:
        try:
            fn()
            print(f"  [OK]   {name}")
            n_pass += 1
        except Exception as exc:
            print(f"  [FAIL] {name}: {exc!r}")
            n_fail += 1

    print(f"\n{n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)
