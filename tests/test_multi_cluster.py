"""Tests fuer Multi-Cluster-Auswertung (v1.0.0+).

Testen die neue ``analyze_all_clusters``-API ohne den vollen Lauf
(verhindert lange Smoke-Test-Zeiten in der Default-Suite).
"""
from __future__ import annotations
import dataclasses
import lzma
from pathlib import Path

import pytest

from modenanalyse_2fe2s import Config

HERE = Path(__file__).parent
TEST_DATA = HERE / "data" / "Cys_2Fe-2S_red3_hpfrq_opt.log.xz"


def test_analyze_all_clusters_field_exists():
    """Config hat das Feld analyze_all_clusters mit Default False."""
    cfg = Config(log_file="dummy.log", output_dir="/tmp/dummy")
    assert hasattr(cfg, "analyze_all_clusters")
    assert cfg.analyze_all_clusters is False


def test_analyze_all_clusters_field_settable():
    """Feld kann auf True gesetzt werden."""
    cfg = Config(log_file="dummy.log", output_dir="/tmp/dummy",
                 analyze_all_clusters=True)
    assert cfg.analyze_all_clusters is True


def test_dataclass_replace_works_with_field():
    """dataclasses.replace funktioniert mit dem neuen Feld -- wichtig
    fuer den Multi-Cluster-Wrapper, der pro Cluster eine Config-Kopie
    erzeugt."""
    cfg = Config(log_file="dummy.log", output_dir="/tmp/dummy",
                 analyze_all_clusters=True)
    cfg2 = dataclasses.replace(cfg, cluster_index=1,
                                output_dir="/tmp/dummy/cluster_1",
                                analyze_all_clusters=False)
    assert cfg2.cluster_index == 1
    assert cfg2.output_dir == "/tmp/dummy/cluster_1"
    assert cfg2.analyze_all_clusters is False


@pytest.mark.skipif(not TEST_DATA.exists(),
                    reason="Test-Logfile nicht vorhanden")
def test_single_cluster_in_cys4_test_data(tmp_path):
    """Im Cys4-Test-Logfile gibt es genau 1 Cluster.

    find_all_clusters auf der entpackten Atomliste muss n=1 ergeben.
    """
    from modenanalyse_2fe2s.geometry import find_all_clusters
    from modenanalyse_2fe2s.logio import scan_log, read_std_orient

    # Test-File ist xz-komprimiert; nach tmp_path entpacken
    log_path = tmp_path / "test.log"
    with lzma.open(TEST_DATA, "rt") as f, open(log_path, "w") as out:
        out.write(f.read())

    cfg = Config(log_file=str(log_path), output_dir=str(tmp_path / "out"))
    so_off, _, _ = scan_log(str(log_path), cfg)
    atoms, _ = read_std_orient(str(log_path), so_off[-1],
                                include_hydrogen=False)
    clusters = find_all_clusters(atoms, cfg)
    assert len(clusters) == 1, (
        f"Erwartet 1 Cluster im Cys4-Modell, gefunden {len(clusters)}")


def test_multi_cluster_with_n_eq_1_falls_back_silently(tmp_path):
    """Bei n=1 mit analyze_all_clusters=True soll das Tool sauber als
    Single-Cluster-Lauf durchlaufen (keine Endlosrekursion, kein Crash).

    Wir testen nur den Fallback-Pfad ohne den vollen Lauf -- die
    Single-Cluster-Logik wird durch monkeypatch unterbrochen.
    """
    if not TEST_DATA.exists():
        pytest.skip("Test-Logfile fehlt")

    from modenanalyse_2fe2s import runner as _runner

    # Logfile entpacken (das Tool kann xz nicht direkt lesen)
    log_path = tmp_path / "test.log"
    with lzma.open(TEST_DATA, "rt") as f, open(log_path, "w") as out:
        out.write(f.read())

    cfg = Config(
        log_file=str(log_path),
        output_dir=str(tmp_path / "out"),
        analyze_all_clusters=True,
    )

    # Stub fuer _run_analysis_single, damit der Test nicht 2 min braucht
    called = []

    def _stub(c):
        called.append(c)
        return 0

    orig = _runner._run_analysis_single
    _runner._run_analysis_single = _stub
    try:
        rc = _runner.run_analysis(cfg)
    finally:
        _runner._run_analysis_single = orig

    assert rc == 0
    # Bei n=1: Stub wird genau einmal mit dem Original-cfg aufgerufen
    assert len(called) == 1
    assert called[0] is cfg


def test_multi_cluster_wrapper_n_eq_2(tmp_path, monkeypatch):
    """Wrapper-Integrationstest: n=2 Cluster -> 2 Subordner + Summary.

    Mockt find_all_clusters (gibt 2 Cluster) und _run_analysis_single
    (zaehlt nur Aufrufe). Verifiziert:
      - 2 Aufrufe an _run_analysis_single
      - Subordner cluster_0/, cluster_1/
      - analyze_all_clusters=False in den rekursiven Aufrufen
        (Endlosrekursions-Schutz)
      - multi_cluster_summary.txt wird geschrieben
    """
    from modenanalyse_2fe2s import runner
    from modenanalyse_2fe2s import geometry as _geom
    from modenanalyse_2fe2s import orca_io as _orca

    def mock_find_all(atoms, cfg):
        return [
            ([1, 2], [3, 4],
             {"fe_fe": 2.7, "fe_s_min": 2.2, "fe_s_max": 2.3}),
            ([5, 6], [7, 8],
             {"fe_fe": 2.74, "fe_s_min": 2.21, "fe_s_max": 2.31}),
        ]

    class _MockPR: pass
    monkeypatch.setattr(_geom, "find_all_clusters", mock_find_all)
    monkeypatch.setattr(_orca, "load_orca_hess", lambda p: _MockPR())
    monkeypatch.setattr(_orca, "parseresult_to_atoms",
                         lambda pr, include_hydrogen: ([], {}))

    calls = []
    def _stub(c):
        calls.append((c.cluster_index, c.output_dir,
                      c.analyze_all_clusters))
        return 0
    monkeypatch.setattr(runner, "_run_analysis_single", _stub)

    fake_hess = tmp_path / "fake.hess"
    fake_hess.touch()
    cfg = Config(
        log_file=str(fake_hess),
        output_dir=str(tmp_path / "out"),
        analyze_all_clusters=True,
    )

    rc = runner.run_analysis(cfg)
    assert rc == 0
    assert len(calls) == 2
    # Cluster-Indizes
    assert calls[0][0] == 0
    assert calls[1][0] == 1
    # Subordner-Konvention
    assert calls[0][1].endswith("cluster_0")
    assert calls[1][1].endswith("cluster_1")
    # Endlosrekursions-Schutz: rekursive Aufrufe haben analyze_all=False
    assert calls[0][2] is False
    assert calls[1][2] is False
    # Summary-Datei wurde geschrieben
    summary = tmp_path / "out" / "multi_cluster_summary.txt"
    assert summary.exists()
    text = summary.read_text(encoding="utf-8")
    assert "Cluster #0" in text
    assert "Cluster #1" in text
    assert "Number of clusters: 2" in text


def test_multi_cluster_wrapper_partial_failure(tmp_path, monkeypatch):
    """Wenn ein Cluster scheitert, soll der Wrapper den Rest fortsetzen
    und am Ende rc=1 zurueckgeben (overall_status)."""
    from modenanalyse_2fe2s import runner
    from modenanalyse_2fe2s import geometry as _geom
    from modenanalyse_2fe2s import orca_io as _orca

    def mock_find_all(atoms, cfg):
        return [
            ([1, 2], [3, 4],
             {"fe_fe": 2.7, "fe_s_min": 2.2, "fe_s_max": 2.3}),
            ([5, 6], [7, 8],
             {"fe_fe": 2.74, "fe_s_min": 2.21, "fe_s_max": 2.31}),
        ]

    class _MockPR: pass
    monkeypatch.setattr(_geom, "find_all_clusters", mock_find_all)
    monkeypatch.setattr(_orca, "load_orca_hess", lambda p: _MockPR())
    monkeypatch.setattr(_orca, "parseresult_to_atoms",
                         lambda pr, include_hydrogen: ([], {}))

    def _stub(c):
        if c.cluster_index == 0:
            raise RuntimeError("Simulierter Fehler")
        return 0
    monkeypatch.setattr(runner, "_run_analysis_single", _stub)

    fake_hess = tmp_path / "fake.hess"
    fake_hess.touch()
    cfg = Config(
        log_file=str(fake_hess),
        output_dir=str(tmp_path / "out"),
        analyze_all_clusters=True,
    )

    rc = runner.run_analysis(cfg)
    # Cluster 0 scheitert, Cluster 1 ist OK -> overall rc=1
    assert rc == 1
    summary = (tmp_path / "out" / "multi_cluster_summary.txt").read_text(
        encoding="utf-8")
    assert "EXCEPTION" in summary
    assert "Status:    OK" in summary
