"""v1.1.0 SS -> SSE rename.

Guards two things: (1) legacy TOML keys (``analyze_ss``, ``ss_chain``) are
still accepted via backward-compatible aliases, and (2) the public API now
exposes the SSE names and no longer the old SS names.
"""
import os
import tempfile
import warnings

from modenanalyse_2fe2s.config import Config


def test_legacy_toml_keys_still_accepted():
    toml = 'log_file = "x.log"\nanalyze_ss = false\nss_chain = "B"\n'
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as fh:
        fh.write(toml)
        path = fh.name
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cfg = Config.from_toml(path)
        assert cfg.analyze_sse is False          # analyze_ss -> analyze_sse
        assert cfg.sse_chain == "B"              # ss_chain    -> sse_chain
        assert any("renamed to" in str(x.message) for x in w)
    finally:
        os.unlink(path)


def test_public_api_uses_sse_names():
    from modenanalyse_2fe2s import core, geometry, embedding
    for name in ("analyze_sse_element", "analyze_all_sse", "SSE_UMAP_METRICS"):
        assert hasattr(core, name), f"core.{name} missing"
    assert hasattr(geometry, "build_sse_center_map")
    assert hasattr(geometry, "build_sse_ca_center_map")
    assert hasattr(embedding, "compute_sse_umap_cluster")
    # the old SS names must be gone
    assert not hasattr(core, "analyze_ss_element")
    assert not hasattr(geometry, "build_ss_center_map")


def test_new_config_has_no_legacy_field_names():
    import dataclasses
    fields = {f.name for f in dataclasses.fields(Config)}
    assert "analyze_sse" in fields and "sse_chain" in fields
    assert "analyze_ss" not in fields and "ss_chain" not in fields
