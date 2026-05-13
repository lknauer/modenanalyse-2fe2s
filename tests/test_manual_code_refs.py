# Copyright (C) 2026 Lukas Knauer, AG Schuenemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
"""Verify that all \\codref{...} entries in Manual.tex point to
functions/classes that actually exist in the source tree.

Purpose
-------
Section 2 of the manual links each non-trivial formula to its
implementation via a small LaTeX macro,

    \\codref{module.function}    →   "core.compute_thermal_amplitude"

If a developer renames or removes such a function in the source
tree without updating the manual, this test fails. The manual
therefore cannot drift silently out of sync with the code.

Run
---
    python3 -m pytest tests/test_manual_code_refs.py -v

What this test does NOT do
--------------------------
- It does not verify that the formula is mathematically equivalent
  to what the function computes. That's a human review task (and
  is exactly the kind of audit we documented in CHANGELOG.md
  under the v1.0.4 entry).
- It does not check the German Anleitung.tex (which mirrors
  Manual.tex 1:1; if a \\codref check there is desired, add it
  analogously).
"""
from __future__ import annotations
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


# Regex strips a single \codref{...} entry. Inside, references can be
# comma-separated, e.g. \codref{core.analyze_his_hn, core.analyze_fe_ligand}.
# The LaTeX file uses backslash-escaped underscores; we strip those.
_CODREF_RE = re.compile(r"\\codref\{([^}]+)\}")


def _extract_codrefs(tex_path: Path) -> list[tuple[str, str]]:
    """Return list of (module, function_or_class) pairs from the tex file.

    Multiple comma-separated refs in one \\codref are split out.
    LaTeX escapes (backslash-underscore) are unescaped.

    Lines that start (after optional whitespace) with `%` are LaTeX
    comments and are skipped -- they may contain example codrefs
    documenting the macro itself (e.g. "\\codref{module.function}"
    in the macro definition's comment block).
    """
    refs: list[tuple[str, str]] = []
    with tex_path.open(encoding="utf-8") as f:
        for line in f:
            # Skip LaTeX comment lines (the macro's own definition comment
            # contains an example \codref that is not a real reference).
            stripped = line.lstrip()
            if stripped.startswith("%"):
                continue
            # Strip any inline comment ('%' that is not preceded by '\')
            # before applying the regex, so an inline example like
            # "...  % see \codref{module.function} for syntax" is ignored.
            inline_comment = re.search(r"(?<!\\)%", line)
            if inline_comment:
                line = line[:inline_comment.start()]
            for match in _CODREF_RE.findall(line):
                # Each match may contain several comma-separated refs
                for item in match.split(","):
                    item = item.strip().replace(r"\_", "_")
                    if not item:
                        continue
                    if "." not in item:
                        pytest.fail(
                            f"\\codref entry has no module-dot-name "
                            f"format: {item!r} (in {tex_path.name})")
                    module, name = item.rsplit(".", 1)
                    refs.append((module.strip(), name.strip()))
    return refs


def _function_exists(module: str, name: str) -> bool:
    """Check if `def name(` or `class name` exists in
    src/modenanalyse_2fe2s/{module}.py."""
    src = ROOT / "src" / "modenanalyse_2fe2s" / f"{module}.py"
    if not src.exists():
        return False
    txt = src.read_text(encoding="utf-8")
    patterns = [
        rf"^def {re.escape(name)}\(",
        rf"^async def {re.escape(name)}\(",
        rf"^class {re.escape(name)}\b",
        # Private helpers used as the canonical implementation:
        rf"^def _{re.escape(name.lstrip('_'))}\(",
    ]
    for pat in patterns:
        if re.search(pat, txt, re.MULTILINE):
            return True
    return False


def test_manual_en_codrefs_exist():
    """Every \\codref{...} entry in Manual.tex must reference a
    function or class that exists in the source tree.
    """
    tex_path = ROOT / "docs" / "Manual.tex"
    assert tex_path.exists(), f"Manual.tex not found at {tex_path}"

    refs = _extract_codrefs(tex_path)
    assert refs, (
        "No \\codref entries found in Manual.tex. Either someone "
        "deleted them all, or the regex broke -- check this test.")

    missing: list[tuple[str, str]] = []
    for module, name in refs:
        if not _function_exists(module, name):
            missing.append((module, name))

    if missing:
        lines = [f"  - {m}.{n}" for m, n in missing]
        pytest.fail(
            f"{len(missing)} \\codref entries in Manual.tex point to "
            f"functions that no longer exist in the source tree. Either "
            f"the code was refactored without updating the manual, or "
            f"the manual contains a typo:\n" + "\n".join(lines))


def test_manual_en_codrefs_nonempty():
    """Sanity: the manual should contain at least 8 \\codref entries
    (one per formula in the theoretical-foundations section)."""
    tex_path = ROOT / "docs" / "Manual.tex"
    refs = _extract_codrefs(tex_path)
    assert len(refs) >= 8, (
        f"Only {len(refs)} \\codref entries found in Manual.tex; "
        f"expected at least 8 for the core formulas (u_rms, OOP, "
        f"dr-classical, dr-HA, M_X, dr-RSS, bend-split, sigma-therm). "
        f"Has the manual lost its formula-to-code annotations?")


def test_manual_en_codrefs_use_known_modules():
    """\\codref entries must point to modules under
    src/modenanalyse_2fe2s/. A typo like 'cores.foo' (instead of
    'core.foo') would otherwise show up only as a missing-function
    failure, which is less helpful."""
    known_modules = {
        p.stem for p in (ROOT / "src" / "modenanalyse_2fe2s").glob("*.py")
        if p.stem != "__init__"}
    tex_path = ROOT / "docs" / "Manual.tex"
    refs = _extract_codrefs(tex_path)
    unknown = sorted({m for m, _ in refs if m not in known_modules})
    if unknown:
        lines = [f"  - {m!r}" for m in unknown]
        pytest.fail(
            f"\\codref in Manual.tex names {len(unknown)} unknown "
            f"module(s). Either a typo or a missing source file:\n"
            + "\n".join(lines) + "\n\n"
            f"Known modules: {sorted(known_modules)}")


# =============================================================================
# Supplement codref tests
# =============================================================================
# The Supplement.tex carries far more \codref entries than the Manual
# (one per derived formula, plus one per algorithm section). Each must
# point to an existing function in the source tree, just like the Manual.

def test_supplement_codrefs_exist():
    """Every \\codref{...} entry in Supplement.tex must reference a
    function or class that exists in the source tree. Same property
    as for Manual.tex, applied to the technical supplement.
    """
    tex_path = ROOT / "docs" / "Supplement.tex"
    if not tex_path.exists():
        pytest.skip("Supplement.tex not present (not yet written).")

    refs = _extract_codrefs(tex_path)
    assert refs, (
        "Supplement.tex exists but contains no \\codref entries. "
        "Either the regex broke or the codrefs have been removed.")

    missing: list[tuple[str, str]] = []
    for module, name in refs:
        if not _function_exists(module, name):
            missing.append((module, name))

    if missing:
        lines = [f"  - {m}.{n}" for m, n in missing]
        pytest.fail(
            f"{len(missing)} \\codref entries in Supplement.tex point "
            f"to functions that no longer exist in the source tree. "
            f"Either the code was refactored without updating the "
            f"supplement, or the supplement contains a typo:\n"
            + "\n".join(lines))


def test_supplement_codrefs_use_known_modules():
    """All \\codref entries in Supplement.tex must reference real
    modules in src/modenanalyse_2fe2s/."""
    tex_path = ROOT / "docs" / "Supplement.tex"
    if not tex_path.exists():
        pytest.skip("Supplement.tex not present.")

    known_modules = {
        p.stem for p in (ROOT / "src" / "modenanalyse_2fe2s").glob("*.py")
        if p.stem != "__init__"}
    refs = _extract_codrefs(tex_path)
    unknown = sorted({m for m, _ in refs if m not in known_modules})
    if unknown:
        lines = [f"  - {m!r}" for m in unknown]
        pytest.fail(
            f"\\codref in Supplement.tex names {len(unknown)} unknown "
            f"module(s). Either a typo or a missing source file:\n"
            + "\n".join(lines) + "\n\n"
            f"Known modules: {sorted(known_modules)}")


# ========================================================================
# v1.0.4 post-Apd1-audit: \fn{} reference verification
# ========================================================================
# In addition to \codref{} tags (which are the primary formula→code
# annotation mechanism), the Supplement also uses \fn{module.function}
# inline references in worked examples, pseudo-code captions, and
# narrative text. These references are NOT covered by the \codref tests
# above, but the same drift problem applies: a renamed function leaves
# a phantom inline reference.
#
# The test below extracts every \fn{module.function} occurrence (with
# proper handling of multi-line LaTeX arguments) and verifies each
# refers to a real function/class/constant in the source tree.
# ========================================================================


def _extract_fn_refs(tex_path: Path) -> list[tuple[str, str]]:
    """Return list of (module, function_or_class) pairs from \\fn{} tags
    in the tex file.

    Only module.function style refs (single dot, module is a known
    .py file in src/) are extracted. Other \\fn{} uses (bare names,
    dataclass-field attributes like ChannelGeometry.r_donor, filenames
    like ``core.py``, historical illustrations marked with ``\\texttt``
    instead of ``\\fn``) are silently ignored — the focus is on real
    inline code references.
    """
    text = tex_path.read_text(encoding="utf-8")
    # Strip LaTeX comments (handle inline % comments too)
    lines = []
    for line in text.split("\n"):
        i = 0
        while i < len(line):
            if line[i] == "%" and (i == 0 or line[i-1] != "\\"):
                break
            i += 1
        lines.append(line[:i])
    text = "\n".join(lines)

    known_modules = {
        p.stem for p in (ROOT / "src" / "modenanalyse_2fe2s").glob("*.py")
        if p.stem != "__init__"
    }

    refs: list[tuple[str, str]] = []
    for m in re.finditer(r"\\fn\{", text):
        # Match the closing brace, supporting nested braces and newlines
        start = m.end()
        depth = 1
        i = start
        while i < len(text) and depth > 0:
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
            i += 1
        inner = text[start:i-1]
        # Normalize whitespace + newlines, strip arguments
        inner = re.sub(r"\s+", " ", inner).strip()
        ref = inner.replace(r"\_", "_").replace(r"\$", "$")
        ref = re.sub(r"\s*\(.*$", "", ref)
        # Only consider module.function style refs
        if ref.count(".") != 1:
            continue
        if ref.endswith(".py"):
            continue
        module, name = ref.split(".", 1)
        if module not in known_modules:
            continue
        refs.append((module, name))
    return refs


def _name_exists_in_module(module: str, name: str) -> bool:
    """Like _function_exists but also matches module-level constants."""
    src_file = ROOT / "src" / "modenanalyse_2fe2s" / f"{module}.py"
    if not src_file.exists():
        return False
    text = src_file.read_text(encoding="utf-8")
    patterns = [
        rf"^def {re.escape(name)}\b",
        rf"^async def {re.escape(name)}\b",
        rf"^class {re.escape(name)}\b",
        rf"^def _{re.escape(name)}\b",
        rf"^class _{re.escape(name)}\b",
        # Module-level assignments (constants like _CACHE_VERSION)
        rf"^{re.escape(name)}\s*[:=]",
    ]
    return any(re.search(p, text, re.MULTILINE) for p in patterns)


def test_supplement_fn_refs_exist():
    """Every inline \\fn{module.function} reference in Supplement.tex
    must point to a real function/class/constant in the source tree.

    Apd1-audit regression: this test was added after the post-release
    audit identified 8 phantom inline \\fn{} refs (e.g. invented
    function names like geometry.detect_mu_s_protonation, or wrong
    module assignments like logio.get_eigvec_orca instead of
    orca_io.get_eigvec_orca). The test catches both classes of error.
    """
    tex_path = ROOT / "docs" / "Supplement.tex"
    if not tex_path.exists():
        pytest.skip("Supplement.tex not present.")

    refs = _extract_fn_refs(tex_path)
    assert refs, (
        "Supplement.tex has no \\fn{module.function} refs. Either the "
        "Supplement was reformulated or the extractor is broken.")

    missing: list[tuple[str, str]] = []
    for module, name in refs:
        if not _name_exists_in_module(module, name):
            missing.append((module, name))

    if missing:
        lines = [f"  - \\fn{{{m}.{n}}}" for m, n in missing]
        pytest.fail(
            f"{len(missing)} inline \\fn{{module.function}} refs in "
            f"Supplement.tex point to names that don't exist in the "
            f"source tree:\n"
            + "\n".join(lines))
