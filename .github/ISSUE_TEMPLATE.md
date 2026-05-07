---
name: Bug report or feature request
about: Help us improve modenanalyse_2fe2s
---

## Type
- [ ] Bug report
- [ ] Feature request
- [ ] Documentation issue
- [ ] Question

## Description
Describe the issue or request clearly.

## For bugs

**Expected behavior:**
What you expected to happen.

**Actual behavior:**
What actually happened. Include error messages and tracebacks if any.

**Reproduction steps:**
1. Configure with `Config(...)` like this: ...
2. Run `run_analysis(cfg)`
3. See error

**Environment:**
- OS: (e.g., Windows 11, macOS 14, Ubuntu 22.04)
- Python version: (`python --version`)
- Package version: (`pip list | grep modenanalyse`)
- Dependencies of interest: (e.g., `pip list | grep -E "scipy|numpy|umap"`)

**Input data:**
- DFT backend: Gaussian-16 / ORCA-6
- File type: `.log` / `.hess`
- Cluster type: [2Fe-2S] Cys₄ / Rieske / mitoNEET / dimer / other: ___
- System size: ~___ atoms, ~___ modes
- (If possible) attach a minimal test file that reproduces the issue.

## For feature requests

**Use case:**
What scientific problem does this feature solve?

**Proposed behavior:**
Sketch how the feature should work from a user's perspective.

**Alternatives considered:**
Any workarounds you've tried.
