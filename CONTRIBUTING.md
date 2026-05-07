# Contributing to modenanalyse_2fe2s

Thank you for your interest in `modenanalyse_2fe2s`! This document
describes the conventions for code, test, and documentation
contributions.

## Bug reports and feature requests

Please use GitHub Issues:
- **Bug:** Step-by-step reproduction, expected vs. observed behavior,
  environment (operating system, Python version,
  `pip list | grep modenanalyse`).
- **Feature:** Use case + why the tool currently can't do this.

## Pull requests

1. Fork the repo, create a feature branch:
   ```
   git checkout -b feat/<short-description>
   ```

2. Write code — please note:
   - **Code style:** PEP 8, max. 100 characters per line.
   - **Docstrings:** NumPy style for all public functions
     (`Parameters`, `Returns`, `Notes`, `Examples`).
   - **Comments:** English preferred; code identifiers (variables,
     functions, classes) always in English.
   - **Constants** in UPPERCASE.

3. Write tests:
   - Every new code path needs a test in `tests/`.
   - Fast unit tests in `test_<module>.py`.
   - Long integration tests with `@pytest.mark.slow` (skipped by default).

4. Run the tests:
   ```
   pytest tests/                  # fast tests (~1 sec)
   pytest tests/ -m slow          # end-to-end smoke tests (~1 min)
   ```

5. Update CHANGELOG:
   - New entry at the top of `CHANGELOG.md` under `## [Unreleased]`.
   - Format as existing entries.

6. Submit pull request with:
   - Meaningful title (e.g., `feat: HA-Reorg for Tyr ligands`).
   - Description of what changed and why.
   - Reference to issue if applicable (`Fixes #42`).

## Version numbers

`modenanalyse_2fe2s` follows [Semantic Versioning](https://semver.org/):
- **MAJOR** (x.0.0): Breaking changes in API or behavior.
- **MINOR** (1.x.0): New features, backward-compatible.
- **PATCH** (1.0.x): Bug fixes, polish, backward-compatible.

## Scientific changes

For changes that affect numerical results (e.g., new reorganization
conventions, different cluster detection, default changes):

- Document the reasoning in CHANGELOG in detail.
- Provide a test log file with the new behavior + before/after comparison.
- List in the manual under the relevant section with `\paragraph{vX.Y.Z}`.

## Repository structure

```
modenanalyse_2fe2s/
├── src/modenanalyse_2fe2s/    # main source code
│   ├── config.py              # Config dataclass
│   ├── runner.py              # run_analysis()
│   ├── core.py                # per-mode analysis
│   ├── reorganization.py      # Marcus-Hush reorganization pipeline
│   ├── pcet_et.py             # geometry adapter
│   ├── geometry.py            # cluster detection
│   ├── embedding.py           # UMAP + HDBSCAN
│   ├── export.py              # Excel writer
│   ├── logio.py               # Gaussian log reader
│   ├── orca_io.py             # ORCA .hess reader
│   └── cli.py                 # command line
├── tests/                     # pytest tests
│   ├── data/                  # test log files (xz-compressed)
│   └── test_*.py
├── docs/                      # manual (LaTeX + PDF)
├── examples/                  # TOML templates
├── CHANGELOG.md               # release notes
├── README.md                  # project overview
├── CITATION.cff               # software citation format
├── LICENSE                    # GPL-3.0
└── pyproject.toml             # build configuration
```

## License

By submitting a pull request, you agree that your contribution will be
released under the GPL-3.0-or-later license.

## Contact

- **Lukas Knauer** — main developer
- **AG Schünemann**, Department of Physics, RPTU Kaiserslautern-Landau
