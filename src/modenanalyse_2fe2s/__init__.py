# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Part of modenanalyse_2fe2s -- see LICENSE in repository root.

# -*- coding: utf-8 -*-
"""
modenanalyse_2fe2s -- Quantum-chemical normal-mode analysis of
[2Fe-2S] iron-sulfur proteins from Gaussian-16 or ORCA-6 frequency
calculations.

Public API
----------

  * ``Config``           : Full configuration dataclass.
  * ``run_analysis(cfg)``: Programmatic entry point.
  * ``__version__``      : Package version (1.0.x).

Usage example
-------------

>>> from modenanalyse_2fe2s import Config, run_analysis
>>> cfg = Config(
...     log_file   = r"D:\\data\\dimer.log",
...     pdb_file   = r"D:\\data\\dimer.pdb",
...     output_dir = r"D:\\data\\results",
...     temp_k     = 40.0,
... )
>>> run_analysis(cfg)

Multi-cluster systems
---------------------

For dimers and multi-[2Fe-2S] systems, set ``analyze_all_clusters=True``
to analyze all clusters automatically into separate subfolders:

>>> cfg = Config(
...     log_file = "dimer.hess",
...     pdb_file = "dimer.pdb",
...     output_dir = "results",
...     analyze_all_clusters = True,
... )
>>> run_analysis(cfg)

Command line
------------

::

    modenanalyse-2fe2s --log-file dimer.log --pdb dimer.pdb \\
                       --output-dir results --temp-k 40

For full documentation see ``docs/QUICKSTART.md`` and ``docs/Manual.pdf``
(English overview), or ``docs/Manual_full_EN.pdf`` (complete English
translation, ~103 pages).
"""

from .config import Config, __version__
from .runner import run_analysis

__all__ = [
    "Config",
    "run_analysis",
    "__version__",
]
