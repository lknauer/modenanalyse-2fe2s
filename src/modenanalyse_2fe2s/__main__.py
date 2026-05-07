# Copyright (C) 2026 Lukas Knauer, AG Schünemann, RPTU Kaiserslautern-Landau
# SPDX-License-Identifier: GPL-3.0-or-later
# Teil von modenanalyse_2fe2s — see LICENSE in the Wurzelverzeichnis.

# -*- coding: utf-8 -*-
"""Erlaubt ``python -m modenanalyse_2fe2s run.toml`` als Aequivalent zur
``modenanalyse-2fe2s``-CLI."""
import sys
from .cli import main

sys.exit(main())
