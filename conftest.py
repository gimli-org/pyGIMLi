"""Pytest setup for headless test runs.

:func:`pygimli.test` prepares the current process for headless testing (Agg
backend, no interactive pyvista viewer), but pytest-xdist workers are separate
interpreter processes. This file is imported in every pytest process, so the
same setup is applied there.

An explicitly set ``MPLBACKEND`` (e.g. ``MPLBACKEND=TkAgg pytest``) wins, then
the workers keep the GUI backend for interactive debugging.
"""

import os

import pygimli as pg

if os.environ.get("MPLBACKEND", "Agg") == "Agg":
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

    # No interactive pyvista viewer while testing.
    pg.rc["pyvista.backend"] = None
