#!/usr/bin/env python
"""Direct current electromagnetics.

This package contains tools, modelling operators, and managers for
Vertical Electric Sounding (VES).
"""

import pygimli as pg
from .vesManager import VESManager
from .vesModelling import VESModelling, VESCModelling, VESRhoModelling

Manager = VESManager
