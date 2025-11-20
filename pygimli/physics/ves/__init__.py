"""
Vertical Electric Sounding (VES).

Tools, modelling operators, and managers for VES.

Main class:

* VESManager - Manager for loading, inverting & plotting

Modelling operators:

* VESModelling - DC modelling operator using block model
* VESRhoModelling - DC modelling operator
* VESCModelling - DC/IP (complex-valued) modelling
"""

from .vesManager import VESManager
from .vesModelling import VESModelling, VESCModelling, VESRhoModelling

Manager = VESManager
