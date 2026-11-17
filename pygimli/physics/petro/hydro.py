#!/usr/bin/env python
"""TODO WRITEME."""

import pygimli as pg


def permeabilityEngelhardtPitter(poro, q=3.5, s=5e-3, mesh=None, meshI=None):
    r"""Empirical model for porosity to hydraulic permeability.

    Postulated for sand and sandstones. :cite:`EngelhardtPit1955`

    .. math::
        k & = 2\cdot 10^7 \frac{\phi^2}{(1-\phi)^2}* \frac{1}{S^2} \\
        S & = q\cdot s \\
        s & = \sum_{i=1}(\frac{P_i}{r_i})

    * :math:`\phi` - poro 0.0 --1.0
    * :math:`q` - (3 for spheres, > 3 shape differ from sphere)
        3.5 sand
    * :math:`s` - in cm^-1 (s = 1/r for particles with homogeneous radii r)
    * :math:`P_i` - Particle ration with radii :math:`r_i` on 1cm^3 Sample
    * :math:`S` - in cm^-1  specific surface in cm^2/cm^3

    Parameters
    ----------
    poro : float | array | callable
        Porosity between 0.0 and 1.0.
    q : float | array | callable [3.5]
        Shape factor (3 for spheres, > 3 for non-spherical grains).
    s : float | array | callable [5e-3]
        Specific surface in cm^-1.
    mesh : :gimliapi:`GIMLI::Mesh` [None]
        Mesh to calculate a value for every cell.
    meshI : :gimliapi:`GIMLI::Mesh` [None]
        Mesh to interpolate the result to.

    Returns
    -------
    k : float | array
        Hydraulic permeability in Darcy.
    """
    poro = pg.solver.parseArgToArray(poro, mesh.cellCount(), mesh)
    q = pg.solver.parseArgToArray(q, mesh.cellCount(), mesh)
    s = pg.solver.parseArgToArray(s, mesh.cellCount(), mesh)

    st = q * s
    k = (
        2e7
        * (poro**2 / (1.0 - poro) ** 2)
        * 1.0
        / st**2
        * pg.physics.Constants.Darcy
    )

    if meshI is not None:
        k = pg.interpolate(mesh, k, meshI.cellCenters())
        k = pg.meshtools.fillEmptyToCellArray(meshI, k)
    return k


if __name__ == "__main__":
    pass
