# coding=utf-8
"""
Tools to assess the quality of unstructured meshes.
---------------------------------------------------

Quality measures are based on the following review article:

Field, D. A. (2000), Qualitative measures for initial meshes. Int. J. Numer.
Meth. Engng., 47: 887–906.
"""

import numpy as np


# Helper functions
def _boundaryLengths(cell):
    """Return boundary lengths of a given cell."""
    return np.array(
        [cell.boundary(i).size() for i in range(cell.boundaryCount())])


def _unitVector(vector):
    """Return the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def _angleBetween(v1, v2):
    """Return the angle (in degrees) between vectors v1 and v2.

    Examples
    --------
    >>> print(_angleBetween((1, 0, 0), (1, 0, 0)))
    0.0
    >>> print(_angleBetween((1, 0, 0), (-1, 0, 0)))
    180.0
    """
    v1_u = _unitVector(v1)
    v2_u = _unitVector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    return np.degrees(angle)


def _cellAngles(cell):
    """Return angles of a triangular cell.

    Examples
    --------
    >>> import pygimli as pg
    >>> mesh = pg.Mesh()
    >>> for pos in (0.,0.), (0.,1.), (1.,0.):
    ...     n = mesh.createNode(pos[0], pos[1], 0.0)
    >>> cell = mesh.createCell((0, 1, 2))
    >>> print(_cellAngles(cell)[0])
    90.0
    """
    if cell.nodeCount() > 3:
        raise Exception("Cell %d is not a triangular cell." % cell.id())

    a = cell.node(0).pos()
    b = cell.node(1).pos()
    c = cell.node(2).pos()

    ab = b - a
    ac = c - a
    bc = c - b

    alpha = _angleBetween(ab, ac)
    beta = _angleBetween(-ab, bc)
    gamma = _angleBetween(-ac, -bc)
    assert np.allclose(gamma, 180.0 - alpha - beta)

    return alpha, beta, gamma


# Quality measures
def eta(cell):
    r"""Return default triangle quality (eta) of a given cell.

    The quality measure relates the area of the triangle (a)
    to its edge lengths (l1, l2, l3).

    .. math::

        \eta = \frac{4\sqrt{3}a}{l_1^2 + l_2^2 + l_3^2}
    """
    return 4 * np.sqrt(3) * cell.size() / np.sum(_boundaryLengths(cell)**2)


def minimumAngle(cell):
    """Return the normalized minimum angle of a given cell."""
    return np.min(_cellAngles(cell)) / 60.


def nsr(cell):
    r"""Return the normalized shape ratio (NSR) for a given cell.

    Also referred to as the radius ratio, as it is described by
    the ratio between the inradius (r) and the circumradius (R).

    .. math::

        \rho = \frac{2r}{R}
    """
    a, b, c = _boundaryLengths(cell)
    r = 2 * cell.size() / (a + b + c)  # inradius
    R = 0.25 * a * b * c / cell.size()  # circumradius

    return 2 * r / R


# Main function
def quality(mesh, measure="eta"):
    """Return the quality of a given triangular mesh.

    Parameters
    ----------
    mesh : mesh object
        Mesh for which the quality is calculated.
    measure : quality measure, str
        Can be either "eta", "nsr", or "minimumAngle".

    Examples
    --------
    >>> # no need to import matplotlib
    >>> import pygimli as pg
    >>> from pygimli.meshtools import polytools as plc
    >>> from pygimli.meshtools import quality
    >>> # Create Mesh
    >>> world = plc.createWorld(start=[-10, 0], end=[10, -10],
    ...                         marker=1, worldMarker=False)
    >>> c1 = plc.createCircle(pos=[0.0, -5.0], radius=3.0, area=.3)
    >>> mesh = pg.meshtools.createMesh([world, c1], quality=21.3)
    >>> # Compute and show quality
    >>> q = quality(mesh, measure="nsr")
    >>> ax, _ = pg.show(mesh, q, cMap="RdYlGn", showMesh=True, cMin=0.5,
    ...                 cMax=1.0, label="Normalized shape ratio")

    See also
    --------
    eta, nsr, minimumAngle
    """

    # Implemented quality measures (Triangular meshses only.)
    measures = {"eta": eta, "nsr": nsr, "minimumAngle": minimumAngle}

    m = measures[measure]
    qualities = [m(cell) for cell in mesh.cells()]
    return qualities
