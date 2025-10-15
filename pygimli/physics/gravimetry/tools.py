"""Tools for gravity and magnetics."""
import numpy as np
import pygimli as pg


def depthWeighting(mesh, z0:float=25, height:float=0, power:float=1.5, normalize:bool=False, cell:bool=False):
    """Return Li&Oldenburg like depth weighting of boundaries or cells.

    To account for inherent ambiguity in potential field methods, a depth
    weighting of the regulariation term is often applied in inversion,
    going back to Li & Oldenburg (1996, eq. 18):

    .. math::
        w(z) = \frac{1}{(z+z_0)^{3/2}}
    
    Here, we reformulate the original equation in order to make it unitless:

    .. math::
        w(z) = \frac{1}{(z/z_0+1)^p}

    Which is the same as before except a factor that can go into the 
    regularization strength. The function describes a decay with depth,
    starting from w(0)=1 and approaching w(z)=0 for z>>z0.
    z is the depth of a boundary or cell below the sensor and is determined by
    its mesh center and the sensor position.
    
    Parameters
    ----------
    z0 : float [25]
        centroid depth (w(z0)=1/2^p)
    height : float
        sensor height to be added to depth
    power : float [1.5]
        exponent
    normalize : bool [False]
        normalize such that the mean weight is 1 (like unweighted)
    cell : bool [True]
        use cell center (for cType=0|2|geostat) instead of boundary center (cType=1)

    Returns
    -------
    depth weighting vector for all interior boundaries or all cells 
    """
    # if cell:
    #     z = np.abs(pg.z(mesh.cellCenter()))
    #     if mesh.dim() == 2 and max(z) == min(z):
    #         z = np.abs(pg.y(mesh.cellCenter()))
    # else:
    #     z = np.abs(pg.z(mesh.innerBoundaryCenters()))
    #     if mesh.dim() == 2 and max(z) == min(z):
    #         z = np.abs(pg.y(mesh.innerBoundaryCenters()))
    if cell:
        cc = mesh.cellCenter()
    else:
        cc = mesh.innerBoundaryCenters()

    z = pg.z(cc)
    if mesh.dim() == 2 and np.isclose(max(z), min(z)):
        z = pg.y(cc)

    weight = 2 / (np.abs(height-z)/z0 + 1)**power 
    if normalize:
        weight /= np.median(weight)  # normalize that maximum is 1

    return weight
