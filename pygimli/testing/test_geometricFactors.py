from math import sqrt
import numpy as np
import pygimli as pg
from pygimli.physics import ert

def geomFact(xABMN, yABMN=None, zABMN=None):
    """Compute geometric factor by hand."""
    if yABMN is None:
        yABMN = np.zeros_like(xABMN)
    if zABMN is None:
        zABMN = np.zeros_like(xABMN)
    dd1 = np.zeros(4)
    dd2 = np.zeros(4)
    # AM, AN, BM, BN
    for i, cp in enumerate([[0, 2], [0, 3], [1, 2], [1, 3]]):
        c, p = cp
        dd1[i] = sqrt((xABMN[c]-xABMN[p])**2 +
                     (yABMN[c]-yABMN[p])**2 +
                     (zABMN[c]-zABMN[p])**2)
        dd2[i] = sqrt((xABMN[c]-xABMN[p])**2 +
                     (yABMN[c]-yABMN[p])**2 +
                     (-zABMN[c]-zABMN[p])**2)
    dI = 1/dd1 + 1/dd2
    return 4*np.pi / (dI[0]-dI[1]-dI[2]+dI[3])
    
def check(data):
    """Check whether both values are matching."""
    k1 = ert.createGeometricFactors(data)
    k2 = geomFact(pg.x(data), pg.y(data), pg.z(data))
    print(k1, k2)
    np.testing.assert_allclose(k1, k2)


if __name__ == '__main__':
    data = ert.DataContainer()
    for i in range(4):
        data.createSensor(pg.Pos(i, 0, 0))

    data.createFourPointData(0, 0, 1, 2, 3)
    data.createFourPointData(1, 0, 1, 2, 3)
    
    # surface 2D
    check(data)
    
    # surface 3D
    for i in range(4):
        data.setSensor(i, pg.Pos(i*2, i, 0))

    check(data)
    
    # crosshole 2D
    data.setSensor(0, [1, -2])
    data.setSensor(1, [1, -3])
    data.setSensor(2, [3, -2])
    data.setSensor(3, [3, -3])
    check(data)
    
    # crosshole 3D
    data.setSensor(0, [1, 2, -2])
    data.setSensor(1, [1, 3, -3])
    data.setSensor(2, [3, 4, -2])
    data.setSensor(3, [3, 5, -3])
    check(data)