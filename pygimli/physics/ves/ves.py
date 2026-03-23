import numpy as np
from .vesModelling import VESModelling

def simulate(res, thk=None, **kwargs):
    """Synthetic forward simulation.

    Parameters
    ----------
    res : array-like
        resistivity vector or thicknesses+resistivities
    thk : array-like
        thickness vector (if not, taken from res entries)
        geometry definition
    ab2 : array-like
        Half distance between the current electrodes A and B.
    mn2 : array-like
        Half distance between the potential electrodes M and N.
        OR
    am : array-like
        Part of data basis. Distances between A and M electrodes.
    bm : array-like
        Part of data basis. Distances between B and M electrodes.
    an : array-like
        Part of data basis. Distances between A and N electrodes.
    bn : array-like
        Part of data basis. Distances between B and N electrodes.
        OR
    data : pg.DataContainerERT
        pyGIMLi data container from which AM/AN/BM,BN are taken
    """
    if thk is None:  # assume d1, .., rho1,
        nLay = (len(res)-1) // 2 + 1
        thk = res[:nLay-1]
        res = res[nLay-1:2*nLay-1]

    assert len(res) == len(thk)+1, "Lengths of res/thk+1 must match!"
    fop = VESModelling(**kwargs)
    return fop.response(np.concatenate([thk, res]))
