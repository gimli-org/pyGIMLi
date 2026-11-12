#!/usr/bin/env python
"""Modelling classes for magnetic resonance sounding."""

# general modules to import according to standards
import pygimli as pg
import numpy as np


class MRS1dBlockQTModelling(pg.core.ModellingBase):
    """
    MRS1dBlockQTModelling - pygimli modelling class for block-mono QT inversion.

    f=MRS1dBlockQTModelling(lay, KR, KI, zvec, t, verbose = False )
    """

    def __init__(self, nlay, K, zvec, t, verbose=False):
        """Set up with number of layers, kernel, z and t vectors."""
        mesh = pg.meshtools.createMesh1DBlock(nlay, 2)  # thk, wc, T2*
        pg.core.ModellingBase.__init__(self, mesh)
        self.K_ = K
        self.zv_ = np.array(zvec)
        self.nl_ = nlay
        self.nq_ = len(K)
        self.t_ = np.array(t)
        self.nt_ = len(t)

    def response(self, par):
        """Yield model response cube as vector."""
        nl = self.nl_
        thk = par[0 : nl - 1]  # (0, nl - 1)
        wc = par[nl - 1 : 2 * nl - 1]  # (nl - 1, 2 * nl - 1)
        t2 = par[2 * nl - 1 : 3 * nl - 1]  # par(2 * nl - 1, 3 * nl - 1)
        zthk = np.cumsum(thk)
        zv = self.zv_
        lzv = len(zv)
        izvec = np.zeros(nl + 1, np.int32)
        rzvec = np.zeros(nl + 1)
        for i in range(nl - 1):
            ii = (zv < zthk[i]).argmin()
            izvec[i + 1] = ii
            if ii <= len(zv):
                rzvec[i + 1] = (zthk[i] - zv[ii - 1]) / (zv[ii] - zv[ii - 1])

        izvec[-1] = lzv - 1
        A = np.zeros((self.nq_, self.nt_), dtype=complex)
        for i in range(nl):
            wcvec = np.zeros(lzv - 1)
            wcvec[izvec[i] : izvec[i + 1]] = wc[i]
            if izvec[i + 1] < lzv:
                wcvec[izvec[i + 1] - 1] = wc[i] * rzvec[i + 1]
            if izvec[i] > 0:
                wcvec[izvec[i] - 1] = wc[i] * (1 - rzvec[i])
            amps = np.dot(self.K_, wcvec)
            for ii, a in enumerate(A):
                a += np.exp(-self.t_ / t2[i]) * amps[ii]

        return np.abs(A).ravel()  # formerly pg as vector


if __name__ == "__main__":
    pass
