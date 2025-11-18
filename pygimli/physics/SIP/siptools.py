#!/usr/bin/env python
"""pygimli functions for dc resistivity / SIP data."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.utils import rndig
from .importData import load

def astausgleich(ab2org, mn2org, rhoaorg):
    """Shift the branches of a dc sounding to generate a matching curve."""
    ab2 = np.asarray(ab2org)
    mn2 = np.asarray(mn2org)
    rhoa = np.asarray(rhoaorg)
    um = np.unique(mn2)
    for i in range(len(um) - 1):
        r0, r1 = [], []
        ac = np.intersect1d(ab2[mn2 == um[i]], ab2[mn2 == um[i + 1]])
        for a in ac:
            r0.append(rhoa[(ab2 == a) * (mn2 == um[i])][0])
            r1.append(rhoa[(ab2 == a) * (mn2 == um[i + 1])][0])

        if len(r0) > 0:
            fak = np.mean(np.array(r0) / np.array(r1))
            print(fak)
            if np.isfinite(fak) and fak > 0.:
                rhoa[mn2 == um[i + 1]] *= fak

    return rhoa  # formerly pg as vector


def loadSIPallData(filename, outnumpy=False):
    """Load SIP data with the columns ab/2,mn/2,rhoa and PHI.

    with the corresponding frequencies in the first row.
    """
    if outnumpy:
        A = np.loadtxt(filename)
        fr = A[0, 3:]
        ab2 = A[1:, 0]
        mn2 = A[1:, 1]
        rhoa = A[1:, 2]
        PHI = A[1:, 3:]
    else:
        A = pg.Matrix()
        pg.loadMatrixCol(A, 'sch/dc.ves')
        ndata = A.cols()
        ab2 = A[0](1, ndata)
        mn2 = A[1](1, ndata)
        rhoa = A[2](1, ndata)
        PHI = pg.Matrix()
        fr = []
        for i in range(3, A.rows()):
            fr.append(A[i][0])
            PHI.push_back(A[i](1, ndata))
    return ab2, mn2, rhoa, PHI, fr


def makeSlmData(ab2, mn2, rhoa=None, filename=None):
    """Generate a pygimli data container from ab/2 and mn/2 array."""
    data = pg.DataContainer()
    data.resize(len(ab2))
    pos = np.unique(np.hstack((ab2, mn2)))

    for elx in np.hstack((-pos[::-1], pos)):
        data.createElectrode(elx, 0., 0.)

    if filename is not None:
        with Path(filename).open('w') as f:
            f.write(str(len(pos) * 2) + '\n#x y z\n')
            for elx in np.hstack((-pos[::-1], pos)):
                f.write(str(elx) + '\t0\t0\n')
                f.write(str(len(ab2)) + '\n#a\tb\tm\tn\tk\trhoa\n')

    lpos = len(pos)
    iab = pos.searchsorted(ab2)
    imn = pos.searchsorted(mn2)

    if (filename is not None) & (rhoa is None):
        rhoa = np.ones(len(ab2))

    for i in range(len(iab)):
        # print -pos[iab[i]], -pos[imn[i]], pos[imn[i]], pos[iab[i]]
        k = (ab2[i]**2 - mn2[i]**2) * np.pi / mn2[i] / 2.0
        if filename is not None:
            f.write(str(lpos - iab[i]) + '\t' + str(lpos + iab[i] + 1) + '\t')
            f.write(str(lpos - imn[i]) + '\t' + str(lpos + imn[i] + 1) + '\t')
            f.write(str(rndig(k, 4)) + '\t' + str(rhoa[i]) + '\n')
        data.createFourPointData(i, int(lpos - iab[i]), int(lpos + iab[i] + 1),
                                 int(lpos - imn[i]), int(lpos + imn[i] + 1))

    if filename is not None:
        f.close()

    data.set('rhoa', pg.asvector(rhoa))
    return data


def showsounding(ab2, rhoa, resp=None, mn2=None, islog=True, xlab=None):
    """Display a sounding curve (rhoa over ab/2) and an additional response."""
    if xlab is None:
        xlab = r'$\rho_a$ in $\Omega$m'

    ab2a = np.asarray(ab2)
    rhoa = np.asarray(rhoa)
    fig, ax = plt.subplots()
    if mn2 is None:
        if islog:
            l1 = ax.loglog(rhoa, ab2, 'rx-', label='observed')
        else:
            l1 = ax.semilogy(rhoa, ab2, 'rx-', label='observed')

        if resp is not None:
            if islog:
                l2 = ax.loglog(resp, ab2, 'bo-', label='simulated')
            else:
                l2 = ax.semilogy(resp, ab2, 'bo-', label='simulated')

            ax.legend((l1, l2), ('obs', 'sim'), loc=0)
    else:
        for unmi in np.unique(mn2):
            if islog:
                l1 = ax.loglog(rhoa[mn2 == unmi], ab2a[mn2 == unmi],
                              'rx-', label='observed')
            else:
                l1 = ax.semilogy(rhoa[mn2 == unmi], ab2a[mn2 == unmi],
                                'rx-', label='observed')

            if resp is not None:
                l2 = ax.loglog(resp[mn2 == unmi], ab2a[mn2 == unmi],
                              'bo-', label='simulated')
                ax.legend((l1, l2), ('obs', 'sim'))

    # ax.axis('tight')
    ax.set_ylim((max(ab2), min(ab2)))
    locs = ax.set_yticks()[0]
    if len(locs) < 2:
        locs = np.hstack((min(ab2), locs, max(ab2)))
    else:
        locs[0] = max(locs[0], min(ab2))
        locs[-1] = min(locs[-1], max(ab2))

    a = []
    for l in locs:
        a.append(f'{rndig(l):g}')

    ax.set_yticks(locs, a)
    locs = ax.get_xticks()[0]
    a = []
    for l in locs:
        a.append(f'{rndig(l):g}')

    ax.set_xticks(locs, a)
    ax.grid(which='both')
    ax.set_xlabel(xlab)
    ax.set_ylabel('AB/2 in m')
    # plt.legend()
    return ax


def showsip1ddata(PHI, fr, ab2, mn2=None, cmax=None, ylab=True, cbar=True):
    """Display SIP phase data as image plot."""
    fig, ax = plt.subplots()
    pal = plt.cm.get_cmap()
    pal.set_under('w')
    pal.set_bad('w')
    if isinstance(PHI, pg.Vector):
        PHI = np.asarray(PHI)

    im = plt.imshow(PHI.reshape((len(ab2), len(fr))),
                  interpolation='nearest', cmap=pal)
    if cmax is None:
        cmax = np.max(PHI)

    im.set_clim((0., cmax))

    ax.xaxis.set_label_position('top')
    ax.set_xlabel('f in Hz')

    a = []
    df = 1
    for f in fr[::df]:
        a.append(f"{rndig(f):g}")

    ax.set_xticks(np.arange(0, len(fr), df), a)
    xtl = ax.get_xticklabels()
    for xtli in xtl:
        xtli.set_rotation('vertical')

    if ylab:
        a = []
        yla = 'AB/2'
        if mn2 is None:
            for i in range(len(ab2)):
                a.append(str(ab2[i]))
        else:
            yla = yla + '-MN/2'
            for i in range(len(ab2)):
                a.append(f'{rndig(ab2[i]):g}{rndig(mn2[i]):g}')

        ax.set_yticks(np.arange(len(ab2)), a)
        ax.set_ylabel(yla + ' in m')

    if cbar:
        fig.colorbar(aspect=40, shrink=0.6)

    plt.ylim((len(ab2) - 0.5, -0.5))
    plt.show()
    plt.ylim((len(ab2) - 0.5, -0.5))
    return


def showsip1dmodel(M, tau, thk, res=None, z=None,
                   cmin=None, cmax=None, islog=True):
    """Display an SIP Debye block model as image."""
    if z is None:
        z = np.cumsum(np.hstack((0., thk)))

    plt.cla()
    pal = plt.cm.get_cmap()
    pal.set_under('w')
    pal.set_bad('w')
    if isinstance(M, pg.Vector):
        M = np.asarray(M)

    if islog:
        M = np.log10(M)

    M = M.reshape((len(z), len(tau)))
    im = plt.imshow(M, interpolation='nearest', cmap=pal)
    if cmax is None:
        cmax = np.max(M)

    if cmax is None:
        cmax = np.max(M)

    im.set_clim((cmin, cmax))

    a = []
    for t in tau[::2]:
        a.append(f"{rndig(t * 1000, 2):g}")

    plt.xticks(np.arange(0, len(tau), 2), a)

    a = []
    for zi in z:
        a.append(str(zi))

    plt.yticks(np.arange(len(z)) - 0.5, a)
    plt.xlabel(r'$\tau$ in ms')
    plt.ylabel('z in m')
    plt.ylim((len(z) - 0.5, -0.5))
    plt.colorbar(orientation='horizontal', aspect=40, shrink=0.6)

    if res is not None:
        xl = plt.xlim()[1]
        for i in range(len(res)):
            plt.text(xl, i, rf' {rndig(res[i], 2):g} $\Omega$m')

    lgm = np.zeros((len(z), 1))
    tch = np.zeros((len(z), 1))
    lgt = np.log(tau)
    if islog:
        M = 10**M

    for n in range(len(M)):
        m = np.abs(M[n])
        tch[n] = np.sum(m)
        lgm[n] = np.exp(np.sum(m * lgt) / np.sum(m))

    tpos = np.interp(np.log(lgm), np.log(tau), np.arange(len(tau)))
    plt.plot(tpos, np.arange(len(z)), 'w*')

    plt.title('logarithmized spectral chargeability')
    plt.show()
    return lgm, tch


class DebyeModelling(pg.core.ModellingBase):
    """Forward operator for Debye decomposition."""

    def __init__(self, fvec, tvec=None, zero=False, verbose=False):

        if tvec is None:
            tvec = np.logspace(-4, 0, 5)

        mesh = pg.meshtools.createMesh1D(len(tvec))

        if zero:
            mesh.cell(0).setMarker(-1)
            mesh.cell(len(tvec) - 1).setMarker(1)

        pg.core.ModellingBase.__init__(self, mesh, verbose)
        self.f_ = pg.asvector(fvec)
        self.t_ = tvec
        self.zero_ = zero

    def response(self, par):
        """Phase spectrum as function of spectral chargeabilities."""
        y = pg.Vector(len(self.f_), 0.0)
        for (t, p) in zip(self.t_, par, strict=False):
            wt = self.f_ * 2.0 * np.pi * t
            y = y + wt / (wt * wt + 1.) * p

        return y


def DebyeDecomposition(fr, phi, maxfr=None, tv=None, verbose=False,
                       zero=False, err=0.25e-3, lam=10., blocky=False):
    """Debye decomposition of a phase spectrum."""
    if maxfr is not None:
        idx = (fr <= maxfr) & (phi >= 0.)
        phi1 = phi[idx]
        fr1 = fr[idx]
        print("using frequencies from ", np.min(fr), " to ", np.max(fr), "Hz")
    else:
        phi1 = phi
        fr1 = fr

    if tv is None:
        tmax = 1. / np.min(fr1) / 2. / np.pi * 4.
        tmin = 1. / np.max(fr1) / 2. / np.pi / 8.
        tvec = np.logspace(np.log10(tmin), np.log10(tmax), 30)
    else:
        tvec = tv

    f = DebyeModelling(fr1, tvec, zero=zero)
    tvec = f.t_
    tm = pg.trans.TransLog()
    start = pg.Vector(len(tvec), 1e-4)
    if zero:
        f.region(-1).setConstraintType(0)  # smoothness
        f.region(0).setConstraintType(1)  # smoothness
        f.region(1).setConstraintType(0)  # min length
        f.regionManager().setInterRegionConstraint(-1, 0, 1.)
        f.regionManager().setInterRegionConstraint(0, 1, 1.)
        f.region(-1).setTransModel(tm)
        f.region(0).setTransModel(tm)
        f.region(1).setTransModel(tm)
        f.region(-1).setModelControl(1000.)
        f.region(1).setModelControl(1000.)
    else:
        f.regionManager().setConstraintType(1)  # smoothness

    inv = pg.Inversion(pg.asvector(phi1 * 1e-3), f, verbose)
    inv.setAbsoluteError(pg.Vector(len(fr1), err))
    inv.setLambda(lam)
    inv.setModel(start)
    inv.setBlockyModel(blocky)
    if zero:
        inv.setReferenceModel(start)
    else:
        inv.setTransModel(tm)

    mvec = inv.run()
    resp = inv.response()

    return tvec, mvec, np.array(resp) * 1e3, idx


class DoubleColeColeModelling(pg.core.ModellingBase):
    """Modelling using two Cole-Cole terms."""

    def __init__(self, mesh, fvec, si=1.0, verbose=False):
        pg.core.ModellingBase.__init__(self, mesh, verbose)
        self.f_ = fvec
        self.si_ = si

    def response(self, par):
        """Phase response response of double Cole Cole model."""
        y = pg.Vector(self.f_.size(), 0.0)
        wti = self.f_ * par[1] * 2.0 * np.pi
        wte = self.f_ * par[4] * 2.0 * np.pi
        for i in range(0, y.size()):
            cpI = 1. / (np.power(wti[i] * 1j, par[2]) + 1.)
            cpE = 1. / (np.power(wte[i] * 1j, par[5]) + 1.)
            y[i] = - np.imag(cpI) * par[0] - np.imag(cpE) * par[3] * self.si_
#            y[i] = - par[0] - np.imag(cpE) * par[3] * self.si_

        return y

# can very probably be deleted
def ReadAndRemoveEM(filename, readsecond:bool=False, doplot:bool=False,
                    dellast:bool=True, ePhi:float=0.5, ePerc:float=1,
                    lam:float=2000):
    """Read res1file and remove EM effects using a double-Cole-Cole model."""
    fr, rhoa, phi, _drhoa, dphi = load(filename),
                                            #   readsecond,
                                            #   dellast=dellast)
    # forward problem
    mesh = pg.meshtools.createMesh1D(1, 6)  # 6 independent parameters
    f = DoubleColeColeModelling(mesh, pg.asvector(fr), phi[2] / abs(phi[2]))
    f.regionManager().loadMap("region.control")
    model = f.createStartVector()

    # inversion
    inv = pg.Inversion(phi, f, True, False)
    inv.setAbsoluteError(phi * ePerc * 0.01 + ePhi / 1000.)
    inv.setRobustData(True)

    # inv.setCWeight(pg.Vector(6, 1.0)) # wozu war das denn gut?
    inv.setMarquardtScheme(0.8)
    inv.setLambda(lam)
    inv.setModel(model)
    erg = inv.run()
    inv.echoStatus()
    chi2 = inv.chi2()
    mod0 = pg.Vector(erg)
    mod0[0] = 0.0  # set IP term to zero to obtain pure EM term
    emphi = f.response(mod0)
    resid = (phi - emphi) * 1000.

    if doplot:
        s = "IP: m= " + str(rndig(erg[0])) + " t=" + str(rndig(erg[1])) + \
            " c =" + str(rndig(erg[2]))
        s += "  EM: m= " + str(rndig(erg[3])) + " t=" + str(rndig(erg[4])) + \
            " c =" + str(rndig(erg[5]))
        fig = plt.figure(1)
        fig.clf()
        ax = plt.subplot(111)
        plt.errorbar(
            fr,
            phi *
            1000.,
            yerr=dphi *
            1000.,
            fmt='x-',
            label='measured')
        ax.set_xscale('log')
        plt.semilogx(fr, emphi * 1000., label='EM term (CC)')
        plt.errorbar(fr, resid, yerr=dphi * 1000., label='IP term')
        ax.set_yscale('log')
        plt.xlim((min(fr), max(fr)))
        plt.ylim((0.1, max(phi) * 1000.))
        plt.xlabel('f in Hz')
        plt.ylabel(r'-$\phi$ in mrad')
        plt.grid(True)
        plt.title(s)
        plt.legend(loc=2)  # ('measured','2-cole-cole','residual'))
        fig.show()

    return np.array(fr), np.array(rhoa), np.array(resid), np.array(
        phi) * 1e3, dphi, chi2, np.array(emphi) * 1e3
