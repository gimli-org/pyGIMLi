#!/usr/bin/env python
"""Frequency Domain Electromagnetics (FDEM) functions and class."""

import numpy as np

import pygimli as pg
from pygimli.viewer.mpl import show1dmodel, drawModel1D
from .hemmodelling import HEMmodelling

freqMaxMin10 = 2**np.arange(10) * 110.
freqMaxMin8 = 2**np.arange(8) * 110.
freqResolveHCP = np.array([387., 1820., 8330., 41500., 133400.])
freqResolveVCX = 5410.
freqResolveHCPOld = np.array([380., 1770., 8300., 41000., 129500.])

fBKS36a = np.array([386, 1817, 8360, 41420, 133200, 5390])
rBKS36a = np.array([7.938, 7.931, 7.925, 7.912, 7.918, 9.055])
fBKS60 = np.array([380, 1773, 8300, 41000, 129500, 5410])
rBKS60 = np.array([7.918, 7.918, 7.957, 8.033, 7.906, 9.042])


def cmapDAERO():
    """Standardized colormap from A-AERO projects (purple=0.3 to red=500)."""
    from matplotlib.colors import LinearSegmentedColormap
    CMY = np.array([
        [127, 255, 31], [111, 255, 47], [95, 255, 63], [79, 255, 79],
        [63, 255, 95], [47, 255, 111], [31, 255, 127], [16, 255, 159],
        [0, 223, 159], [0, 191, 159], [0, 159, 207], [0, 127, 175],
        [0, 95, 175], [0, 63, 175], [0, 47, 175], [0, 31, 191], [0, 0, 255],
        [0, 0, 159], [15, 0, 127], [47, 0, 143], [79, 0, 143], [111, 0, 143],
        [143, 0, 127], [159, 31, 63], [175, 47, 31], [207, 63, 0],
        [223, 111, 0], [231, 135, 0], [239, 159, 0], [255, 191, 47],
        [239, 199, 63], [223, 207, 79], [207, 239, 111]], dtype=float)
    RGB = 1.0 - CMY/255
    return LinearSegmentedColormap.from_list('D-AERO', RGB)


def xfplot(ax, DATA, x, freq, everyx=5, orientation='horizontal', aspect=30,
           label=None, cMap="Spectral_r"):
    """Plot a matrix according to x and frequencies."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    nt = list(range(0, len(x), everyx))
    im = ax.matshow(DATA.T, interpolation='nearest', cmap=cMap)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xticks(nt)
    ax.set_xticklabels(["%g" % xi for xi in x[nt]])
    ax.set_yticks(list(range(0, len(freq) + 1, 2)))
    ax.set_yticklabels(["%g" % freq[i] for i in range(0, len(freq), 2)])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('f [Hz]')
    ax.xaxis.set_label_position('top')
    ax.set_aspect("auto")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='5%', pad=0.3)
    pg.plt.colorbar(im, ax=ax, cax=cax, orientation=orientation, aspect=aspect)
    if label is not None:
        cax.set_title(label)
    # plt.colorbar(im, ax=ax, orientation=orientation, aspect=aspect)
    return im


class FDEM2dFOPold(pg.core.ModellingBase):
    """Old variant of 2D FOP (to be deleted)."""

    def __init__(self, data, nlay=2, verbose=False):
        """Initialize with data and number of layers."""
        pg.core.ModellingBase.__init__(self, verbose)
        self.nlay = nlay
        self.FOP1d = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        self.mesh_ = pg.meshtools.createMesh1D(self.nx, 2 * nlay - 1)
        self.setMesh(self.mesh_)

    def response(self, model):
        """Yield forward model response."""
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        resp = pg.Vector(0)
        for modi in modA:
            resp = pg.cat(resp, self.FOP1d.response(modi))

        return resp


class FDEM2dFOP(pg.core.ModellingBase):
    """FDEM 2d-LCI modelling class based on BlockMatrices."""

    def __init__(self, data, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers."""
        super(FDEM2dFOP, self).__init__(verbose)
        self.nlay = nlay
        self.header = {}
        self.pos, self.z, self.topo = None, None, None
        self.FOP = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        npar = 2 * nlay - 1
        self.mesh1d = pg.meshtools.createMesh1D(self.nx, npar)
        self.mesh_ = pg.meshtools.createMesh1D(self.nx, 2 * nlay - 1)
        self.setMesh(self.mesh_)

        # self.J = NDMatrix(self.nx, self.nf*2, npar)
        self.J = pg.matrix.BlockMatrix()
        self.FOP1d = []
        for i in range(self.nx):
            self.FOP1d.append(pg.core.FDEM1dModelling(
                nlay, data.freq(), data.coilSpacing, -data.height))
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, self.nf * 2 * i, npar * i)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())

    def response(self, model):
        """Cut together forward responses of all soundings."""
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        resp = pg.Vector(0)
        for modi in modA:
            resp = pg.cat(resp, self.FOP.response(modi))

        return resp

    def createJacobian(self, model):
        """Create Jacobian matrix by creating individual Jacobians."""
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class HEM1dWithElevation(pg.core.ModellingBase):
    """Airborne FDEM modelling including variable bird height."""

    def __init__(self, frequencies, coilspacing, nlay=2, verbose=False):
        """Set up class by frequencies and geometries."""
        pg.core.ModellingBase.__init__(self, verbose)
        self.nlay_ = nlay  # real layers (actually one more!)
        self.FOP_ = pg.core.FDEM1dModelling(nlay + 1, frequencies,
                                            coilspacing, self.height)
        self.mesh_ = pg.meshtools.createMesh1D(nlay, 2)  # thicknesses & res
        self.mesh_.cell(0).setMarker(2)
        self.setMesh(self.mesh_)

    def response(self, model):
        """Return forward response for a given model."""
        thk = model(0, self.nlay)  # all thicknesses including bird height
        res = model(self.nlay - 1, self.nlay * 2)
        res[0] = 10000.
        return self.FOP_.response(pg.cat(thk, res))


class FDEM():
    """Class for managing Frequency Domain EM data and their inversions."""

    def __init__(self, x=None, freqs=None,
                 coilSpacing=None, inphase=None, outphase=None,
                 filename=None, scaleFreeAir=False):
        r"""Initialize data class and load data. Provide filename or data.

        If filename is given, data is loaded, overwriting settings.

        Parameters
        ----------
        x: array
            Array of measurement positions

        freq: array
            Measured frequencies

        coilSpacing : float
            Distance between 2 two coils

        inphase : array
            real part of :math:`|amplitude| * \exp^{i phase}`

        outphase : array
            imaginary part of :math:`|amplitude| * \exp^{i phase}`

        filename : str
            Filename to read from. Supported: .xyz (MaxMin), .txt (Emsys)

        scaleFreeAir : bool
            Scale inphase and outphase data by free air (primary) solution
        """
        if isinstance(x, str) and freqs is None:  # string/file init
            filename = x

        self.x = x
        self.frequencies = freqs
        self.coilSpacing = coilSpacing
        self.scaling = "ppm"

        self.IP = inphase
        self.OP = outphase
        self.ERR = None

        self.height = 1.0  # standard height for MaxMin/Promys devices
        self.fop = None  # better apply MethodManger base interface
        self.transData, self.transRes, self.transThk = None, None, None

        if filename:
            # check if filename extension is TXT or CSV
            fl = filename.lower()
            if fl.endswith('.txt') or fl.endswith('.csv'):
                try:
                    self.importMaxMinData(filename)
                except Exception:
                    self.importEmsysAsciiData(filename)
            else:
                self.importMaxMinData(filename)

        if np.any(self.frequencies):
            self.isActiveFreq = self.frequencies > 0.0
            self.activeFreq = np.nonzero(self.isActiveFreq)[0]

        if scaleFreeAir:
            freeAirSolution = self.FOP().freeAirSolution()
            self.IP /= freeAirSolution
            self.OP /= freeAirSolution

    def __repr__(self):
        """Create string representation."""
        if self.x is None:
            part1 = "<FDEMdata: sounding with {:d} frequencies".format(
                len(self.frequencies))
        else:
            part1 = "<FDEMdata: {:d} soundings with {:d} frequencies".format(
                len(self.x), len(self.frequencies))
        if self.coilSpacing is not None:
            cs = self.coilSpacing
            if hasattr(cs, '__iter__'):
                if len(cs) > 1:
                    part2 = "coil spacing is {:f}-{:f} m>".format(min(cs),
                                                                  max(cs))
            else:
                part2 = "coil spacing is {:f} m".format(cs)
            return part1 + ' , ' + part2
        else:
            return part1

    def importEmsysAsciiData(self, filename):
        """Import data from emsys text export file.

        columns: no, pos(1-3), separation(4), frequency(6), error(8),
        inphase (9-11), outphase (12-14),
        reads: positions, data, frequencies, error and geometry
        """
        cols = (1, 4, 6, 8, 9, 12, 15, 16)
        xx, sep, f, pf, ip, op, hmod, q = np.loadtxt(filename, skiprows=1,
                                                     usecols=cols, unpack=True)
        err = q / pf * 100.  # percentage of primary field

        if len(np.unique(sep)) > 1:
            print("Warning! Several coil spacings present in file!")

        self.coilSpacing = np.median(sep)
        f = np.round_(f)
        self.frequencies, mf, nf = np.unique(f, True, True)
        x, mx, nx = np.unique(xx, True, True)
        self.IP = np.ones((len(x), len(f))) * np.nan
        self.OP = np.ones((len(x), len(f))) * np.nan
        self.ERR = np.ones((len(x), len(f))) * np.nan

        for i in range(len(f)):
            self.IP[nx[i], nf[i]] = ip[i]
            self.OP[nx[i], nf[i]] = op[i]
            self.ERR[nx[i], nf[i]] = err[i]

    def readHEMData(self, filename, takeevery=1, choosevcp=True):
        """Read RESOLVE type airborne EM data from .XYZ file."""
        self.header = {}
        keyword = ''
        i = 0
        with open(filename) as f:
            for i, line in enumerate(f):
                if line[0] == '/':
                    line = line[1:].strip('\n').replace(',', '').replace('AND',
                                                                         '')
                    try:
                        result = [float(co) for co in line.split()]
                    except ValueError:
                        result = line.split()
                    if len(result) == 1:
                        result = result[0]
                    if keyword:
                        if isinstance(keyword, list):
                            for kw, res in zip(keyword, result):
                                self.header[kw] = res
                        else:
                            self.header[keyword] = result
                        keyword = ''
                    else:
                        keyword = result
                else:
                    break
            line = f.readline()
            print(line)
#            tmp = np.genfromtxt(fname=f, autostrip=True, comments='/',
#                skip_header=0, dtype=float, names=1, case_sensitive='lower',
#                missing_values='*', filling_values=-9999, skip_footer=1)
        tmp = np.genfromtxt(
            fname=filename, autostrip=True, comments='/',
            skip_header=i+1, dtype=float, names=True, case_sensitive='lower',
            missing_values='*', filling_values=-9999, skip_footer=1)
        # read properties from header
        if choosevcp:
            ivcp = np.nonzero(np.array(self.header['COILGEOMETRY']) == 1)[0]
        else:
            ivcp = range(len(self.header['FREQUENCY']))
        self.frequencies = np.array(self.header['FREQUENCY'])[ivcp]
        self.coilSpacing = np.array(self.header['COILSEPERATION'])[ivcp]

        # read properties from data block
        names = tmp.dtype.names
        if 'lon' in names and 'lat' in names:
            utm = pg.utils.getUTMProjection(zone=32)
            x, y = utm(tmp['lon'], tmp['lat'])
        else:
            x, y = tmp['x'], tmp['y']

        self.pos = np.column_stack((x, y))[::takeevery]
        dx = np.sqrt(np.diff(self.pos[:, 0])**2 + np.diff(self.pos[:, 1])**2)
        self.x = np.hstack((0., np.cumsum(dx)))
        self.z = tmp['h_laser'][::takeevery]
        self.topo = tmp['topo'][::takeevery]
        IP = np.column_stack([tmp['real_'+str(i+1)] for i in ivcp])
        OP = np.column_stack([tmp['quad_'+str(i+1)] for i in ivcp])
        # better do a decimation or running average here
        self.IP = IP[::takeevery, :]
        self.OP = OP[::takeevery, :]
        self.isActiveFreq = self.frequencies > 0.0
        self.activeFreq = np.nonzero(self.isActiveFreq)[0]

    def importIPXData(self, filename, verbose=False):
        """Import MaxMin IPX format with pos, data, frequencies & geometry."""
        delim = None
        fid = open(filename)
        aline = ''
        i = 0  # just in case there is no header
        for i, aline in enumerate(fid):
            if aline.split()[0][0].isdigit():  # number found
                break
            elif aline.find('COIL') > 0:  # [:6] == '/ COIL':
                self.coilSpacing = float(aline.replace(':', ': ').split()[-2])
            elif aline.find('FREQ') > 0:  # [:6] == '/ FREQ':
                mya = aline[aline.find(':') + 1:].replace(',', ' ').split()
                myf = [float(aa) for aa in mya if aa[0].isdigit()]
                self.frequencies = np.array(myf)

        fid.close()

        if verbose:
            print("CS=", self.coilSpacing, "F=", self.frequencies)
        if aline.find(',') > 0:
            delim = ','

        nf = len(self.frequencies)
        if verbose:
            print("delim=", delim, "nf=", nf)
        A = np.loadtxt(filename, skiprows=i, delimiter=delim, comments='/').T
        x, y, self.IP, self.OP = A[0], A[1], A[
            2:nf * 2 + 2:2].T, A[3:nf * 2 + 2:2].T
        if max(x) == min(x):
            self.x = y
        else:
            self.x = x

    def importMaxMinData(self, filename, verbose=False):
        """Import MaxMin ASCII export (.txt) data."""
        with open(filename) as fid:
            lines = fid.readlines()

        self.coilSpacing = 99.9
        f, re, im, err, cond = [], [], [], [], []
        x, RE, IM, ERR, COND = [], [], [], [], []
        for i, line in enumerate(lines):
            stline = line.split()
            if line.startswith("Coil Sep"):
                self.coilSpacing = float(stline[-1])
            if len(stline) > 3 and stline[3].find("Stn") >= 0:
                x.append(float(stline[4]))
                if len(re) > 0:
                    RE.append(np.array(re))
                    IM.append(np.array(im))
                    ERR.append(np.array(err))
                    COND.append(np.array(cond))
                    f, re, im, err, cond = [], [], [], [], []

            if len(stline) > 0 and stline[0] == "MAX1":  # data
                f.append(float(stline[1]))
                re.append(float(stline[3]))
                im.append(float(stline[5]))
                err.append(float(stline[7]))
                cond.append(float(stline[9]))

        if len(re) > 0:
            RE.append(np.array(re))
            IM.append(np.array(im))
            ERR.append(np.array(err))
            COND.append(np.array(cond))

        self.x = np.array(x)
        self.frequencies = np.array(f)
        self.IP = np.array(RE)
        self.OP = np.array(IM)
        self.ERR = np.array(ERR)

    def deactivate(self, fr):
        """Deactivate a single frequency."""
        fi = np.nonzero(np.absolute(self.frequencies / fr - 1.) < 0.1)
        self.isActiveFreq[fi] = False
        self.activeFreq = np.nonzero(self.isActiveFreq)[0]

    def freq(self):
        """Return active (i.e., non-deactivated) frequencies."""
        return self.frequencies[self.activeFreq]

    def FOP(self, nlay=2, useHEM=1):  # createFOP deciding upon block or smooth
        """Forward modelling operator using a block discretization.

        Parameters
        ----------
        nlay : int
            Number of blocks
        """
        if useHEM:
            return HEMmodelling(nlay, self.height, f=self.freq(),
                                r=self.coilSpacing, scaling=self.scaling)
        else:
            return pg.core.FDEM1dModelling(nlay, self.freq(), self.coilSpacing,
                                           -self.height)

    def FOPsmooth(self, zvec):
        """Forward modelling operator using fixed layers (smooth inversion).

        Parameters
        ----------
        zvec : array
        """
        return pg.FDEM1dRhoModelling(zvec, self.freq(), self.coilSpacing,
                                     -self.height)

    def selectData(self, xpos=0):
        """Select sounding at a specific position or by number.

        Retrieve inphase, outphase and error(if exist) vector from index
        or near given position

        Parameters
        ----------
        xpos : int | float
            index (int) or position (float) along profile to choose

        Returns
        -------
            IP : array
            OP : array
            ERR : array or None (if no error is specified)
        """
        # check for index
        if isinstance(xpos, int) and (xpos < len(self.x)) and (xpos >= 0):
            n = xpos
        else:
            n = np.argmin(np.absolute(self.x - xpos))

        self.height = self.z[n]
        ip = self.IP[n, self.activeFreq]
        op = self.OP[n, self.activeFreq]
        err = None
        if self.ERR is not None:
            err = self.ERR[n, self.activeFreq]

        return ip, op, err

    def error(self, xpos=0):
        """Return error as vector."""
        _, _, err = self.selectData(xpos)
        return err

    def datavec(self, xpos=0):
        """Extract data vector (stack in and out phase) for given pos/no."""
        ip, op, _ = self.selectData(xpos)
        return np.hstack((ip, op))

    def errorvec(self, xpos=0, minvalue=0.0):
        """Extract error vector for a give position or sounding number."""
        _, _, err = self.selectData(xpos)
        return np.tile(np.maximum(err * 0.7071, minvalue), 2)
#        return pg.asvector(np.tile(np.maximum(err * 0.7071, minvalue), 2))

    def invBlock(self, xpos=0, nlay=2, noise=1.0, show=True, stmod=30.,
                 lam=1000., lBound=0., uBound=0., verbose=False, **kwargs):
        """Create and return Gimli inversion instance for block inversion.

        Parameters
        ----------
        xpos : array
            position vector

        nLay : int
            Number of layers of the model to be determined OR
            vector of layer numbers OR forward operator

        noise : float
            Absolute data err in percent

        stmod : float or pg.Vector
            Starting model

        lam : float
            Global regularization parameter lambda.

        lBound : float
            Lower boundary for the model

        uBound : float
            Upper boundary for the model. 0 means no upper booundary

        verbose : bool
            Be verbose
        """
        # self.transThk = pg.trans.TransLog()
        # self.transRes = pg.trans.TransLogLU(lBound, uBound)
        # self.transData = pg.trans.Trans()
        self.transData = pg.trans.TransSymLog(tol=0.1)
        self.transLog = pg.trans.TransLog()

        useHEM = kwargs.pop("useHEM", True)
        # EM forward operator
        if isinstance(nlay, pg.core.FDEM1dModelling):
            self.fop = nlay
        else:
            self.fop = self.FOP(nlay, useHEM=useHEM)

        dataVec = self.datavec(xpos)

        # self.fop.region(0).setTransModel(self.transThk)
        # self.fop.region(1).setTransModel(self.transRes)

        if isinstance(noise, float):
            errorVec = pg.Vector(len(dataVec), noise)
        else:
            errorVec = pg.asvector(noise)

        # independent EM inversion

        if isinstance(stmod, float):  # real model given
            model = pg.Vector(nlay * 2 - 1, stmod)
            model[0] = 2.
        else:
            if len(stmod) == nlay * 2 - 1:
                model = stmod
            else:
                model = pg.Vector(nlay * 2 - 1, 30.)

            print("Model", model)
        if 1:
            from pygimli.frameworks import MarquardtInversion
            self.inv = MarquardtInversion(fop=self.fop, verbose=verbose,
                                          debug=True)
            self.inv.dataTrans = self.transData
            self.inv.modelTrans = self.transLog
            # self.dataTrans = self.transData
            self.model1d = self.inv.run(dataVec, np.abs(errorVec/dataVec),
                                        lam=lam, startModel=model, **kwargs)
            response = self.inv.response
        else:
            self.inv = pg.core.RInversion(dataVec, self.fop, self.transData,
                                          verbose)
            self.inv.setAbsoluteError(errorVec)
            self.inv.setLambda(lam)
            self.inv.setMarquardtScheme(0.8)
            self.inv.setDeltaPhiAbortPercent(0.5)
            self.inv.setModel(model)
            self.model1d = self.inv.run()
            response = self.inv.response()

        if show:
            self.plotData(xpos=xpos, response=response)

        return self.model1d

    def plotData(self, xpos=0, response=None, error=None, ax=None,
                 marker='bo-', rmarker='rx-', clf=True, addlabel='', nv=2):
        """Plot data as curves at given position."""
        ip, op, err = self.selectData(xpos)

        if error is not None and err is not None:
            error = err

        fr = self.freq()

        if ax is None:
            _, ax = pg.plt.subplots(nrows=1, ncols=nv)
            ipax = ax[-2]
        else:
            ipax = ax[0]

        markersize = 4

        if error is not None:
            markersize = 2

        ipax.semilogy(ip, fr, marker, label='obs' + addlabel,
                      markersize=markersize)

        if error is not None and len(error) == len(ip):
            ipax.errorbar(ip, fr, xerr=error)

        # ipax.set_axis('tight')

        if error is not None:
            ipax.ylim((min(fr) * .98, max(fr) * 1.02))

        ipax.grid(True)
        ipax.set_xlabel('inphase [ppm]')
        ipax.set_ylabel('f [Hz]')

        if response is not None:
            rip = np.asarray(response)[:len(ip)]
            ipax.semilogy(rip, fr, rmarker, label='syn' + addlabel)

        ipax.legend(loc='best')

        opax = None

        if ax is None:
            opax = pg.plt.subplot(1, nv, nv)
        else:
            opax = ax[-1]

        opax.semilogy(op, fr, marker, label='obs' + addlabel,
                      markersize=markersize)

        if error is not None and len(error) == len(ip):
            opax.errorbar(op, fr, xerr=error)

        if response is not None:
            rop = np.asarray(response)[len(ip):]
            opax.semilogy(rop, fr, rmarker, label='syn' + addlabel)

#        opax.set_axis('tight')
        if error is not None:
            opax.ylim((min(fr) * .98, max(fr) * 1.02))

        opax.grid(True)
        opax.set_xlabel('outphase [ppm]')
        opax.set_ylabel('f [Hz]')
        opax.legend(loc='best')
        # plt.subplot(1, nv, 1)
        return ax

    # def plotDataOld(self, xpos=0, response=None,
    #                 marker='bo-', rmarker='rx-', clf=True):
    #     """Plot data as curves at given position."""
    #     ip, op = self.selectData(xpos)
    #     fr = self.freq()

    #     if clf:
    #         plt.clf()

    #     plt.subplot(121)
    #     plt.semilogy(ip, fr, marker, label='obs')
    #     plt.axis('tight')
    #     plt.grid(True)
    #     plt.xlabel('inphase [%]')
    #     plt.ylabel('f [Hz]')

    #     if response is not None:
    #         rip = np.asarray(response)[:len(ip)]
    #         plt.semilogy(rip, fr, rmarker, label='syn')

    #     plt.legend(loc='best')

    #     plt.subplot(122)
    #     plt.semilogy(op, fr, marker, label='obs')

    #     if response is not None:
    #         rop = np.asarray(response)[len(ip):]
    #         plt.semilogy(rop, fr, rmarker, label='syn')

    #     plt.axis('tight')
    #     plt.grid(True)
    #     plt.xlabel('outphase [%]')
    #     plt.ylabel('f [Hz]')
    #     plt.legend(loc='best')
    #     plt.show()

    #     return

    def showModelAndData(self, model, xpos=0, response=None, figsize=(8, 6)):
        """Show both model and data with response in subfigures."""
        fig, ax = pg.plt.subplots(1, 3, figsize=figsize)

        model = np.asarray(model)
        nlay = int((len(model) + 1) / 2)

        thk = model[:nlay - 1]
        res = model[nlay - 1: 2 * nlay - 1]

        drawModel1D(ax[0], thk, res, plotfunction='semilogx')

        self.plotData(xpos, response, ax=ax[1:3], clf=False)
        return fig, ax

    def plotAllData(self, orientation='horizontal', aspect=1000,
                    outname=None, show=False, figsize=(11, 8), everyx=None):
        """Plot data along a profile as image plots for IP and OP."""
        if self.x is None:
            raise Exception("No measurement position array x given")

        freq = self.freq()
        nr = 2

        if self.ERR is not None:
            nr = 3

        if everyx is None:
            everyx = len(self.x) // 10

        _, ax = pg.plt.subplots(ncols=1, nrows=nr, figsize=figsize)
        xfplot(ax[0], self.IP[:, self.activeFreq], self.x, freq,
               orientation=orientation, aspect=aspect, everyx=everyx,
               label='inphase percent')
        # ax[0].set_title('inphase percent')

        xfplot(ax[1], self.OP[:, self.activeFreq], self.x, freq,
               orientation=orientation, aspect=aspect, everyx=everyx,
               label='outphase percent')
        # ax[1].set_title('outphase percent')

        if self.ERR is not None:
            xfplot(ax[2], self.ERR[:, self.activeFreq], self.x, freq,
                   orientation=orientation)
            ax[2].set_title('error percent')

        if outname is not None:
            pg.plt.savefig(outname)

        if show:
            pg.plt.show()

        return ax

    def plotModelAndData(self, model, xpos, response,
                         modelL=None, modelU=None):
        """Plot both model and data in subfigures."""
        self.plotData(xpos, response, nv=3)
        show1dmodel(model, color='blue')
        if modelL is not None and modelU is not None:
            pass
            # draw1dmodelErr(model, modelL, modelU)  # !!!!

        return

    def FOP2d(self, nlay):
        """2d forward modelling operator."""
        return FDEM2dFOP(self, nlay)

    def inv2D(self, nlay, lam=100., resL=1., resU=1000., thkL=1.,
              thkU=100., minErr=1.0):
        """2d LCI inversion class."""
        if isinstance(nlay, int):
            modVec = pg.Vector(nlay * 2 - 1, 30.)
            cType = 0  # no reference model
        else:
            modVec = nlay
            cType = 10  # use this as referencemodel
            nlay = (len(modVec) + 1) / 2

        # init forward operator
        self.f2d = self.FOP2d(nlay)

        # transformations
        self.transData = pg.trans.Trans()
        self.transThk = pg.trans.TransLogLU(thkL, thkU)
        self.transRes = pg.trans.TransLogLU(resL, resU)

        for i in range(nlay - 1):
            self.f2d.region(i).setTransModel(self.transThk)

        for i in range(nlay - 1, nlay * 2 - 1):
            self.f2d.region(i).setTransModel(self.transRes)

        # set constraints
        self.f2d.region(0).setConstraintType(cType)
        self.f2d.region(1).setConstraintType(cType)

        # collect data vector
        datvec = pg.Vector(0)

        for i in range(len(self.x)):
            datvec = pg.cat(datvec, self.datavec(i))

        # collect error vector
        if self.ERR is None:
            error = 1.0
        else:
            error = []
            for i in range(len(self.x)):
                err = np.maximum(self.ERR[i][self.activeFreq] * 0.701, minErr)
                error.extend(err)

        # generate starting model by repetition
        model = np.repeat(modVec, len(self.x))
        INV = pg.core.Inversion(datvec, self.f2d, self.transData)
        INV.setAbsoluteError(error)
        INV.setLambda(lam)
        INV.setModel(model)
        INV.setReferenceModel(model)

        return INV


if __name__ == "__main__":
    import argparse  # better get an argparser from method manager

    parser = argparse.ArgumentParser(description="usage: %prog [options] fdem")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="be verbose", default=False)
    parser.add_argument("-n", "--nLayers", dest="nlay",
                        help="number of layers", type=int, default="4")
    parser.add_argument("-x", "--xPos", dest="xpos", help="position/no to use",
                        type=float, default="0")
    parser.add_argument("-l", "--lambda", dest="lam",
                        help="init regularization", type=float, default="30")
    parser.add_argument("-e", "--error", dest="err", help="error estimate",
                        type=float, default="1")
    parser.add_argument("datafile")
    options = parser.parse_args()

    if options.verbose:
        __verbose__ = True

    print(options)
    fdem = FDEM(options.datafile)
    print(fdem)
    datafile = options.datafile

    numlay = options.nlay
    xvector = options.xpos
    name = datafile.lower().rstrip('.xyz')
    fdem = FDEM(datafile)
    print(fdem)
    fdem.deactivate(56320.)  # do not use highest frequency
    fdem.plotAllData(outname=name + '-alldata.pdf')
    # that's bullshit (that's what we have the class for)
    inv = fdem.invBlock(xpos=xvector, lam=options.lam, nlay=options.nlay,
                        noise=options.err, verbose=False)
    mymodel = np.asarray(inv.run())
    inv.echoStatus()
    print("thk = ", mymodel[:numlay - 1])
    print("res = ", mymodel[numlay - 1:])
    figure, axes = fdem.showModelAndData(mymodel, xvector, inv.response())
    INV = fdem.inv2D(options.nlay)
    INV.run()
#    fig.savefig(name+str(xpos)+'-result.pdf', bbox_inches='tight')
    pg.plt.show()
