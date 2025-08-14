"""ERT manager (derived) with FD or TD IP inversion."""
import numpy as np
import pygimli as pg
from .ertManager import ERTManager
from .ertModelling import ERTModelling
from .ipModelling import DCIPMModelling


class ERTIPManager(ERTManager):
    """Method manager for ERT including induced polarization (IP).

    This class should be use for any single IP data, which can
    be a single-frequency frequency-domain (FD) amplitude and
    phase, or a time-domain (TD) IP chargeability (one gate or an
    integral value).
    """

    def __init__(self, *args, **kwargs):
        """Initialize DC part of it (parent class).

        Parameters
        ----------
        fd : bool
            Frequency-domain, otherwise time-domain
        """
        self.isfd = kwargs.pop("fd", False)
        super().__init__(*args, **kwargs)

    def invertTDIP(self, ipdata='ip', **kwargs):
        """IP inversion in time domain."""
        if isinstance(ipdata, str):
            ipdata = self.data[ipdata]
        if max(ipdata) > 1:  # mV/V
            ipdata /= 1000
        mesh0 = pg.Mesh(self.paraDomain)
        mesh0.setCellMarkers(mesh0.cellCount())
        fopIP = DCIPMModelling(self.fop, mesh0, self.model,
                               response=self.inv.response)
        fopIP.createRefinedForwardMesh(True)
        self.invIP = pg.Inversion(fop=fopIP, verbose=True)
        self.invIP.modelTrans = pg.trans.TransLogLU(0.0, 1.0)
        relErr = kwargs.pop("relativeError", 0.03)
        absErr = kwargs.pop("absoluteError", 0.001)
        errorIP = pg.Vector(self.data.size(), relErr) + absErr / pg.abs(ipdata)
        kwargs.setdefault("lam", 100)
        kwargs.setdefault("startModel", pg.median(ipdata))
        kwargs.setdefault("verbose", True)
        self.modelIP = self.invIP.run(ipdata, errorIP, **kwargs)

    def invertFDIP(self, ipdata="ip", **kwargs):
        """IP inversion in frequency domain."""
        data = kwargs.pop("data", self.data)
        if isinstance(ipdata, str):
            ipdata = self.data[ipdata]

        if "iperr" in kwargs:
            iperr = kwargs["iperr"]
        elif data.haveData("iperr"):
            iperr = data["iperr"]
        else:
            iperr = 1

        if max(ipdata) > 3.15:  # mrad
            ipdata /= 1000
        if np.max(iperr) >= 1:
            iperr /= 1000

        complexData = pg.utils.toComplex(data["rhoa"], ipdata)
        datavec = pg.cat(complexData.real, complexData.imag)
        errvec = pg.cat(data["err"], iperr / (pg.abs(ipdata)+1e-4))

        self.fopC = ERTModelling(sr=kwargs.pop("sr", False), verbose=True)
        self.fopC.setComplex(True)
        self.fopC.setData(self.data)
        self.fopC.setMesh(self.mesh)#, ignoreRegionManager=True)
        self.fopC.setDefaultBackground()
        ipmodel = pg.RVector(len(self.model), np.median(self.model) * 0.001)
        kwargs.setdefault("startModel", pg.cat(self.model, ipmodel))
        self.invIP = pg.Inversion(fop=self.fopC)  # , verbose=True, debug=True)
        self.invIP.dataTrans = kwargs.pop("dataTrans",
                                          "log" if min(datavec) > 0 else "lin")
        self.invIP.modelTrans = "log"
        kwargs.setdefault("verbose", True)
        kwargs.setdefault("isReference", True)  # "final phase improvement"
        self.modelC = pg.utils.toComplex(self.invIP.run(
            datavec, relativeError=errvec, **kwargs))
        self.modelIP = np.angle(self.modelC)  # mrad

    def showDCModel(self, **kwargs):
        """Explicitly show absolute of complex-valued inversion."""
        if self.isfd:
            return self.showResult(np.abs(self.modelC), **kwargs)
        else:
            return self.showModel(**kwargs)

    def showIPModel(self, **kwargs):
        """"Show IP model."""
        kwargs.setdefault("logScale", False)
        if self.isfd:
            kwargs.setdefault("label", r"$\phi$ (mrad)")
            kwargs.setdefault("cMap", "viridis")
        else:
            kwargs.setdefault("label", r"$m$ (mV/V)")
            kwargs.setdefault("cMap", "magma_r")

        return self.showModel(self.modelIP*1000, **kwargs)

    def showResults(self, reskw={}, ipkw={}, **kwargs):
        """Show DC and IP results.

        Parameters
        ----------
        reskw : dict
            keyword arguments for showing resistivity image
        ipkw : dict
            keyword arguments for showing IP image
        **kwargs : dict
            keyword arguments for showing resistivity image
        """
        _, ax = pg.plt.subplots(nrows=2, sharex=True)
        kwargs.setdefault("orientation", "vertical")
        kwargs.setdefault("xlabel", "x (m)")
        kwargs.setdefault("ylabel", "z (m)")
        reskw.setdefault("ax", ax[0])
        super().showResult(**reskw, **kwargs)
        ipkw.setdefault("ax", ax[1])
        ipkw.setdefault("logScale", False)
        ipkw.setdefault("cMin", 0)
        self.showIPModel(**ipkw, **kwargs)
        return ax

    def invert(self, *args, **kwargs):
        """Carry out DC and IP inversion."""
        super().invert(*args, **kwargs)  # DC first (not needed for FD)
        self.invertIP(**kwargs)

    def invertIP(self, **kwargs):
        """Invert IP data according to FD/TD settings."""
        if self.isfd:
            self.invertFDIP(**kwargs)
        else:
            self.invertTDIP(**kwargs)

    def invertDC(self, *args, **kwargs):
        # Needed if we want to do ERT first without the IP and do IP later
        super().invert(*args, **kwargs)

    def simulate(self, mesh, res, m, scheme=None, **kwargs):
        """."""
        from pygimli.physics.ert import ERTModelling
        data = scheme or pg.DataContainerERT(self.data)
        if hasattr(res[0], '__iter__'):  # ndim == 2
            if len(res[0]) == 2:  # res seems to be a res map
                resVec = pg.solver.parseArgToArray(res, mesh.cellCount(), mesh)
        elif len(res) == mesh.cellCount():
            resVec = res
        else:
            resVec = res[mesh.cellMarkers()]

        if hasattr(m[0], '__iter__'):  # ndim == 2
            if len(m[0]) == 2:  # res seems to be a res map
                mVec = pg.solver.parseArgToArray(m, mesh.cellCount(), mesh)
        elif len(m) == mesh.cellCount():
            mVec = m
        else:
            mVec = m[mesh.cellMarkers()]

        print(mesh, len(resVec), len(mVec))
        mesh0 = pg.Mesh(mesh)
        mesh0.setCellMarkers(mesh0.cellCount())
        fopDC = ERTModelling()
        fopDC.setData(scheme)
        fopDC.setMesh(mesh0)
        data["rhoa"] = fopDC.response(resVec)
        fopDC.createJacobian(resVec)
        fopIP = DCIPMModelling(fopDC, mesh, resVec)
        data["ma"] = fopIP.response(mVec)  # SI
        data["ip"] = data["ma"] * 1000  # mV/V
        return data

    def saveResult(self, folder=None, *args, **kwargs):
        """Save all results in given or date-based folder."""
        super().saveResult(folder=folder, **kwargs, ip=self.modelIP*1000)
