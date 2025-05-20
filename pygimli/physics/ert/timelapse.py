"""Timelapse ERT manager class."""
import os.path
from glob import glob
from datetime import datetime, timedelta
import numpy as np
import pygimli as pg
from pygimli.physics import ert
from .processing import combineMultipleData


# move general timelapse stuff to method-independent class
# class Timelapse():
#     mask
#     chooseTime
# class TimelapseERT(Timelapse)


class TimelapseERT():
    """Class for crosshole ERT data manipulation.

    Note that this class is to be split into a hierarchy of classes for general
    timelapse data management, timelapse ERT and crosshole ERT.
    You can load data, filter them data in the temporal or measuring axis, plot
    data, run inversion and export data and result files.

    """

    def __init__(self, filename=None, **kwargs):
        """Initialize class and possibly load data.

        Parameters
        ----------
        filename : str
            filename to load data, times, RHOA and ERR from
        data : DataContainerERT
            The data with quadrupoles for all
        times : np.array of datetime objects
            measuring times
        DATA : 2d np.array (data.size(), len(times))
            all apparent resistivities
        ERR : 2d np.array (data.size(), len(times))
            all apparent relative errors
        bhmap : array
            map electrode numbers to borehole numbers
        mesh : array
            mesh for inversion
        """
        self.data = kwargs.pop("data", None)
        self.DATA = kwargs.pop("DATA", [])
        self.ERR = kwargs.pop("ERR", [])
        self.times = kwargs.pop("times", [])
        self.mesh = kwargs.pop("mesh", None)
        self.pd = None
        self.name = "new"
        self.models = []
        self.responses = []
        self.chi2s = []
        self.model = None
        self.mgr = ert.ERTManager()
        if self.mesh is not None:
            self.mgr.setMesh(self.mesh)

        if filename is not None:
            if isinstance(filename, str):
                results = kwargs.pop("results", None)
                self.load(filename, **kwargs)                
                if results is not None:
                    self.loadResults(results)
            else:
                self.DATA = filename

        if np.any(self.DATA):
            if isinstance(self.DATA[0], pg.DataContainerERT):
                self.data, self.DATA, self.ERR = combineMultipleData(self.DATA)
            self.mask()

        if "name" in kwargs:
            self.name = kwargs["name"]

        nt = 0
        if np.any(self.DATA):
            nt = self.DATA.shape[1]

        if len(self.times) != nt:  # default: days from now
            self.times = datetime.now() + np.arange(nt) * timedelta(days=1)

    def __repr__(self):  # for print function
        """Return string representation of the class."""
        out = ['Timelapse ERT data:', self.data.__str__()]
        if np.any(self.DATA):
            out.append("{} time steps".format(self.DATA.shape[1]))
            if np.any(self.times):
                out[-1] += " from " + self.times[0].isoformat(" ", "minutes")
                out[-1] += " to " + self.times[-1].isoformat(" ", "minutes")

        return "\n".join(out)

    def load(self, filename, **kwargs):
        """Load or import data (or data files using *)."""
        if os.path.isfile(filename):
            self.data = ert.load(filename)
            if os.path.isfile(filename[:-4]+".rhoa"):
                self.DATA = np.loadtxt(filename[:-4]+".rhoa")
            if os.path.isfile(filename[:-4]+".err"):
                self.ERR = np.loadtxt(filename[:-4]+".err")
            if os.path.isfile(filename[:-4]+".times"):
                timestr = np.loadtxt(filename[:-4]+".times", dtype=str)
                self.times = np.array(
                    [datetime.fromisoformat(s) for s in timestr])
        elif "*" in filename:
            DATA = [ert.load(fname) for fname in glob(filename)]
            self.data, self.DATA, self.ERR = combineMultipleData(DATA)

        self.name = filename[:-4].replace("*", "All")

    def saveData(self, filename=None, masknan=True):
        """Save all data as datacontainer, times, rhoa and error arrays."""
        filename = filename or self.name
        if filename.endswith(".shm"):
            filename = filename[:-4]

        self.data.save(filename+".shm", "a b m n k")
        DATA = self.DATA.copy()
        if masknan:
            DATA[DATA.mask] = np.nan

        np.savetxt(filename+".rhoa", DATA, fmt="%6.2f")
        if np.any(self.ERR):
            np.savetxt(filename+".err", self.ERR, fmt="%6.2f")

        with open(filename+".times", "w", encoding="utf-8") as fid:
            for d in self.times:
                fid.write(d.isoformat()+"\n")

        self.name = filename

    def timeIndex(self, t):  #
        """Return index of closest timestep in times to t.

        Parameters
        ----------
        t : int|str|datetime
            datetime object or string
        """
        if isinstance(t, str):  # convert into datetime
            t = datetime.fromisoformat(t) # check others
        if isinstance(t, datetime): # detect closest point
            return np.argmin(np.abs(self.times-t))
        elif isinstance(t, (int, np.int32, np.int64)):
            return t
        elif hasattr(t, "__iter__"):
            return np.array([self.timeIndex(ti) for ti in t], dtype=int)
        else:
            raise TypeError("Unknown type", type(t))

    def filter(self, tmin=0, tmax=None, t=None, select=None, kmax=None):
        """Filter data set temporally or data-wise.

        Parameters
        ----------
        tmin, tmax : int|str|datetime
            minimum and maximum times to keep
        t : int|str|datetime
            time to remove
        select : list[int]
            times to select
        kmax : float
            maximum geometric factor to allow
        """
        if np.any(self.DATA):
            if select is not None:
                ind = self.timeIndex(select)
            else:
                tmin = self.timeIndex(tmin)  # converts dt/str to int
                if tmax is None:
                    tmax = self.DATA.shape[1]
                else:
                    tmax = self.timeIndex(tmax)

                ind = np.arange(tmin, tmax)
                if t is not None:
                    ind = np.setxor1d(ind, self.timeIndex(t))

            self.DATA = self.DATA[:, ind]
            if np.any(self.ERR):
                self.ERR = self.ERR[:, ind]

            if np.any(self.times):
                self.times = self.times[ind]

        if kmax is not None:
            ind = np.nonzero(np.abs(self.data["k"]) < kmax)[0]
            self.data["valid"] = 0
            self.data.markValid(ind)
            self.data.removeInvalid()
            if np.any(self.DATA):
                self.DATA = self.DATA[ind, :]

            if np.any(self.ERR):
                self.ERR = self.ERR[ind, :]

    def mask(self, rmin=0.1, rmax=1e6, emax=None):
        """Mask data (i.e. remove them from inversion).

        Parameters
        ----------
        rmin, rmax : float
            minimum and maximum apparent resistivity
        emax : float
            maximum error
        """
        self.DATA = np.ma.masked_invalid(self.DATA)
        self.DATA = np.ma.masked_outside(self.DATA, rmin, rmax)
        if emax is not None and np.any(self.ERR):
            self.DATA.mask = np.bitwise_or(self.DATA.mask, self.ERR > emax)

    def automask(self, dmax=0.3, nc=5, robustData=True):
        """Automatic outlier masking using dist to smoothed curve.
        
        Parameters
        ----------
        dmax : float
            maximum relative deviation
        nc : int
            number of coefficient pairs
        """
        from pygimli.frameworks import harmfit
        tt = np.array([ti.toordinal() for ti in self.times])
        for data in self.DATA:
            ddata = data[~data.mask].data
            if len(ddata) > 3:
                hf = harmfit(ddata, tt[~data.mask], verbose=False,
                            nc=nc, robustData=robustData, resample=tt)[0]
                misfit = np.abs(data/hf - 1)
                data.mask[misfit > dmax] = True

    def showData(self, v="rhoa", t=None, **kwargs):
        """Show data as pseudosections (single-hole) or cross-plot (crosshole)

        Parameters
        ----------
        v : str|array ["rhoa]
            array or field to plot
        t : int|str|datetime
            time to choose (can also be first argument)
        x, y : str|array ["a", "m"]
            values to use for x and y axes
        crossplot : bool [x and y given]
            force AB-MN crossplot
        kwargs : dict
            forwarded to ert.show or showDataContainerAsMatrix
        """
        if isinstance(v, (int, str)) and t is None:  # obviously t meant
            t = v
            v = "rhoa"

        kwargs.setdefault("cMap", "Spectral_r")
        if t is not None:
            t = self.timeIndex(t)
            rhoa = self.DATA[:, t].copy()
            v = rhoa.data
            v[rhoa.mask] = np.nan

        if kwargs.pop("crossplot", "x" in kwargs and "y" in kwargs):
            x = kwargs.pop("x", ["a", "b"])
            y = kwargs.pop("y", ["m", "n"])
            return pg.viewer.mpl.showDataContainerAsMatrix(
                self.data, x, y, v, **kwargs)
        else:
            if not "style" in kwargs:
                kwargs.setdefault("x", "a")
                kwargs.setdefault("y", "m")

            return self.data.show(v, **kwargs)

    def showTimeline(self, ax=None, **kwargs):
        """Show data timeline.

        Parameters
        ----------
        ax : mpl.Axes|None
            matplotlib axes to plot (otherwise new)
        a, b, m, n : int
            tokens to extract data from
        """
        if ax is None:
            _, ax = pg.plt.subplots(figsize=[8, 5])
        good = np.ones(self.data.size(), dtype=bool)
        lab = kwargs.pop("label", "ABMN") + ": "
        popit = []
        for k, v in kwargs.items():
            if self.data.haveData(k):
                good = np.bitwise_and(good, self.data[k] == v)
                popit.append(k)

        for k in popit:
            kwargs.pop(k)

        abmn = [self.data[tok] for tok in "abmn"]
        for i in np.nonzero(good)[0]:
            lab1 = lab + " ".join([str(tt[i]) for tt in abmn])
            ax.semilogy(self.times, self.DATA[i, :], "x-", label=lab1, **kwargs)

        ax.grid(True)
        ax.legend()
        ax.set_xlabel("time")
        ax.set_ylabel("resistivity (Ohmm)")
        return ax

    def fitReciprocalErrorModel(self, **kwargs):
        """Fit all data by analysing normal/reciprocal data.

        Parameters
        ----------
        show : bool
            show temporal behaviour of absolute & relative errors
        kwargs passed on to ert.fitReciprocalErrorModel)
            nBins : int
                number of bins to subdivide data (4 < data.size()//30 < 30)
            rel : bool [False]
                fit relative instead of absolute errors

        Returns
        -------
        p, a : array
            relative (p) and absolute (a) errors for every time step
        """
        data = self.data.copy()
        p = np.zeros(self.DATA.shape[1])
        a = np.zeros_like(p)
        show = kwargs.pop("show", False)  # avoid show single fits
        rel = kwargs.get("rel", True)
        for i, rhoa in enumerate(self.DATA.T):
            if isinstance(rhoa, np.ma.MaskedArray):
                rhoa = rhoa.data

            data['rhoa'] = rhoa
            try:
                pp, aa = ert.fitReciprocalErrorModel(data, **kwargs)
                if not rel:
                    pp, aa = aa, pp

                p[i], a[i] = pp, aa
            except:
                print("Could not get reciprocal model for timestep", i)

        if show:
            _, ax = pg.plt.subplots(nrows=2, sharex=True)
            ax[0].plot(self.times, p*100)
            ax[1].plot(self.times, a)
            ax[0].set_ylabel("relative error (%)")
            ax[1].set_ylabel("absolute error (Ohm)")
            ax[0].grid()
            ax[1].grid()

        return p, a

    def generateDataPDF(self, **kwargs):
        """Generate a pdf with all data as timesteps in individual pages.

        Iterates through times and calls showData into multi-page pdf
        """
        kwargs.setdefault("verbose", False)
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(self.name+'-data.pdf') as pdf:
            fig = pg.plt.figure(figsize=kwargs.pop("figsize", [5, 5]))
            for i in range(self.DATA.shape[1]):
                ax = fig.subplots()
                self.showData(t=i, ax=ax, **kwargs)
                ax.set_title(str(i)+": "+ self.times[i].isoformat(" ", "minutes"))
                fig.savefig(pdf, format='pdf')
                fig.clf()

    def generateTimelinePDF(self, key="a", filename=None, **kwargs):
        """Generate multipage PDF with timeline data.

        Parameters
        ----------
        key : str ['a']
            data key to sort measurements after
        filename : str [name+'time'+key+'.pdf']
            output pdf filename
        kwargs : dict
            passed to showTimeline
        """
        from matplotlib.backends.backend_pdf import PdfPages
        if filename is None:
            filename = self.name+'-time-'+key+'.pdf'
        with PdfPages(filename) as pdf:
            fig = pg.plt.figure()
            for a in np.unique(self.data[key]):
                ax = fig.subplots()
                self.showTimeline(ax=ax, a=a)
                fig.savefig(pdf, format='pdf')
                fig.clf()

    def chooseTime(self, t=None):
        """Return data for specific time.

        Parameters
        ----------
        t : int|str|datetime
        """
        if not isinstance(t, (int, np.int32, np.int64)):
            t = self.timeIndex(t)

        rhoa = self.DATA[:, t].copy()
        arhoa = np.abs(rhoa.data)
        arhoa[rhoa.mask] = np.nanmedian(arhoa)
        data = self.data.copy()
        data["rhoa"] = arhoa
        data.estimateError()
        data["err"][rhoa.mask] = 1e8
        self.data = data
        return data

    def createMesh(self, **kwargs):
        """Generate inversion mesh."""
        self.mesh = ert.createInversionMesh(self.data, **kwargs)
        self.mgr.setMesh(mesh=self.mesh)
        if kwargs.pop("show", False):
            print(self.mesh)
            pg.show(self.mesh, markers=True, showMesh=True)

    def invert(self, t=None, reg=None, regTL=None, **kwargs):
        """Run inversion for a specific timestep or all subsequently.

        Parameter
        ---------
        t : int|datetime|str|array
            time index, string or datetime object, or array of any of these
        reg : dict
            regularization options (setRegularization) for all inversions
        regTL : dict
            regularization options for timesteps inversion only
        **kwargs : dict
            keyword arguments passed to ERTManager.invert
        """
        if t is not None:
            t = self.timeIndex(t)

        if self.mesh is None:
            self.createMesh()
        self.mgr.fop.setVerbose(False)
        if isinstance(reg, dict):
            self.mgr.inv.setRegularization(**reg)
        if t is None:  # all
            t = np.arange(len(self.times))

        t = np.atleast_1d(t)
        models = []
        responses = []
        self.chi2s = []
        startModel = kwargs.pop("startModel", 100)
        creep = kwargs.pop("creep", False)
        kwargs.setdefault("isReference", True)
        for i, ti in enumerate(np.atleast_1d(t)):
            self.mgr.setData(self.chooseTime(ti))
            self.model = self.mgr.invert(startModel=startModel, **kwargs)
            if i == 0 or creep:
                startModel = self.model.copy()

            models.append(self.model)
            responses.append(self.mgr.inv.response)
            self.chi2s.append(self.mgr.inv.chi2())
            if i == 0 and isinstance(regTL, dict):
                kwargs.update(regTL)
                # self.mgr.inv.setRegularization(**regTL)

        if len(t) == 1:
            self.mgr.showResult()

        self.models = np.array(models)
        self.responses = np.array(responses)
        self.pd = self.mgr.paraDomain

    def fullInversion(self, scalef=1.0, **kwargs):
        """Full (4D) inversion."""
        DATA = [self.chooseTime(ti) for ti in range(len(self.times))]
        fop = pg.frameworks.MultiFrameModelling(ert.ERTModelling, scalef=scalef)
        fop.setData(DATA)
        if self.mesh is None:
            self.createMesh()

        fop.setMesh(self.mesh)
        print(fop.mesh())  # important to call mesh() function once!
        dataVec = np.concatenate([data["rhoa"] for data in DATA])
        errorVec = np.concatenate([data["err"] for data in DATA])
        startModel = fop.createStartModel(dataVec)
        self.inv = pg.Inversion(fop=fop, startModel=startModel, verbose=True)
        self.inv.modelTrans = 'log'
        fop.createConstraints(C=kwargs.pop("C", None))
        kwargs.setdefault("maxIter", 10)
        kwargs.setdefault("verbose", True)
        kwargs.setdefault("startModel", startModel)
        model = self.inv.run(dataVec, errorVec, **kwargs)
        self.models = np.reshape(model, [len(DATA), -1])
        self.responses = np.reshape(self.inv.response, [DATA[0].size(), -1])
        misfits = np.log(self.responses) - np.log(self.DATA)
        misfits /= np.reshape(self.data["err"], [-1, 1])
        self.chi2s = np.mean(misfits**2, axis=0)
        self.pd = fop.paraDomain
        return model

    def showFit(self, **kwargs):
        """Show data, model response and misfit."""
        _, ax = pg.plt.subplots(nrows=3, figsize=(10, 6), sharex=True, sharey=True)
        kwargs.setdefault("verbose", False)
        _, cb = self.showData(ax=ax[0], **kwargs)
        self.showData(self.mgr.inv.response, ax=ax[1],
                      cMin=cb.vmin, cMax=cb.vmax, **kwargs)
        misfit = self.mgr.inv.response / self.data["rhoa"] * 100 - 100
        self.showData(misfit, ax=ax[2], cMin=-10,
                      cMax=10, cMap="bwr", **kwargs)
        return ax

    def showAllModels(self, ncols=2, **kwargs):
        """Show all models as subplots.

        Parameters
        ----------
        ncols : int [2]
            number of columns
        **kwargs : dict
            keyword arguments passed to pg.show
        """
        nT = self.DATA.shape[1]
        showRatio = kwargs.pop("ratio", False)
        nrows = int(np.ceil(nT/ncols))
        fig, ax = pg.plt.subplots(nrows=nrows, ncols=ncols,
                                sharex=True, sharey=True,
                                figsize=kwargs.pop("figsize", [8, 5]))
        ratiokw = kwargs.copy()
        kwargs.setdefault("cMin", np.min(self.models))
        kwargs.setdefault("cMax", np.max(self.models))
        kwargs.setdefault("cMap", "Spectral_r")
        kwargs.setdefault("logScale", True)
        models = self.models.copy()
        if showRatio:
            models = self.models / self.models[0]
            rmax = np.maximum(np.max(models), 1/np.min(models))
            ratiokw["cMax"] = kwargs.pop("rMax", rmax)
            ratiokw["cMin"] = 1 / ratiokw["cMax"]
            ratiokw["cMap"] = "bwr"
            ratiokw["logScale"] = True

        for i, model in enumerate(models):
            if i == 0 or not showRatio:
                pg.show(self.pd, self.models[i], ax=ax.flat[i], **kwargs)
            else:
                pg.show(self.pd, model, ax=ax.flat[i], **ratiokw)

        fig.tight_layout()
        return ax

    def generateModelPDF(self, **kwargs):
        """Generate a multi-page pdf with the model results.

        Parameters
        ----------
        **kwargs : keyword arguments passed to pg.show()
        """
        kwargs.setdefault('label', pg.unit('res'))
        kwargs.setdefault('cMap', pg.utils.cMap('res'))
        kwargs.setdefault('logScale', True)
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(self.name+'-model.pdf') as pdf:
            fig = pg.plt.figure(figsize=kwargs.pop("figsize", [8, 5]))
            for i, model in enumerate(self.models):
                ax = fig.subplots()
                pg.show(self.pd, model, ax=ax, **kwargs)
                ax.set_title(str(i)+": " + self.times[i].isoformat(" ", "minutes"))
                fig.savefig(pdf, format='pdf')
                fig.clf()

    def generateRatioPDF(self, **kwargs):
        """Generate a multi-page pdf with the model results.

        Parameters
        ----------
        creep : bool [False]
            Use preceding time step as reference (default is baseline)
        cMax : float [2]
            maximum of color scale
        cMax : float [1/cMax]
            minimum of color scale
        logScale : bool [True]
            logarithmic color scale
        cMap : str ['bwr']
            colormap
        """
        kwargs.setdefault('label', 'ratio')
        kwargs.setdefault('cMap', 'bwr')
        kwargs.setdefault('logScale', True)
        kwargs.setdefault("cMax", 2.0)
        kwargs.setdefault("cMin", 1/kwargs["cMax"])
        basemodel = self.models[0]
        creep = kwargs.pop("creep", False)
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(self.name+'-ratio.pdf') as pdf:
            fig = pg.plt.figure(figsize=kwargs.pop("figsize", [8, 5]))
            for i, model in enumerate(self.models[1:]):
                ax = fig.subplots()
                pg.show(self.pd, model/basemodel, ax=ax, **kwargs)
                ax.set_title(str(i)+": " + self.times[i+1].isoformat(" ", "minutes") + "/" +
                             self.times[i].isoformat(" ", "minutes"))
                fig.savefig(pdf, format='pdf')
                fig.clf()
                if creep:
                    basemodel = model

    def exportVTK(self, name=None, oneforall=False):
        """Generate output vtk(s) for postprocessing."""
        name = name or self.name
        if name.endswith(".vtk"):
            name = name[:-4]

        vtk = self.pd.copy()
        if oneforall:
            for i, model in enumerate(self.models):
                vtk[f"model{i}"] = model

            vtk.exportVTK(name+"_results.vtk")
        else:
            for i, model in enumerate(self.models):
                vtk["resistivity"] = model
                vtk.exportVTK(name+f"_result{i}.vtk")

    def saveResults(self, basename=None):
        """Save inversion results."""
        if basename is None:
            basename = self.name
        
        self.pd.save(basename+".bms")
        np.savetxt(basename+".result", self.models, fmt="%.1f")

    def loadResults(self, basename=None):
        """Load inversion results."""
        if basename is None or basename is True:
            basename = self.name
        
        self.pd = pg.load(basename+".bms")
        self.models = np.loadtxt(basename+".result")


if __name__ == "__main__":
    pass
