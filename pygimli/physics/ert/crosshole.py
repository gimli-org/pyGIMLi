"""Crosshole ERT inversion."""
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import pygimli as pg
import pygimli.meshtools as mt
# from pygimli.physics import ert
from .timelapse import TimelapseERT


class CrossholeERT(TimelapseERT):
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
        self.bhmap = kwargs.pop("bhmap", None)
        super().__init__(filename=filename, **kwargs)
        if self.bhmap is None:
            self.determineBHmap()

    def __repr__(self):  # for print function
        """Return string representation of the class."""
        mys = 'Crosshole '
        if self.bhmap is not None:
            mys += '({:d} boreholes) '.format(len(np.unique(self.bhmap)))
        return mys + super().__repr__()

    def determineBHmap(self):  # XH-specific
        """Auto-determine borehole map from xy positions."""
        xy = pg.x(self.data)*999+pg.y(self.data)*999  # maybe needs a data.sortX() first
        d0 = np.diff(np.hstack([0, np.nonzero(np.diff(xy))[0]+1, len(xy)]))
        self.bhmap = np.concatenate([np.ones(dd, dtype=int)*i for i, dd in enumerate(d0)])

    def load(self, filename, **kwargs):
        """Load or import data."""  # TL-ERT
        self.bhmap = kwargs.pop("bhmap", None)
        super().load(filename, **kwargs)
        if self.bhmap is None:
            self.determineBHmap()

        for tok in "abmn":
            self.data["n"+tok] = self.bhmap[self.data[tok]]

    def showData(self, v="rhoa", x="a", y="m", t=None, **kwargs):
        """Show data.

        Show data as pseudosections (single-hole) or cross-plot (crosshole)

        Parameters
        ----------
        v : str|array ["rhoa]
            array or field to plot
        x, y : str|array ["a", "m"]
            values to use for x and y axes
        t : int|str|datetime
            time to choose
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

        if len(np.unique(self.bhmap)) == 1 or "style" in kwargs:
            return self.data.show(v, **kwargs)
        else:
            ax, cb = pg.viewer.mpl.showDataContainerAsMatrix(
                self.data, x=x, y=y, v=v, **kwargs)
            xx = np.nonzero(np.diff(self.bhmap))[0] + 1
            ax.set_xticks(xx)
            ax.set_yticks(xx)
            return ax, cb

    def extractSubset(self, nbh, plane=None, name=None):
        """Extract a subset (slice) by borehole number.

        Returns a CrossholeERT instance with reduced boreholes

        Parameters
        ----------
        nbh : int|array
            borehole(s) to extract data from
        name : str
            name to give the new instance
        plane : bool [None]
            reduce to 2D (automatic if max. 2 boreholes are used)
        """
        xh2 = pg.DataContainerERT(self.data)
        good = pg.Vector(xh2.size(), 1)
        good = np.ones(xh2.size(), dtype=bool)
        if isinstance(nbh, int):
            nbh = [nbh]
        for tok in ["a", "b", "m", "n"]:
            bla = np.zeros(xh2.size(), dtype=bool)
            for nn in nbh:
                bla = bla | (self.data["n"+tok] == nn - 1)

            good = good & bla

        xh2["valid"] *= 0  # all invalid
        xh2.markValid(np.nonzero(good)[0])
        xh2.removeInvalid()
        xh2.removeUnusedSensors()
        if plane is None:  # decide upon number
            plane = len(nbh) < 3
        if plane:
            p0 = self.data.sensor(np.nonzero(self.bhmap==nbh[0]-1)[0][0])
            dists = [p0.dist(self.data.sensor(np.nonzero(
                self.bhmap==nn-1)[0][0])) for nn in nbh]
            for i in range(xh2.sensorCount()):
                xh2.setSensor(i, [dists[self.bhmap[i]], 0,
                                  xh2.sensorPosition(i).z()])
        name = name or self.name + "xh" + "".join([str(i) for i in nbh])
        return CrossholeERT(data=xh2, DATA=self.DATA[good], times=self.times,
                            name=name)

    def createMesh(self, ibound=2, obound=10, quality=None, show=False,
                   threeD=None, ref=0.25, area=1):
        """Create a 2D mesh around boreholes.

        Parameters
        ----------
        ibound : float
            inner boundary in m
        ibound : float
            outer boundary in m
        quality : float
            triangle or tetgen quality
        threeD : bool|None
            create 3D model (None-automatic)
        ref : float
            electrode refinement in m
        """
        data = self.data
        xmin, xmax = min(pg.x(data)), max(pg.x(data))
        ymin, ymax = min(pg.y(data)), max(pg.y(data))
        zmin, zmax = min(pg.z(data)), max(pg.z(data))
        if threeD is None:
            threeD = (ymin != ymax) and (zmin != zmax)

        ztop = np.minimum(0, zmax+ibound)
        if threeD:
            world = mt.createWorld(start=[xmin-obound, ymin-obound, zmin-obound],
                                end=[xmax+obound, ymax+obound, 0], marker=1)
            box = mt.createCube(start=[xmin-ibound, ymin-ibound, zmin-ibound],
                                end=[xmax+ibound, ymax+ibound, ztop],
                                marker=2, area=area)
            if quality is None:
                quality = 1.3
        else:  # 2D
            zodown = zmin - obound
            zidown = zmin - ibound
            if zmin == zmax:  # z coordinate already in y
                ztop = np.minimum(0, ymax+ibound)
                zodown = ymin - obound
                zidown = ymin - ibound

            world = mt.createWorld(start=[xmin-obound, zodown],
                                   end=[xmax+obound, 0.], marker=1)
            box = mt.createRectangle(start=[xmin-ibound, zidown],
                                     end=[xmax+ibound, ztop], marker=2,
                                     area=area)
            if quality is None:
                quality = 34.4

        geo = world + box
        for pos in data.sensors():
            if threeD:
                geo.createNode(pos, marker=-99)
                geo.createNodeWithCheck(pos - pg.Pos(0, 0, ref))  # refinement
            else:
                if zmin == zmax:
                    geo.createNode([pos.x(), pos.y()], marker=-99)
                    geo.createNode([pos.x(), pos.y()-ref])
                else:
                    geo.createNode([pos.x(), pos.z()], marker=-99)
                    geo.createNode([pos.x(), pos.z()-ref])

        self.mesh = mt.createMesh(geo, quality=quality)
        self.mgr.setMesh(self.mesh)
        if show:
            pg.show(self.mesh, markers=True, showMesh=True)

    # def showFit(self, **kwargs):
    #     """Show data, model response and misfit."""
    #     _, ax = plt.subplots(ncols=3, figsize=(10, 6), sharex=True, sharey=True)
    #     _, cb = self.showData(ax=ax[0], verbose=False)
    #     self.showData(self.mgr.inv.response, ax=ax[1],
    #                   cMin=cb.vmin, cMax=cb.vmax, verbose=False)
    #     misfit = self.mgr.inv.response / self.data["rhoa"] * 100 - 100
    #     self.showData(misfit, ax=ax[2], cMin=-10, cMax=10, cMap="bwr", verbose=0)
    #     return ax



if __name__ == "__main__":
    pass
