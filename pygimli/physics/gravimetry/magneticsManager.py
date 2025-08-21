#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Method Manager for Magnetics."""
import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.viewer import pv
from pygimli.frameworks import MeshMethodManager
from .MagneticsModelling import MagneticsModelling
from .tools import depthWeighting


class MagManager(MeshMethodManager):
    """ Magnetics Manager.
    """
    def __init__(self, data=None, **kwargs):
        """ Create Magnetics Manager instance.
        """
        self.DATA = kwargs.pop("DATA", None)
        self.x = kwargs.pop("x", None)
        self.y = kwargs.pop("y", None)
        self.z = kwargs.pop("z", None)
        self.igrf = kwargs.pop("igrf", None)
        self.mesh_ = kwargs.pop("mesh", None)
        self.cmp = kwargs.pop("cmp", None)
        self.dem = kwargs.pop("dem", None)
        self.line = None

        # self.inv_ = pg.frameworks.Inversion()
        if isinstance(self.dem, str):
            from pygimli.utils import DEM
            self.dem = DEM(self.dem)

        if isinstance(data, str):
            self.DATA = np.genfromtxt(data, names=True)
            self.x = self.DATA["x"]
            self.y = self.DATA["y"]
            self.z = self.DATA["z"]
            if self.cmp is None:
                self.cmp = [t for t in self.DATA.dtype.names
                            if t.startswith("B") or t.startswith("T")]
            if self.igrf is None:
                import utm
                lat, lon = utm.to_latlon(np.mean(self.x), np.mean(self.y), 33, 'U')
                pg.info(f"Center of data: {lat}, {lon}")
                self.igrf = [lat, lon]

        super().__init__(**kwargs)
        if self.mesh_ is not None:
            if isinstance(self.mesh_, str):
                self.mesh_ = pg.load(self.mesh_)

            self.setMesh(self.mesh_)


    def __repr__(self):
        """String representation of the Magnetics Manager."""
        out = f"MagManager with {len(self.x)} data points."
        out += f"\n{len(self.cmp)} active components: "+", ".join(self.cmp)
        if self.DATA is not None:
            out += f"\nData fields: {self.DATA.dtype.names}"
        if self.mesh_ is not None:
            out += f"\nMesh: {self.mesh_.cellCount()} cells, {self.mesh_.nodeCount()} nodes."
        if self.igrf is not None:
            out += f"\nIGRF: {self.igrf}"
        if self.dem is not None:
            out += f"\nDEM: {self.dem.__repr__()}"

        return out

    def showData(self, cmp=None, **kwargs):
        """ Show data.
        """
        cmp = cmp or self.cmp
        nc = kwargs.pop("ncols", max(2 if len(cmp) > 1 else 1, len(cmp)))
        nr = (len(cmp)+nc-1) // nc
        fig, ax = pg.plt.subplots(nr, nc, sharex=True, sharey=True,
                                  squeeze=False, figsize=(7, len(self.cmp)*1+3))
        axs = np.atleast_1d(ax.flat)
        kwargs.setdefault("cmap", "bwr")
        ori = kwargs.pop("orientation", "horizontal")
        for i, c in enumerate(cmp):
            fld = self.DATA[c]
            vv = max(-np.min(fld)*1., np.max(fld)*1.)
            sc = axs[i].scatter(self.x, self.y, c=fld,
                                vmin=-vv, vmax=vv, **kwargs)
            axs[i].set_title(c)
            axs[i].set_aspect(1.0)
            fig.colorbar(sc, ax=ax.flat[i], orientation=ori)

        return ax

    def detectLines(self, **kwargs):
        """Detect lines in data.

        Keyword arguments
        -----------------
        mode: str|float|array
            'x'/'y': along coordinate axis
            spacing vector: by given spacing
            float: minimum distance
        axis: str='x'
            Axis to use for line detection.
        show: bool=False
            Show detected lines.
        """
        from pygimli.utils import detectLines
        if self.x is None or self.y is None:
            pg.error("No x and y coordinates available for line detection.")
            return

        self.line = detectLines(self.x, self.y, **kwargs)
        return self.line


    def showLineData(self, line=None, cmp=None, **kwargs):
        """ Show data for a specific line.

        Parameters
        ----------
        line: int|array
            Line number or array of line numbers to show.
        cmp: list
            List of components to show.
        """
        if cmp is None:
            cmp = self.cmp

        for l in np.atleast_1d(line):
            x = self.x[self.line==l]
            y = self.y[self.line==l]
            t = np.hstack([0, np.cumsum(np.sqrt(np.diff(x)**2+ np.diff(y)**2))])
            for c in cmp:
                pg.plt.plot(t, self.DATA[c][self.line==l], label=c+f" (line {l})")

        pg.plt.xlabel("X")
        pg.plt.ylabel("Data")
        pg.plt.legend()

    def createGrid(self, dx:float=50, depth:float=800, bnd:float=0):
        """ Create a grid.

        TODO
        ----
            * check default values, make them more sensible
            and depending on data

        Arguments
        ---------
        dx: float=50
            Grid spacing in x and y direction.
        depth: float=800
            Depth of the grid in z direction.
        bnd: float=0
            Boundary distance to extend the grid in x and y direction.

        Returns
        -------
        mesh: :gimliapi:`GIMLI::Mesh`
            Created 3D structured grid.
        """
        x = np.arange(min(self.x)-bnd, max(self.x)+bnd+.1, dx)
        y = np.arange(min(self.y)-bnd, max(self.y)+bnd+.1, dx)
        z = np.arange(-depth, .1, dx)
        self.mesh_ = mt.createGrid(x=x, y=y, z=z)
        self.fop.setMesh(self.mesh_)
        return self.mesh_


    def createMesh(self, bnd:float=0, area:float=1e5, depth:float=800,
                   quality:float=1.3, addPLC:pg.Mesh=None, addPoints:bool=True):
        """ Create an unstructured 3D mesh.

        TODO
        ----
            * check default values, make them more sensible
            and depending on data

        Arguments
        ---------
        bnd: float=0
            Boundary distance to extend the mesh in x and y direction.
        area: float=1e5
            Maximum area constraint for cells.
        depth: float=800
            Depth of the mesh in z direction.
        quality: float=1.3
            Quality factor for mesh generation.
        addPLC: :gimliapi:`GIMLI::Mesh`
            PLC mesh to add to the mesh.
        addPoints: bool=True
            Add points from self.x and self.y to the mesh.

        Returns
        -------
        mesh: :gimliapi:`GIMLI::Mesh`
            Created 3D unstructured mesh.
        """
        geo = mt.createCube(start=[min(self.x)-bnd, min(self.x)-bnd, -depth],
                            end=[max(self.x)+bnd, max(self.y)+bnd, 0])
        if addPoints is True:
            for xi, yi in zip(self.x, self.y):
                geo.createNode([xi, yi, 0])
        if addPLC:
            geo += addPLC

        self.mesh_ = mt.createMesh(geo, quality=quality, area=area)
        self.fop.setMesh(self.mesh_)
        return self.mesh_


    def createForwardOperator(self, verbose=False):
        """ Create forward operator (computationally extensive!).
        """
        # points = np.column_stack([self.x, self.y, -np.abs(self.z)])
        points = np.column_stack([self.x, self.y, self.z])
        self.fwd = MagneticsModelling(points=points,
                                      cmp=self.cmp, igrf=self.igrf, verbose=verbose)
        return self.fwd


    def inversion(self, noise_level=2, noisify=False, **kwargs):
        """Run Inversion (requires mesh and FOP).

        Arguments
        ---------
        noise_level: float|array
            absolute noise level (absoluteError)
        noisify: bool
            add noise before inversion
        relativeError: float|array [0.01]
            relative error to stabilize very low data
        depthWeighting: bool [True]
            apply depth weighting after Li&Oldenburg (1996)
        z0: float
            skin depth for depth weighting
        mul: array
            multiply constraint weight with

        Keyword arguments
        -----------------
        startModel: float|array=0.001
            Starting model (typically homogeneous)
        relativeError: float=0.001
            Relative error to stabilize very low data.
        lam: float=10
            regularization strength
        verbose: bool=True
            Be verbose
        symlogThreshold: float [0]
            Threshold for symlog data trans.
        limits: [float, float]
            Lower and upper parameter limits.
        C: int|Matrix|[float, float, float] [1]
            Constraint order.
        C(,cType): int|Matrix|[float, float, float]=C
            Constraint order, matrix or correlation lengths.
        z0: float=25
            Skin depth for depth weighting.
        depthWeighting: bool=True
            Apply depth weighting after Li&Oldenburg (1996).
        mul: float=1
            Multiply depth weighting constraint weight with this factor.
        **kwargs:
            Additional keyword arguments for the inversion.

        Returns
        -------
        model: np.array
            Model vector (also saved in self.inv.model).
        """
        dataVec = np.concatenate([self.DATA[c] for c in self.cmp])
        if noisify:
            dataVec += np.random.randn(len(dataVec)) * noise_level

        # self.inv_ = pg.Inversion(fop=self.fwd, verbose=True)
        self.inv.setForwardOperator(self.fwd)
        kwargs.setdefault("startModel", 0.001)
        kwargs.setdefault("relativeError", 0.001)
        kwargs.setdefault("lam", 10)
        kwargs.setdefault("verbose", True)

        thrs = kwargs.pop("symlogThreshold", 0)
        if thrs > 0:
            self.inv.dataTrans = pg.trans.TransSymLog(thrs)

        limits = kwargs.pop("limits", [0, 0.1])
        self.inv.setRegularization(limits=limits)
        C = kwargs.pop("C", 1)
        cType = kwargs.pop("cType", C)

        if hasattr(C, "__iter__"):
            self.inv.setRegularization(correlationLengths=C)
            cType = -1
        elif isinstance(C, pg.core.MatrixBase):
            self.inv.setRegularization(C=C)
        else:
            self.inv.setRegularization(cType=C)

        z0 = kwargs.pop("z0", 25)  # Oldenburg&Li(1996)
        dw = kwargs.pop("depthWeighting", True)
        if np.any(dw):
            pg.info("Using depth")
            if dw is True:
                pg.info("Compute depth weighting with z0 = ", z0)
                dw = depthWeighting(self.mesh_, cell=not(cType==1), z0=z0)

            cw = self.fwd.regionManager().constraintWeights()
            if len(cw) > 0 and len(dw) == len(cw):
                # dw *= cw
                print(min(dw), max(dw))
            else:
                print("lengths not matching!")

            dw *= kwargs.pop("mul", 1)
            self.inv.setConstraintWeights(dw)

        model = self.inv.run(dataVec, absoluteError=noise_level, **kwargs)
        return model


    def showDataFit(self, **kwargs):
        """ Show data, model response and misfit.

        Keyword arguments
        -----------------
        cmap : str ['bwr']
            colormap
        maxField : float
            colorbar for field
        maxMisfit : float
            colorbar for (error-weighted) misfit    
        vmin/vmax
        """
        nc = len(self.cmp)
        _, ax = pg.plt.subplots(ncols=3, nrows=nc, figsize=(12, 3*nc),
                                sharex=True, sharey=True, squeeze=False)
        vals = np.reshape(self.inv.dataVals, [nc, -1])
        resp = np.reshape(self.inv.response, [nc, -1])
        errs = np.reshape(self.inv.errorVals, [nc, -1])  # relative!
        misf = (vals - resp) / np.abs(errs *  vals)
        fkw = kwargs.copy()
        fkw.setdefault('cmap', "bwr")
        mm = fkw.pop('maxField', np.max(np.abs(vals)))
        fkw.setdefault('vmin', -mm)
        fkw.setdefault('vmax', mm)
        mkw = kwargs.copy()
        mmis = fkw.pop('maxMisfit', 3)
        mkw.setdefault('cmap', "bwr")
        mkw.setdefault('vmin', -mmis)
        mkw.setdefault('vmax', mmis)
        for i in range(nc):
            ax[i, 0].scatter(self.x, self.y, c=vals[i], **fkw)
            ax[i, 1].scatter(self.x, self.y, c=resp[i], **fkw)
            ax[i, 2].scatter(self.x, self.y, c=misf[i], **mkw)


    def show3DModel(self, label:str=None, trsh:float=0.025,
                    synth:pg.Mesh=None, invert:bool=False,
                    position:str="yz", elevation:float=10, azimuth:float=25,
                    zoom:float=1.2, **kwargs):
        """ Standard 3D view.

        Arguments
        ---------
        label: str='sus'
            Label for the mesh data to visualize.

        trsh: float=0.025
            Threshold for the mesh data to visualize.
        synth: :gimliapi:`GIMLI::Mesh`=None
            Synthetic model to visualize in wireframe.
        invert: bool=False
            Invert the threshold filter.
        position: str="yz"
            Camera position, e.g. "yz", "xz", "xy".
        elevation: float=10
            Camera elevation angle.
        azimuth: float=25
            Camera azimuth angle.
        zoom: float=1.2
            Camera zoom factor.

        Keyword arguments
        -----------------
        cMin: float=0.001
            Minimum color value for the mesh data.
        cMax: float=None
            Maximum color value for the mesh data. If None, it is set to the
            maximum value of the mesh data.
        cMap: str="Spectral_r"
            Colormap for the mesh data visualization.
        logScale: bool=False
            Use logarithmic scale for the mesh data visualization.
        **kwargs:
            Additional keyword arguments for the pyvista plot.

        Returns
        -------
        pl: pyvista Plotter
            Plot widget with the 3D model visualization.
        """
        if label is None:
            ## ever happen that this is a string?
            label = self.inv.model

        if not isinstance(label, str):
            self.mesh_["sus"] = np.array(label)
            label = "sus"

        kwargs.setdefault("cMin", 0.001)
        kwargs.setdefault("cMax", max(self.mesh_[label]))
        kwargs.setdefault("cMap", "Spectral_r")
        kwargs.setdefault("logScale", False)

        flt = None
        pl, _ = pg.show(self.mesh_, style="wireframe", hold=True,
                        alpha=0.1)
        # mm = [min(self.mesh_[label]), min(self.mesh_[label])]
        if trsh > 0:
            flt = {"threshold": {'value':trsh, 'scalars':label,'invert':invert}}
            pv.drawModel(pl, self.mesh_, label=label, style="surface",
                        filter=flt, **kwargs)

        pv.drawMesh(pl, self.mesh_, label=label, style="surface", **kwargs,
                    filter={"slice": {'normal':[-1, 0, 0], 'origin':[0, 0, 0]}})

        if synth:
            pv.drawModel(pl, synth, style="wireframe")

        pl.camera_position = position
        pl.camera.azimuth = azimuth
        pl.camera.elevation = elevation
        pl.camera.zoom(zoom)
        pl.show()
        return pl


if __name__ == "__main__":
    pass
