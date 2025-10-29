#!/usr/bin/env python3
"""Class for managing first arrival travel time inversions."""
from pathlib import Path
import numpy as np

import pygimli as pg
from pygimli.frameworks import MeshMethodManager

from pygimli.utils import getSavePath
from . modelling import TravelTimeDijkstraModelling, FatrayDijkstraModelling
from . plotting import drawFirstPicks


class TravelTimeManager(MeshMethodManager):
    """Manager for refraction seismics (traveltime tomography).

    TODO Document main members and use default MethodManager interface
    e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
    """

    def __init__(self, data=None, **kwargs):
        """Create an instance of the Traveltime manager.

        Parameters
        ----------
        data: :gimliapi:`GIMLI::DataContainer` | str
            You can initialize the Manager with data or give them a dataset
            when calling the inversion.
        """
        self.useFatray = kwargs.pop("fatray", False)
        self.frequency = kwargs.pop("frequency", 100.)
        self.secNodes = kwargs.pop("secNodes", 2)

        super().__init__(data=data, **kwargs)

        self.inv.dataTrans = pg.trans.Trans()


    @property
    def velocity(self):
        """Return velocity vector (the inversion model)."""
        # we can check here if there was an inversion run
        return self.fw.model  # shouldn't it be the inverse?


    def createForwardOperator(self, **kwargs):
        """Create default forward operator for Traveltime modelling.

        Your want your Manager use a special forward operator you can add them
        here on default Dijkstra is used.
        """
        if self.useFatray:
            fop = FatrayDijkstraModelling(frequency=self.frequency, **kwargs)
        else:
            fop = TravelTimeDijkstraModelling(verbose=self.verbose)
        return fop


    def load(self, fileName):
        """Load any supported data file."""
        self.data = pg.physics.traveltime.load(fileName)
        return self.data


    def createMeshMovedToMeshManager(self, data=None, **kwargs):
        """Create default inversion mesh.

        Inversion mesh for traveltime inversion does not need boundary region.

        Args
        ----
        data: DataContainer
            Data container to read sensors from.

        Keyword Args
        ------------
        Forwarded to `:py:func:pygimli.meshtools.createParaMesh`
        """
        pg._error('do not use. will be removed soon') #25032025
        data = data or self.data

        if not hasattr(data, 'sensors'):
            pg.critical('Please provide a data container for mesh generation')

        return pg.meshtools.createParaMesh(data.sensors(), paraBoundary=0,
                                           boundary=0, **kwargs)


    def checkData(self, data):
        """Return data from container."""
        if isinstance(data, pg.DataContainer):
            if not data.haveData('t'):
                pg.critical('DataContainer has no "t" values.')
            return data['t']

        return data


    def checkError(self, err, dataVals):
        """Return relative error."""
        if isinstance(err, pg.DataContainer):
            if not err.haveData('err'):
                pg.error('DataContainer has no "err" values. Fallback to 3%')
                return np.ones(err.size()) * 0.03
            return err['err'] / dataVals

        return err


    def applyMesh(self, mesh, secNodes=None, ignoreRegionManager=False):
        """Apply mesh, i.e. set mesh in the forward operator class."""
        if secNodes is None:
            secNodes = self.secNodes

        self.fop._refineSecNodes = secNodes
        if secNodes > 0 and ignoreRegionManager:
            mesh = self.fop.createRefinedFwdMesh(mesh)

        self.fop.setMesh(mesh, ignoreRegionManager=ignoreRegionManager)
        self.fop.jacobian().clear()


    def simulate(self, mesh=None, scheme=None, slowness=None, vel=None, seed=None,
                 secNodes=2, noiseLevel=0.0, noiseAbs=0.0, **kwargs):
        """Simulate traveltime measurements.

        Perform the forward task for a given mesh, a slowness distribution (per
        cell) and return data (traveltime) for a measurement scheme.

        Parameters
        ----------
        mesh : :gimliapi:`GIMLI::Mesh`
            Mesh to calculate for or use the last known mesh.
        scheme: :gimliapi:`GIMLI::DataContainer`
            Data measurement scheme needs 's' for shot and 'g' for geophone
            data token.
        slowness : array(mesh.cellCount()) | array(N, mesh.cellCount())
            Slowness distribution for the given mesh cells can be:

            * a single array of len mesh.cellCount()
            * a matrix of N slowness distributions of len mesh.cellCount()
            * a slowness map [[marker0, slow0], [marker1, slow1], ...]
        vel : array(mesh.cellCount()) | array(N, mesh.cellCount())
            Velocity distribution for the mesh cells (overwrites slowness!).
        secNodes: int [2]
            Number of refinement nodes to increase accuracy of the forward
            calculation.
        noiseLevel: float [0.0]
            Add relative noise to the simulated data. noiseLevel*100 in %
        noiseAbs: float [0.0]
            Add absolute noise to the simulated data in ms.
        seed: int [None]
            Seed the random generator for the noise.

        Keyword Arguments
        -----------------
        returnArray: [False]
            Return only the calculated times.
        verbose: [self.verbose]
            Overwrite verbose level.
        **kwargs
            Additional kwargs ...

        Returns
        -------
        t : array(N, data.size()) | DataContainer
            The resulting simulated travel time values (returnArray=True)
            or DataContainer containing them in t field (returnArray=False).
            Either one column array or matrix in case of slowness matrix.
        """
        verbose = kwargs.pop('verbose', self.verbose)

        fop = self.fop
        scheme = scheme or self.data
        fop.data = scheme
        fop.verbose = verbose

        if mesh is not None:
            self.applyMesh(mesh, secNodes=secNodes, ignoreRegionManager=True)

        if vel is not None:
            slowness = 1/vel

        if slowness is None:
            pg.critical("Need some slowness or velocity distribution for"
                        " simulation.")

        if len(slowness) == self.fop.mesh().cellCount():
            t = fop.response(slowness)
            if verbose:
                print('min/max t:', min(t), max(t))
        else:
            print(self.fop.mesh())
            print("slowness: ", slowness)
            pg.critical("Simulate called with wrong slowness array.")

        ret = pg.DataContainer(scheme)

        if noiseLevel > 0 or noiseAbs > 0:
            if not ret.allNonZero('err'):
                err = noiseAbs + t * noiseLevel
                ret['err'] = err

            pg.verbose("Absolute error estimates (min:max) {}:{}".format(
                min(ret['err']), max(ret['err'])))

            t += pg.randn(ret.size(), seed=seed) * ret['err']

        if kwargs.pop('returnArray', False) is True:
            return t

        ret['t'] = t
        return ret


    def invert(self, data=None, useGradient=True, vTop=500, vBottom=5000,
               secNodes=None, **kwargs):
        """Invert data.

        Parameters
        ----------
        data : pg.DataContainer()
            Data container with at least SensorIndices 's g' (shot/geophone) &
            data values 't' (traveltime in s) and 'err' (absolute error in s)
        useGradient: bool [True]
            Use gradient starting model typical for refraction cases.
            For crosshole tomography geometry you should set this to False for
            a non-gradient (e.g. homogeneous) starting model.
        vTop: float
            Top velocity for gradient starting model.
        vBottom: float
            Bottom velocity for gradient starting model.
        secNodes: int [2]
            Number of secondary nodes for accuracy of forward computation.

        Returns
        -------
        model:
            Mapped (for paradomain) velocity model.

        Keyword Arguments
        -----------------
        ** kwargs:
            Inversion related arguments:
            See :py:mod:`pygimli.frameworks.MeshMethodManager.invert`
        """
        mesh = kwargs.pop('mesh', None)
        if secNodes is not None:
            self.secNodes = secNodes

        if 'limits' in kwargs and kwargs['limits'][0] > 1:
            tmp = kwargs['limits'][0]
            kwargs['limits'][0] = 1.0 / kwargs['limits'][1]
            kwargs['limits'][1] = 1.0 / tmp
            pg.verbose('Switching velocity limits to slowness limits.',
                       kwargs['limits'])

        if useGradient:
            self.fop._useGradient = [vTop, vBottom]
        else:
            self.fop._useGradient = None

        ### invert return mapped models
        slowness = super().invert(data, mesh, **kwargs)
        velocity = 1.0 / slowness
        velocity.isParaModel = slowness.isParaModel
        self.fw.model = 1.0 / self.fw.model #C42 self.fw only hold non-mapped model
        # that needs to be compatible to self.fw.mesh
        return velocity


    def showFit(self, axs=None, firstPicks=True, **kwargs):
        """Show data fit as first-break picks or apparent velocity."""
        if firstPicks:
            kwargs.setdefault("linestyle", "None")
            ax, _ = self.showData(firstPicks=True, **kwargs)
            drawFirstPicks(ax, self.fop.data, self.inv.response, marker=None)
        else:
            super().showFit(axs=axs, **kwargs)


    def getRayPaths(self, model=None):
        """Compute ray paths.

        If model is not specified, the last calculated Jacobian is used.

        Parameters
        ----------
        model : array
            Velocity model for which to calculate and visualize ray paths (the
            default is model for last Jacobian calculation in self.velocity).

        Returns
        -------
        list of two-column array holding x and y positions
        """
        if model is not None:
            # if self.fop.jacobian().size() == 0 or model != self.model:
            self.fop.createJacobian(1/model)

        shots = self.fop.data.id("s")
        recei = self.fop.data.id("g")

        segs = []
        nodes = self.fop._core.mesh().positions(withSecNodes=True)

        for s, g in zip(shots, recei, strict=False):
            wi = self.fop.way(s, g)
            points = nodes[wi]
            segs.append(np.column_stack((pg.x(points), pg.y(points))))

        return segs


    def drawRayPaths(self, ax, model=None, rayPaths=None, **kwargs):
        """Draw the the ray paths for model or last model.

        If model is not specified, the last calculated Jacobian is used.

        Parameters
        ----------
        rayPaths : list of np.array
            x/y column array with ray point positions
        model : array
            Velocity model for which to calculate and visualize ray paths (the
            default is model for last Jacobian calculation in self.velocity).
        ax : matplotlib.axes object
            To draw the model and the path into.
        **kwargs : type
            Additional arguments passed to LineCollection (alpha, linewidths,
            color, linestyles).

        Returns
        -------
        lc : matplotlib.LineCollection
        """
        rayPaths = rayPaths or self.getRayPaths(model=model)

        _ = kwargs.setdefault("color", "w")
        _ = kwargs.setdefault("alpha", 0.5)
        _ = kwargs.setdefault("linewidths", 0.8)

        from matplotlib.collections import LineCollection
        lc = LineCollection(rayPaths, **kwargs)
        ax.add_collection(lc)

        return lc


    def showRayPaths(self, model=None, ax=None, **kwargs):
        """Show the model with ray paths for given model.

        If not model specified, the last calculated Jacobian is taken.

        Parameters
        ----------
        model : array
            Velocity model for which to calculate and visualize ray paths (the
            default is model for last Jacobian calculation in self.velocity).
        ax : matplotlib.axes object
            To draw the model and the path into.
        **kwargs : type
            forward to drawRayPaths

        Returns
        -------
        ax : matplotlib.axes object
        cb : matplotlib.colorbar object (only if model is provided)

        Examples
        --------
        >>> import pygimli as pg
        >>> from pygimli.physics import traveltime as tt
        >>>
        >>> x, y = 8, 6
        >>> mesh = pg.createGrid(x, y)
        >>> data = tt.createRAData([(0, 0)] + [(x, i) for i in range(y)],
        ...                       shotDistance=y+1)
        >>> data["t"] = 1.0
        >>> mgr = tt.Manager(data)
        >>> mgr.applyMesh(mesh, secNodes=10)
        >>> ax, cb = mgr.showRayPaths(showMesh=True, diam=0.1)
        """
        if model is None:
            if self.fop.jacobian().size() == 0:
                self.fop.mesh()  # initialize any meshs .. just to be sure is 1
                model = pg.Vector(self.fop.regionManager().parameterCount(),
                                  1.0)
            else:
                model = self.model

        ax, cbar = self.showModel(ax=ax, model=model,
                                  showMesh=kwargs.pop('showMesh', None),
                                  diam=kwargs.pop('diam', None))

        self.drawRayPaths(ax, model=model, **kwargs)

        return ax, cbar


    def rayCoverage(self):
        """Ray coverage, i.e. summed raypath lengths."""
        J = self.fop.jacobian()
        return J.transMult(np.ones(J.rows()))


    def createTraveltimefield(self, v=None, startPos=None, withSec=False):
        """Compute a single traveltime field."""
        startPos = startPos or self.data.sensor(0)
        fop = self.fop
        mesh = fop.mesh()
        Di = fop.dijkstra
        slowPerCell = fop.createMappedModel(1/v, 1e16)
        Di.setGraph(fop._core.createGraph(slowPerCell))
        Di.setStartNode(mesh.findNearestNode(startPos))
        dist = Di.distances()
        if withSec:
            return dist
        else:
            return dist[:mesh.nodeCount()]


    def standardizedCoverage(self):
        """Standardized coverage vector (0|1) using neighbor info."""
        coverage = self.rayCoverage()
        C = self.fop.constraintsRef()
        return np.sign(np.absolute(C.transMult(C * coverage)))


    def showCoverage(self, ax=None, name='coverage', **kwargs):
        """Show the ray coverage in log-scale."""
        if ax is None:
            fig, ax = pg.plt.subplots()

        cov = self.rayCoverage()
        return pg.show(self.fop.paraDomain,
                       pg.log10(cov+min(cov[cov > 0])*.5), ax=ax,
                       coverage=self.standardizedCoverage(), **kwargs)


    def saveResult(self, folder=None, size=(16, 10), verbose=False, **kwargs):
        """Save the results in a specified (or date-time derived) folder.

        Saved items are:
            * Resulting inversion model
            * Velocity vector
            * Coverage vector
            * Standardized coverage vector
            * Mesh (bms and vtk with results)

        Args
        ----
        path: str[None]
            Path to save into. If not set the name is automatically created
        size: (float, float) (16,10)
            Figure size.

        Keyword Args
        ------------
        Will be forwarded to showResults

        Returns
        -------
        str:
            Name of the result path.
        """
        subfolder = self.__class__.__name__
        path = Path(getSavePath(folder, subfolder))

        if verbose:
            pg.info(f'Saving refraction data to: {path}')

        np.savetxt(path / 'velocity.vector', self.velocity)
        np.savetxt(path / 'velocity-cov.vector', self.rayCoverage())
        np.savetxt(path / 'velocity-scov.vector', self.standardizedCoverage())

        m = pg.Mesh(self.paraDomain)

        m['Velocity'] = self.paraModel(self.velocity)
        m['Coverage'] = self.rayCoverage()
        m['S_Coverage'] = self.standardizedCoverage()
        m.exportVTK(path / 'velocity')
        m.saveBinaryV2(path / 'velocity-pd')
        self.fop.mesh().save(path / 'velocity-mesh')

        np.savetxt(path / 'chi.txt', self.inv.chi2History)

        if m.dim() == 2:
            fig, ax = pg.plt.subplots()
            self.showResult(ax=ax, cov=self.standardizedCoverage(), **kwargs)
            fig.set_size_inches(size)
            fig.savefig(path / 'velocity.pdf', bbox_inches='tight')
            pg.plt.close(fig)

        return path
