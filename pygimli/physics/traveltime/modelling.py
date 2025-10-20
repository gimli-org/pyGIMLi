#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Modelling classes for managing first arrival travel-time problems"""

import numpy as np
import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli.viewer.mpl import createColorBar  # , updateColorBar
from .utils import createGradientModel2D, shotReceiverDistances
from .plotting import drawVA


class TravelTimeDijkstraModelling(MeshModelling):
    """Forward modelling class for traveltime using Dijktras method."""

    def __init__(self, **kwargs):
        secNodes = kwargs.pop("secNodes", 3)
        super().__init__(**kwargs)

        self._core = pg.core.TravelTimeDijkstraModelling()
        self._core.setRegionManager(self.regionManagerRef())

        self._useGradient = None  # assumed to be [vTop, vBot] if set
        self._refineSecNodes = secNodes
        # self._refineSecNodes = kwargs.pop("secNodes", 3)  # inactive!
        self.jacobian = self._core.jacobian
        self.setThreadCount = self._core.setThreadCount
        # self.createJacobian = self.dijkstra.createJacobian
        self.setJacobian(self._core.jacobian())

    @property
    def dijkstra(self):
        """Return current Dijkstra graph associated to mesh and model."""
        return self._core.dijkstra()

    # def regionManagerRef(self):
    #     """Region manager reference (core Dijkstra has an own!)."""
    #     return self._core.regionManagerRef()

    def createRefinedFwdMesh(self, mesh):
        """Refine the current mesh for higher accuracy.

        This is called automatic when accesing self.mesh() so it ensures any
        effect of changing region properties (background, single).
        """
        pg.info("Creating refined mesh (secnodes: {0}) to "
                "solve forward task.".format(self._refineSecNodes))
        self.meshNoSec = pg.Mesh(mesh)
        m = mesh.createMeshWithSecondaryNodes(self._refineSecNodes)
        pg.verbose(m)
        return m

    def setMeshPost(self, mesh):
        """Set mesh after forward operator has been initalized."""
        self._core.setMesh(mesh, ignoreRegionManager=True)

    def setDataPost(self, data):
        """Set data after forward operator has been initalized."""
        self._core.setData(data)

    def createGraph(self, slowness):
        """Create Dijkstra graph."""
        return self._core.createGraph(slowness)

    def createStartModel(self, dataVals):
        """Create a starting model from data values (gradient or constant)."""
        sm = None

        if self._useGradient is not None:
            [vTop, vBot] = self._useGradient  # something strange here!!!
            pg.info('Create gradient starting model. {0}: {1}'.format(vTop,
                                                                      vBot))
            sm = createGradientModel2D(self.data,
                                       self.paraDomain,
                                       vTop, vBot)
        else:
            dists = shotReceiverDistances(self.data, full=True)
            aSlow = 1. / (dists / dataVals)

            # pg._r(self.regionManager().parameterCount())
            sm = pg.Vector(self.regionManager().parameterCount(),
                           pg.math.median(aSlow))
            pg.info('Create constant starting model:', sm[0])

        return sm

    def createJacobian(self, par):
        """Create Jacobian (way matrix)."""
        if not self.mesh():
            pg.critical("no mesh")

        return self._core.createJacobian(par)

    def response(self, par):
        """Return forward response (simulated traveltimes)."""
        if not self.mesh():
            pg.critical("no mesh")
        return self._core.response(par)

    def way(self, s, g):
        """Return node indices for the way from the shot to the receiver.

        The index is based on the given data, mesh and last known model.
        """
        return self._core.way(s, g)

    def drawModel(self, ax, model, **kwargs):
        """Draw the model."""
        kwargs.setdefault('label', pg.unit('vel'))
        kwargs.setdefault('cMap', pg.utils.cMap('vel'))

        return super().drawModel(ax=ax, model=model,
                                 logScale=kwargs.pop('logScale', True),
                                 **kwargs)

    def drawData(self, ax, data=None, **kwargs):
        """Draw the data (as apparent velocity crossplot by default).

        Parameters
        ----------
        data: pg.DataContainer
        """
        if hasattr(data, '__iter__'):
            kwargs['vals'] = data
            data = self.data
        elif data is None:
            data = self.data

        if kwargs.pop('firstPicks', False):
            pg.physics.traveltime.drawFirstPicks(ax, data, **kwargs)
            return ax
        else:
            kwargs.setdefault('usePos', False)
            gci = drawVA(ax, data, **kwargs)
            cBar = createColorBar(gci, **kwargs)

            return gci, cBar


class FatrayDijkstraModellingInterpolate(TravelTimeDijkstraModelling):
    """Shortest-path (Dijkstra) based travel time with fat ray jacobian."""

    def __init__(self, frequency=100., **kwargs):
        super().__init__(**kwargs)
        self.frequency = frequency
        self.iMat = pg.matrix.SparseMapMatrix()
        self.J = pg.Matrix()
        self.setJacobian(self.J)
        self._core.setJacobian(self.J)
        self.sensorNodes = None

    def createJacobian(self, slowness):
        """Generate Jacobian matrix using fat-ray after Jordi et al. (2016)."""
        # mesh = self.mesh()
        mesh = self.meshNoSec  # change back with pgcore=1.5
        self.J.resize(self.data.size(), mesh.cellCount())
        self.sensorNodes = [mesh.findNearestNode(pos)
                            for pos in self.data.sensorPositions()]
        if (self.iMat.cols() != mesh.nodeCount() or
                self.iMat.rows() != mesh.cellCount()):
            self.iMat = mesh.interpolationMatrix(mesh.cellCenters())

        Di = self.dijkstra
        slowPerCell = self.createMappedModel(slowness, 1e16)
        Di.setGraph(self._core.createGraph(slowPerCell))
        numN = mesh.nodeCount()
        data = self.data
        numS = data.sensorCount()
        Tmat = pg.Matrix(numS, numN)
        Dmat = pg.Matrix(numS, numS)
        for i, node in enumerate(self.sensorNodes):
            Di.setStartNode(node)
            Tmat[i] = Di.distances()[:numN]  # change back with pgcore=1.5
            Dmat[i] = Tmat[i][self.sensorNodes]

        self.FresnelWeight = pg.Matrix(data.size(), len(slowness))
        for i in range(data.size()):
            iS = int(data("s")[i])
            iG = int(data("g")[i])
            tsr = Dmat[iS][iG]  # shot-receiver travel time
            dt = self.iMat * (Tmat[iS] + Tmat[iG]) - tsr
            weight = np.maximum(1 - 2 * self.frequency * dt, 0.0)  # 1 on ray
            wa = weight  # * np.sqrt(self.mesh().cellSizes())
            if np.sum(wa) > 0:  # not if all values are zero
                wa /= np.sum(wa)

            self.J[i] = wa * tsr / slowness
            self.FresnelWeight[i] = wa

        self.setJacobian(self.J)
        self._core.setJacobian(self.J)


class FatrayDijkstraModellingMidpoint(TravelTimeDijkstraModelling):
    """Shortest-path (Dijkstra) based travel time with fat ray jacobian."""

    def __init__(self, frequency=100., verbose=False):
        super().__init__(verbose)
        self.frequency = frequency
        self.mids = None
        self.J = pg.Matrix()
        self.sensorNodes = None
        self.setJacobian(self.J)
        self._core.setJacobian(self.J)

    def setMesh(self, mesh, **kwargs):  # secondaryNodes=3):
        """Set mesh and create additional secondary Nodes in cell centers."""
        super().setMesh(mesh, **kwargs)  # ignoreRegionManager=True)
        print(self.mesh(), self.mesh().secondaryNodeCount())
        self.mids = pg.IVector()
        for c in self.mesh().cells():
            n = self.mesh().createSecondaryNode(c.center())
            c.addSecondaryNode(n)
            self.mids.push_back(n.id())

        print(self.mesh())

    def createJacobian(self, slowness):
        """Generate Jacobian matrix using fat-ray after Jordi et al. (2016)."""
        self.J.resize(self.data.size(), self.mesh().cellCount())
        self.sensorNodes = [self.mesh.findNearestNode(pos)
                            for pos in self.data.sensorPositions()]
        Di = self.dijkstra()
        slowPerCell = self.createMappedModel(slowness, 1e16)
        Di.setGraph(self.createGraph(slowPerCell))
        numN = self.mesh.nodeCount()
        data = self.data
        numS = data.sensorCount()
        Tmat = pg.Matrix(numS, numN)
        Dmat = pg.Matrix(numS, numS)
        pg.debug(self.mesh())
        pg.debug(self.mesh().nodeCount(), max(self.mids))
        for i, node in enumerate(self.sensorNodes):
            Di.setStartNode(node)
            dist0 = Di.distances()
            dist = Di.distances(withSecNodes=True)
            print("dist len ", len(dist0), len(dist))
            Tmat[i] = dist[self.mids]
            Dmat[i] = Tmat[i][self.sensorNodes]

        for i in range(data.size()):
            iS = int(data("s")[i])
            iG = int(data("g")[i])
            tsr = Dmat[iS][iG]  # shot-receiver travel time
            dt = Tmat[iS] + Tmat[iG] - tsr
            weight = np.maximum(1 - 2 * self.frequency * dt, 0.0)  # 1 on ray
            wa = weight  # * np.sqrt(self.mesh().cellSizes())
            if np.sum(wa) > 0:  # not if all values are zero
                wa /= np.sum(wa)

            self.J[i] = wa * tsr / slowness

        self.setJacobian(self.J)
        self._core.setJacobian(self.J)


FatrayDijkstraModelling = FatrayDijkstraModellingInterpolate
# FatrayDijkstraModelling = FatrayDijkstraModellingMidpoint
