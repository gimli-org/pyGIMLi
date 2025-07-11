#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic modelling proxies.
"""
import numpy as np
import pygimli as pg

DEFAULT_STYLES = {
    'Default': {
        'color': 'C0',
        'lw': 1.5,
        'linestyle': '-'
    },
    'Data': {
        'color': 'C0',  # blueish
        'lw': 1.5,
        'linestyle': ':',
        'marker': 'o',
        'ms': 6
    },
    'Response': {
        'color': 'C0',  # blueish
        'lw': 2.0,
        'linestyle': '-',
        'marker': 'None',
        'alpha': 0.4
    },
    'Error': {
        'color': 'C3',  # reddish
        'lw': 0,
        'linestyle': '-',
        'elinewidth': 2,
        'alpha': 0.5
    },
}


class Modelling(pg.core.ModellingBase):
    """Abstract Forward Operator.

    Abstract Forward Operator that is or can use a Modelling instance.
    Can be seen as some kind of proxy Forward Operator.

    Todo
    ----
        * Doc:
            - describe members (model transformation,
                                dictionary of region properties)
        * think about splitting all mesh related into MeshModelling
        * clarify difference: setData(array|DC), setDataContainer(DC),
                              setDataValues(array)
        * clarify dataSpace(comp. ModelSpace):
            The unique spatial or temporal origin of a datapoint
            (time, coordinates, receiver/transmitter indices, counter)
            - Every inversion needs, dataValues and dataSpace
            - DataContainer contain, dataValues and dataSpace
            - initialize both with initDataSpace(), initModelSpace
        * createJacobian should also return J
    """
    def __init__(self, **kwargs):
        """Initialize.

        Attributes
        ----------
        fop : pg.frameworks.Modelling

        data : pg.DataContainer

        modelTrans : [pg.trans.TransLog()]

        Keyword Arguments
        -----------------
        **kwargs:
            fop: Modelling

        """
        self._fop = None  # pg.frameworks.Modelling .. not needed .. remove it
        self._data = None  # dataContainer
        self._modelTrans = None

        self.fop = kwargs.pop('fop', None)
        # super(Modelling, self).__init__(**kwargs)
        super().__init__(**kwargs)

        self._regionProperties = {}
        self._interRegionCouplings = []
        self._regionsNeedUpdate = False
        self._regionChanged = True
        self._regionManagerInUse = False

        self.modelTrans = pg.trans.TransLog()  # Model transformation operator


    def __hash__(self):
        """ Create a hash for Method Manager.
        """
        # ^ pg.utils.dirHash(self._regionProperties)
        if self._data is not None:
            return pg.utils.strHash(str(type(self))) ^ hash(self._data)
        else:
            return pg.utils.strHash(str(type(self)))


    def Sx(self, x):
        """ Right-hand side multiplication of Jacobian J*x.

        By default, uses self.jacobian().mult(x)
        Overwrite for efficient use with gradient-type inversion.
        """
        return self.jacobian().mult(x)


    def STy(self, y):
        """ Right-hand side multiplication of Jacobian J.T*y.

        By default, uses self.jacobian().transMult(x)
        Overwrite for efficient use with gradient-type inversion.
        """
        return self.jacobian().transMult(y)


    def __call__(self, *args, **kwargs):
        """ Call forward operator.
        """
        return self.response(*args, **kwargs)


    @property
    def fop(self):
        """ Forward operator.
        """
        return self._fop


    @fop.setter
    def fop(self, fop):
        """ Set forward operator.
        """
        if fop is not None:
            if not isinstance(fop, pg.frameworks.Modelling):
                pg.critical('Forward operator needs to be an instance of '
                            'pg.modelling.Modelling but is of type:', fop)

            self._fop = fop


    @property
    def data(self):
        """ Return data.
        """
        return self._data


    @data.setter
    def data(self, d):
        """ Set data (short property setter).
        """
        self.setData(d)


    def setData(self, data):
        """ Set data (actual version).
        """
        if isinstance(data, pg.DataContainer):
            self.setDataContainer(data)
        else:
            self._data = data


    def setDataPost(self, data):
        """ Called when the dataContainer has been set successfully.
        """
        pass


    def setDataContainer(self, data):
        """ Set Data container.
        """
        if self.fop is not None:
            pg.critical('in use?')
            self.fop.setData(data)
        else:
            super().setData(data)
            self._data = data

        self.setDataPost(self.data)


    @property
    def modelTrans(self):
        """ Return model transformation.
        """
        self._applyRegionProperties()
        if self.regionManager().haveLocalTrans():
            return self.regionManager().transModel()
        return self._modelTrans


    @modelTrans.setter
    def modelTrans(self, tm):
        """ Set model transformation.
        """
        if isinstance(tm, str):
            if tm.lower() == "log":
                tm = pg.trans.TransLog()
            elif tm.lower() == "linear" or tm.lower() == "lin":
                tm = pg.trans.Trans()
            else:  # something like "10-1000"
                raise Exception("Could not use transformation" + tm)

        self._modelTrans = tm


    def regionManager(self):
        """ Region manager.
        """
        self._regionManagerInUse = True
        # initialize RM if necessary
        super().regionManager()
        # set all local properties
        self._applyRegionProperties()
        return super().regionManager()


    @property
    def parameterCount(self):
        """ Return parameter count.
        """
        pC = self.regionManager().parameterCount()
        if pC == 0:
            pg.warn("Parameter count is 0")

        return pC


    def ensureContent(self):
        """ Whatever this is.
        """
        pass


    def initModelSpace(self, **kwargs):
        """TODO."""
        pass


    def createDefaultStartModel(self, dataVals):
        """ Create the default startmodel as the median of the data values.
        """
        pg.critical("'don't use me")


    def createStartModel(self, dataVals=None):
        """ Create the default startmodel as the median of the data values.

        Overwriting might be a good idea.
        Its used by inversion to create a valid startmodel if there are
        no starting values from the regions.
        """
        if dataVals is not None:
            mv = pg.math.median(dataVals)
            pg.info("Use median(data values)={0}".format(mv))
            sm = pg.Vector(self.parameterCount, mv)
        else:
            sm = self.regionManager().createStartModel()
        return sm


    def clearRegionProperties(self):
        """ Clear all region parameter.
        """
        self._regionChanged = True
        self._regionProperties = {}


    def regionProperties(self, regionNr=None):
        """ Return dictionary of all properties for region number regionNr.
        """
        if regionNr is None:
            return self._regionProperties

        try:
            return self._regionProperties[regionNr]
        except KeyError:
            print(self._regionProperties)
            pg.error("no region for region #:", regionNr)


    def setRegionProperties(self, regionNr, **kwargs):
        """Set region properties. regionNr can be '*' for all regions.

            startModel=None, limits=None, trans=None,
            cType=None, zWeight=None, modelControl=None,
            background=None, fix=None, single=None,
            correlationLengths=None, dip=None, strike=None

        Parameters
        ----------
        regionNr : int, [ints], '*'
            Region number, list of numbers, or wildcard "*" for all.
        startModel : float
            starting model value
        limits : [float, float]
            lower and upper limit for value using a barrier transform
        trans : str
            transformation for model barrier: "log", "cot", "lin"
        cType : int
            constraint (regularization) type
        zWeight : float
            relative weight for vertical boundaries
        background : bool
            exclude region from inversion completely (prolongation)
        fix : float
            exclude region from inversion completely (fix to value)
        single : bool
            reduce region to one unknown
        correlationLengths : [floats]
            correlation lengths for geostatistical inversion (x', y', z')
        dip : float [0]
            angle between x and x' (first correlation length)
        strike : float [0]
            angle between y and y' (second correlation length)
        """
        if regionNr == '*':
            for regionNr in self.regionManager().regionIdxs():
                self.setRegionProperties(regionNr, **kwargs)
            return
        elif isinstance(regionNr, (list, tuple)):
            for r in regionNr:
                self.setRegionProperties(r, **kwargs)
            return

        pg.verbose(f'Set property for region: {regionNr}: {kwargs}')
        if regionNr not in self._regionProperties:
            self._regionProperties[regionNr] = {'startModel': None,
                                                'modelControl': 1.0,
                                                'zWeight': 1.0,
                                                'cType': None,  # RM defaults
                                                'limits': [0, 0],
                                                'trans': 'Log',  # RM defauts
                                                'background': None,
                                                'single': None,
                                                'fix': None,
                                                'correlationLengths': None,
                                                'dip': None,
                                                'strike': None,
                                                }

        for key in list(kwargs.keys()):
            val = kwargs.pop(key)
            if val is not None:
                if self._regionProperties[regionNr][key] != val:
                    self._regionsNeedUpdate = True
                    self._regionProperties[regionNr][key] = val

        if len(kwargs) > 0:
            pg.warn('Unhandled region properties:', kwargs)

    def setInterRegionCoupling(self, region1, region2, weight=1.0):
        """Set the weighting for constraints across regions."""
        if region1 == "*":
            region1 = self.regionManager().regionIdxs()
        else:
            region1 = [region1]

        if region2 == "*":
            region2 = self.regionManager().regionIdxs()
        else:
            region2 = [region2]

        for reg1 in region1:
            for reg2 in region2:
                if reg1 != reg2 and \
                    (not self._regionProperties[reg1]['background'] and
                        not self._regionProperties[reg2]['background']):
                    self._interRegionCouplings.append([reg1, reg2, weight])

        self._regionsNeedUpdate = True

    def _applyRegionProperties(self):
        """Apply the region properties from dictionary into the region man."""
        if self._regionsNeedUpdate is False:
            return

        # call super class her because self.regionManager() calls allways
        #  __applyRegionProperies itself
        rMgr = super().regionManager()
        for rID, vals in self._regionProperties.items():
            if vals['fix'] is not None:
                if rMgr.region(rID).fixValue() != vals['fix']:
                    vals['background'] = True
                    rMgr.region(rID).setFixValue(vals['fix'])
                    self._regionChanged = True

            if vals['background'] is not None:
                if rMgr.region(rID).isBackground() != vals['background']:
                    rMgr.region(rID).setBackground(vals['background'])
                    self._regionChanged = True

            if vals['single'] is not None:
                if rMgr.region(rID).isSingle() != vals['single']:
                    rMgr.region(rID).setSingle(vals['single'])
                    self._regionChanged = True

            if vals['startModel'] is not None:
                rMgr.region(rID).setStartModel(vals['startModel'])

            if vals['trans'] is not None:
                rMgr.region(rID).setModelTransStr_(vals['trans'])

            if vals['cType'] is not None:
                if rMgr.region(rID).constraintType() != vals['cType']:
                    self.clearConstraints()
                    rMgr.region(rID).setConstraintType(vals['cType'])

            if vals['zWeight'] is not None:
                rMgr.region(rID).setZWeight(vals['zWeight'])
                self.clearConstraints()

            rMgr.region(rID).setModelControl(vals['modelControl'])

            if vals['limits'][0] != 0:
                rMgr.region(rID).setLowerBound(vals['limits'][0])

            if vals['limits'][1] > vals['limits'][0]:
                rMgr.region(rID).setUpperBound(vals['limits'][1])

            if vals['correlationLengths'] is not None:
                self.clearConstraints()
            if vals['dip'] is not None:
                self.clearConstraints()
            if vals['strike'] is not None:
                self.clearConstraints()

        for r1, r2, w in self._interRegionCouplings:
            rMgr.setInterRegionConstraint(r1, r2, w)

        self._regionsNeedUpdate = False

    def setDataSpace(self, **kwargs):
        """Set data space, e.g., DataContainer, times, coordinates."""
        if self.fop is not None:
            pg.critical('in use?')
            self.fop.setDataSpace(**kwargs)
        else:
            data = kwargs.pop('dataContainer', None)
            if isinstance(data, pg.DataContainer):
                self.setDataContainer(data)
            else:
                print(data)
                pg.critical("nothing known to do? "
                            "Implement me in derived classes")

    def estimateError(self, data, **kwargs):
        """Create data error fallback when the data error is not known.

        Should be implemented method-specific.
        """
        raise Exception("Needed?? Implement me in derived classes")
        # data = data * (pg.randn(len(data)) * errPerc / 100. + 1.)
        # return data

    def drawModel(self, ax, model, **kwargs):
        """Draw a model into a given axis."""
        if self.fop is not None:
            pg.critical('in use?')
            self.fop.drawModel(ax, model, **kwargs)
        else:
            print(kwargs)
            raise Exception("No yet implemented")

    def drawData(self, ax, data, **kwargs):
        """Draw data into a given axis."""
        if self.fop is not None:
            self.fop.drawData(ax, data, **kwargs)
        else:
            print(kwargs)
            raise Exception("No yet implemented")


class LinearModelling(Modelling):
    """Modelling class for linearized problems with a given matrix."""

    def __init__(self, A):
        """Initialize by storing the (reference to the) matrix."""
        super().__init__()
        self.A = A
        self.setJacobian(self.A)

    def response(self, model):
        """Linearized forward modelling by matrix-vector product."""
        return self.A * model

    def createJacobian(self, model):
        """Do not compute a jacobian (linear)."""
        pass

    @property
    def parameterCount(self):
        """Define the number of parameters from the matrix size."""
        return self.A.cols()


class Block1DModelling(Modelling):
    """General forward operator for 1D layered models.

    Model space: [thickness_i, parameter_jk],
    with i = 0 - nLayers-1, j = (0 .. nLayers), k=(0 .. nPara)
    """

    def __init__(self, nPara=1, nLayers=4, **kwargs):
        """Constructor.

        Parameters
        ----------
        nLayers : int [4]
            Number of layers.
        nPara : int [1]
            Number of parameters per layer
            (e.g. nPara=2 for resistivity and phase)
        """
        self._nLayers = 0
        super(Block1DModelling, self).__init__(**kwargs)
        self._withMultiThread = True
        self._nPara = nPara  # number of parameters per layer

        self.initModelSpace(nLayers)

    @property
    def nPara(self):
        """Number of parameters."""
        return self._nPara

    @property
    def nLayers(self):
        """Number of layers."""
        return self._nLayers

    @nLayers.setter
    def nLayers(self, nLayers):
        """Set number of layers."""
        return self.initModelSpace(nLayers)

    def initModelSpace(self, nLayers):
        """Set number of layers for the 1D block model."""
        if nLayers == self._nLayers:
            return
        self._nLayers = nLayers

        if nLayers < 2:
            pg.critical("Number of layers need to be at least 2")

        mesh = pg.meshtools.createMesh1DBlock(nLayers, self._nPara)
        self.clearRegionProperties()
        self.setMesh(mesh)
        # setting region 0 (layers) and 1..nPara (values)
        for i in range(1 + self._nPara):
            self.setRegionProperties(i, trans='log')

        if self._withMultiThread:
            self.setMultiThreadJacobian(2*nLayers - 1)

        # self._applyRegionProperties()

    def drawModel(self, ax, model, **kwargs):
        """Draw model into a given axis."""
        pg.viewer.mpl.drawModel1D(ax=ax,
                                  model=model,
                                  plot='loglog',
                                  xlabel=kwargs.pop('xlabel',
                                                    'Model parameter'),
                                  **kwargs)
        return ax

    def drawData(self, ax, data, err=None, label=None, **kwargs):
        """Default data view.

        Modelling creates the data and should know best how to draw them.

        Probably ugly and you should overwrite it in your derived forward
        operator.
        """
        nData = len(data)
        yVals = range(1, nData+1)
        ax.loglog(data, yVals,
                  label=label,
                  **DEFAULT_STYLES.get(label, DEFAULT_STYLES['Default'])
                  )

        if err is not None:
            ax.errorbar(data, yVals,
                        xerr=err*data,
                        label='Error',
                        **DEFAULT_STYLES.get('Error',
                                             DEFAULT_STYLES['Default'])
                        )

        ax.set_ylim(max(yVals), min(yVals))
        ax.set_xlabel('Data')
        ax.set_ylabel('Data Number')
        return ax


class MeshModelling(Modelling):
    """Modelling class with a mesh discretization."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._axs = None
        self._meshNeedsUpdate = True
        self._baseMesh = None               ## keep a copy of the original mesh
        # optional p2 refinement for forward task
        self._refineP2 = False
        self._refineH2 = True
        self._pd = None
        self._C = None # custom Constraints matrix


    def __hash__(self):
        """ Unique hash for caching.
        """
        return super().__hash__() ^ hash(self.mesh())


    @property
    def paraDomain(self):
        """ Return parameter (inverse) mesh.
        """
        # We need our own copy here because its possible that we want to use
        # the mesh after the fop was deleted
        if not self.mesh():
            pg.critical('paraDomain needs a mesh')

        self._pd = pg.Mesh(self.regionManager().paraDomain())
        return self._pd


    def setCustomConstraints(self, C):
        """ Set custom constraints matrix for lazy evaluation.

        To remove them set it to 'None' again.
        """
        self._C = C


    def createConstraints(self):
        """ Create constraint matrix.
        """
        # just ensure there is valid mesh
        # self.mesh()

        foundGeoStat = False
        for reg, props in self.regionProperties().items():
            if not props['background'] and \
                props['correlationLengths'] is not None or \
                    props['dip'] is not None or props['strike'] is not None:

                cL = props.get('correlationLengths') or 5
                dip = props.get('dip') or 0
                strike = props.get('strike') or 0

                pg.info('Creating GeostatisticConstraintsMatrix for region' +
                        f' {reg} with: I={cL}, dip={dip}, strike={strike}')

                if foundGeoStat is True:
                    pg.critical('Only one global GeostatisticConstraintsMatrix'
                                'possible at the moment.')

                # keep a copy of C until refcounting in the core works
                self._C = pg.matrix.GeostatisticConstraintsMatrix(
                    mesh=self.paraDomain, I=cL, dip=dip, strike=strike,
                )
                foundGeoStat = True
                self.setConstraints(self._C)

        if foundGeoStat is False:
            super().createConstraints()

        return self.constraints()


    def paraModel(self, model):
        """ Return parameter model, i.e. model mapped back with cell markers.
        """
        mod = model[self.paraDomain.cellMarkers()]
        if isinstance(mod, np.ndarray):
            mod = pg.Vector(mod)
            # Else next line fails as np.array does not allow set attributes.

        mod.isParaModel = True
        return mod


    def ensureContent(self):
        """ Internal function to ensure there is a valid initialized mesh.

        Initialization means the cell marker are recounted and/or there was a
        mesh refinement or boundary enlargement, all to fit the needs for the
        method-depending forward problem.
        """
        # Need to call this once to be sure the mesh is initialized when needed
        self.mesh()


    def setMeshPost(self, mesh):
        """ Interface to be called when the mesh has been set successfully.

        Might be overwritten by child classes.
        """
        pass


    def createRefinedFwdMesh(self, mesh):
        """ Refine the current mesh for higher accuracy.

        This is called automatic when accessing self.mesh() so it ensures any
        effect of changing region properties (background, single).
        """
        m = pg.Mesh(mesh)
        if self._refineH2:
            pg.info("Creating refined mesh (H2) to solve forward task.")
            m = m.createH2()

        if self._refineP2:
            pg.info("Creating refined mesh (P2) to solve forward task.")
            m = m.createP2()

        pg.info("Mesh for forward task:", m)
        return m


    def createFwdMesh_(self):
        """ Create forward mesh.
        """
        pg.info("Creating forward mesh from region infos.")
        m = pg.Mesh(self.regionManager().mesh())

        regionIds = self.regionManager().regionIdxs()
        for iId in regionIds:
            pg.verbose("\tRegion: {0}, Parameter: {1}, PD: {2},"
                       " Single: {3}, Background: {4}, Fixed: {5}".format(
                           iId,
                           self.regionManager().region(iId).parameterCount(),
                           self.regionManager().region(iId).isInParaDomain(),
                           self.regionManager().region(iId).isSingle(),
                           self.regionManager().region(iId).isBackground(),
                           self.regionManager().region(iId).fixValue()))

        m = self.createRefinedFwdMesh(m)
        self.setMeshPost(m)

        self._regionChanged = False
        super().setMesh(m, ignoreRegionManager=True)

        if self._C is not None:
            pg.info('Set custom constraints matrix.')
            self.setConstraints(self._C)


    def mesh(self):
        """ Returns the currently used mesh.
        """
        self._applyRegionProperties()

        if self._regionManagerInUse and self._regionChanged is True:
            self.createFwdMesh_()

        return super().mesh()


    def setMesh(self, mesh, ignoreRegionManager=False):
        """ Set mesh and specify whether region manager can be ignored.
        """
        # pg._b('setMesh', id(mesh), mesh, ignoreRegionManager)
        # keep a copy, just in case
        self._baseMesh = mesh

        if ignoreRegionManager is False:
            self._regionManagerInUse = True

        # Modelling without region manager
        if ignoreRegionManager is True or not self._regionManagerInUse:
            self._regionManagerInUse = False
            super(Modelling, self).setMesh(mesh, ignoreRegionManager=True)

            self.setMeshPost(mesh)
            return

        # copy the mesh to the region manager who renumber cell markers
        self.clearRegionProperties()
        self.regionManager().setMesh(mesh)
        self.setDefaultBackground()


    def setDefaultBackground(self):
        """ Set the lowest region to background if several exist.
        """
        regionIds = self.regionManager().regionIdxs()
        pg.info("Found {} regions.".format(len(regionIds)))
        if len(regionIds) > 1:
            bk = pg.sort(regionIds)[0]
            pg.info("Region with smallest marker set to background "
                    "(marker={0})".format(bk))
            self.setRegionProperties(bk, background=True)


    def drawModel(self, ax, model, **kwargs):
        """Draw the model as mesh-based distribution."""
        mod = None
        # TODO needs to be checked if mapping is always ok (region example)
        # is (len(model) == self.paraDomain.cellCount() or \
        if hasattr(model, "isParaModel") and model.isParaModel is False:
            # pg._y(model.isParaModel)
            mod = self.paraModel(model)
        elif hasattr(model, "isParaModel") and model.isParaModel is True:
            # pg._g(model.isParaModel)
            mod = model
        elif len(model) == self.paraDomain.nodeCount():
            # why nodeCount? a field as model result?
            # pg._b('node count')
            mod = model
        elif len(model) == self.paraDomain.cellCount():
            # pg._b('cell count')
            mod = model
        else:
            mod = self.paraModel(model)

        if ax is None:
            if self._axs is None:
                self._axs, _ = pg.show()
            ax = self._axs

        if hasattr(ax, '__cBar__'):
            # we assume the axes already holds a valid mappable and we only
            # update the model data
            cBar = ax.__cBar__
            kwargs.pop('label', None)
            kwargs.pop('cMap', None)
            pg.viewer.mpl.setMappableData(cBar.mappable, mod, **kwargs)
        else:
            diam = kwargs.pop('diam', None)

            ax, cBar = pg.show(mesh=self.paraDomain,
                               data=mod,
                               label=kwargs.pop('label', 'Model parameter'),
                               logScale=kwargs.pop('logScale', False),
                               ax=ax,
                               **kwargs
                               )
            if diam is not None:
                pg.viewer.mpl.drawSensors(ax, self.data.sensors(), diam=diam,
                                          edgecolor='black', facecolor='white')

        return ax, cBar


class PetroModelling(MeshModelling):
    """Combine petrophysical relation with the modelling class f(p).

    Combine petrophysical relation :math:`p(m)` with a modelling class
    :math:`f(p)` to invert for the petrophysical model :math:`p` instead
    of the geophysical model :math:`m`.

    :math:`p` be the petrophysical model, e.g., porosity, saturation, ...
    :math:`m` be the geophysical model, e.g., slowness, resistivity, ...
    """

    def __init__(self, fop, petro, **kwargs):
        """Save forward class and transformation, create Jacobian matrix."""
        self._f = fop
        # self._f createStartModel might be called and depends on the regionMgr
        self._f.regionManager = self.regionManager
        # self.createRefinedFwdMesh depends on refinement strategy of self._f
        self.createRefinedFwdMesh = self._f.createRefinedFwdMesh

        super(PetroModelling, self).__init__(fop=None, **kwargs)
        # petroTrans.fwd(): p(m), petroTrans.inv(): m(p)
        self._petroTrans = petro  # class defining p(m)

        self._jac = pg.matrix.MultRightMatrix(self._f.jacobian())
        self.setJacobian(self._jac)

    @property
    def petro(self):
        """Petrophysical model transformation."""
        return self._petroTrans

    def setMeshPost(self, mesh):
        """Set mesh after init."""
        self._f.setMesh(mesh, ignoreRegionManager=True)

    def setDataPost(self, data):
        """Set data after init."""
        self._f.setData(data)

    def createStartModel(self, data):
        """Use inverse transformation to get m(p) for the starting model."""
        sm = self._f.createStartModel(data)
        pModel = self._petroTrans.inv(sm)
        return pModel

    def response(self, model):
        """Use transformation to get p(m) and compute response f(p)."""
        tModel = self._petroTrans.fwd(model)
        ret = self._f.response(tModel)
        return ret

    def createJacobian(self, model):
        r"""Fill the individual jacobian matrices.

        J = dF(m) / dm = dF(m) / dp  * dp / dm
        """
        tModel = self._petroTrans.fwd(model)

        self._f.createJacobian(tModel)
        self._jac.A = self._f.jacobian()
        self._jac.r = self._petroTrans.deriv(model)  # set inner derivative
        # print(self._jac.A.rows(), self._jac.A.cols())
        # print(self._jac.r)
        # pg._r("create Jacobian", self, self._jac)
        self.setJacobian(self._jac)  # to be sure .. test if necessary


# 220817 to be changed later!!
# class JointModelling(Modelling):
class JointModelling(MeshModelling):
    """Cumulative (joint) forward operator."""

    def __init__(self, fopList):
        """Initialize with lists of forward operators."""
        super().__init__()
        self.fops = fopList
        self.jac = pg.matrix.BlockMatrix()

        # self.modelTrans = self.fops[0].modelTrans
        self.modelTrans = pg.trans.TransLogLU()
        self.fops[0].regionManager()
        self.setRegionManager(self.fops[0].regionManagerRef())

    def createStartModel(self, data):
        """Use inverse transformation to get m(p) for the starting model."""
        sm = self._f.createStartModel(data)
        pModel = self._petroTrans.inv(sm)
        return pModel

    def response(self, model):
        """Concatenate responses for all fops."""
        resp = []
        for f in self.fops:
            resp.extend(f.response(model))
        return resp

    def createJacobian(self, model):
        """Fill the individual Jacobian matrices."""
        self.initJacobian()
        for f in self.fops:
            f.createJacobian(model)

    def setData(self, data):
        """Distribute list of data to the forward operators."""
        if len(data) != len(self.fops):
            pg.critical("Please provide data for all forward operators")

        self._data = data
        nData = 0
        for i, fi in enumerate(self.fops):
            fi.setData(data[i])
            self.jac.addMatrix(fi.jacobian(), nData, 0)
            nData += data[i].size()  # update total vector length
        self.setJacobian(self.jac)

    def setMesh(self, mesh, **kwargs):  # to be removed from here
        """Set the parameter mesh to all fops."""
        for fi in self.fops:
            fi.setMesh(mesh)

# 220817 to be implemented!!
# class JointMeshModelling(JointModelling):
#    def __init__(self, fopList):
        # super().__init__(self, fopList)
        # self.setRegionManager(self.fops[0].regionManagerRef())


class LCModelling(Modelling):
    """2D Laterally constrained (LC) modelling.

    2D Laterally constrained (LC) modelling based on BlockMatrices.
    """

    def __init__(self, fop, **kwargs):
        """Parameters: fop class ."""
        super(LCModelling, self).__init__()
        self._singleRegion = False
        self._fopTemplate = fop
        self._fopKwargs = kwargs
        self._fops1D = []
        self._mesh = None
        self._nSoundings = 0
        self._parPerSounding = 0
        self._jac = None

        self.soundingPos = None

    def setDataBasis(self, **kwargs):
        """Set homogeneous data basis.

        Set a common data basis to all forward operators.
        If you want individual you need to set them manually.
        """
        for f in self._fops1D:
            f.setDataBasis(**kwargs)

    def initModelSpace(self, nLayers):
        """Initialize model space."""
        for i, f in enumerate(self._fops1D):
            f.initModelSpace(nLayers)

    def createDefaultStartModel(self, models):
        """Create default starting model."""
        sm = pg.Vector()
        for i, f in enumerate(self._fops1D):
            sm = pg.cat(sm, f.createDefaultStartModel(models[i]))
        return sm

    def response(self, par):
        """Cut together forward responses of all soundings."""
        mods = np.asarray(par).reshape(self._nSoundings, self._parPerSounding)

        resp = pg.Vector(0)
        for i in range(self._nSoundings):
            r = self._fops1D[i].response(mods[i])
            # print("i:", i, mods[i], r)
            resp = pg.cat(resp, r)

        return resp

    def createJacobian(self, par):
        """Create Jacobian matrix by creating individual Jacobians."""
        mods = np.asarray(par).reshape(self._nSoundings, self._parPerSounding)

        for i in range(self._nSoundings):
            self._fops1D[i].createJacobian(mods[i])

    def createParametrization(self, nSoundings, nLayers=4, nPar=1):
        """Create LCI mesh and suitable constraints informations.

        Parameters
        ----------
        nLayers : int
            Numbers of depth layers

        nSoundings : int
            Numbers of 1D measurements to laterally constrain

        nPar : int
            Numbers of independent parameter types,
            e.g., nPar = 1 for VES (invert for resisitivies),
            nPar = 2 for VESC (invert for resisitivies and phases)
        """
        nCols = (nPar+1) * nLayers - 1  # fail for VES-C
        self._parPerSounding = nCols
        self._nSoundings = nSoundings

        self._mesh = pg.meshtools.createMesh2D(range(nCols + 1),
                                               range(nSoundings + 1))
        self._mesh.rotate(pg.RVector3(0, 0, -np.pi/2))

        cm = np.ones(nCols * nSoundings) * 1

        if not self._singleRegion:
            for i in range(nSoundings):
                for j in range(nPar):
                    cm[i * self._parPerSounding + (j+1) * nLayers-1:
                       i * self._parPerSounding + (j+2) * nLayers-1] += (j+1)

        self._mesh.setCellMarkers(cm)
        self.setMesh(self._mesh)

        pID = self.regionManager().paraDomain().cellMarkers()
        cID = [c.id() for c in self._mesh.cells()]
        # print(np.array(pID))
        # print(np.array(cID))
        # print(self.parameterCount
        perm = [0]*self.parameterCount
        for i in range(len(perm)):
            perm[pID[i]] = cID[i]

        # print(perm)
        self.regionManager().permuteParameterMarker(perm)
        # print(self.regionManager().paraDomain().cellMarkers())

    def initJacobian(self, dataVals, nLayers, nPar=None):
        """Initialize Jacobian matrix.

        Parameters
        ----------
        dataVals : ndarray | RMatrix | list
            Data values of size (nSounding x Data per sounding).
            All data per sounding need to be equal in length.
            If they don't fit into a matrix use list of sounding data.
        """
        nSoundings = len(dataVals)

        if nPar is None:
            # TODO get nPar Infos from fop._fopTemplate
            nPar = 1

        self.createParametrization(nSoundings, nLayers=nLayers, nPar=nPar)

        if self._jac is not None:
            self._jac.clear()
        else:
            self._jac = pg.matrix.BlockMatrix()

        self.fops1D = []
        nData = 0

        for i in range(nSoundings):
            kwargs = {}
            for key, val in self._fopKwargs.items():
                if hasattr(val, '__iter__'):
                    if len(val) == nSoundings:
                        kwargs[key] = val[i]
                else:
                    kwargs[key] = val

            f = None
            if issubclass(self._fopTemplate, pg.frameworks.Modelling):
                f = self._fopTemplate(**kwargs)
            else:
                f = type(self._fopTemplate)(self.verbose, **kwargs)

            f.setMultiThreadJacobian(self._parPerSounding)

            self._fops1D.append(f)

            nID = self._jac.addMatrix(f.jacobian())
            self._jac.addMatrixEntry(nID, nData, self._parPerSounding * i)
            nData += len(dataVals[i])

        self._jac.recalcMatrixSize()
        # print("Jacobian size:", self.J.rows(), self.J.cols(), nData)
        self.setJacobian(self._jac)

    def drawModel(self, ax, model, **kwargs):
        """Draw models as stitched 1D model section."""
        mods = np.asarray(model).reshape(self._nSoundings,
                                         self._parPerSounding)
        pg.viewer.mpl.showStitchedModels(mods, ax=ax, useMesh=True,
                                         x=self.soundingPos,
                                         **kwargs)


class ParameterModelling(Modelling):
    """Model with symbolic parameter names instead of numbers."""

    def __init__(self, func=None, **kwargs):
        """Initialize, optionally with given function."""
        self.function = None
        self._params = {}
        self.dataSpace = None  # x, t freqs, or whatever
        self.defaultModelTrans = 'lin'
        self._model = None # last active to use the FOP as function for interpolation

        super(ParameterModelling, self).__init__(**kwargs)

        if func is not None:
            self._initFunction(func)

    @property
    def params(self):
        """Return number of parameters."""
        return self._params

    @property
    def model(self):
        return self._model

    def _initFunction(self, funct):
        """Init any function and interpret possible args and kwargs."""
        self.function = funct
        # the first varname is supposed to be f or freqs
        self.dataSpaceName = funct.__code__.co_varnames[0]
        pg.debug('data space:', self.dataSpaceName)

        args = funct.__code__.co_varnames[1:funct.__code__.co_argcount]
        for varname in args:
            if varname != 'verbose':
                pg.debug('add parameter:', varname)
                self._params[varname] = 0.0

        # nPara = len(self._params.keys())  # not used!

        for i, [k, p] in enumerate(self._params.items()):
            self.addParameter(k, id=i, cType=0,
                              single=True,
                              trans=self.defaultModelTrans,
                              startModel=1)

    def response(self, params):
        """Compute and return model response."""
        if np.isnan([*params]).any():
            print(params)
            pg.critical('invalid params for response')
        if self.dataSpace is None:
            pg.critical('no data space given')

        ## Store model, which can be handy
        self._model = np.array(params)
        ret = self.function(self.dataSpace, *params)
        return ret

    def __call__(self, dataSpace):
        return self.function(np.asarray(dataSpace), *self._model.tolist())

    def setRegionProperties(self, k, **kwargs):
        """Set Region Properties by parameter name."""
        if isinstance(k, int) or (k == '*'):
            super(ParameterModelling, self).setRegionProperties(k, **kwargs)
        else:
            self.setRegionProperties(self._params[k], **kwargs)

    def addParameter(self, name, id=None, **kwargs):
        """Add a parameter."""
        if id is None:
            id = len(self._params)
        self._params[name] = id
        self.regionManager().addRegion(id)
        self.setRegionProperties(name, **kwargs)
        return id

    def drawModel(self, ax=None, model=None):
        """"""
        if model is None:
            model = self._model

        label = ''
        for k, p in self._params.items():
            label += k + "={0} ".format(pg.utils.prettyFloat(model[p]))
        pg.info("Model: ", label)


class PriorModelling(MeshModelling):
    """Forward operator for grabbing values out of a mesh (prior data)."""

    def __init__(self, mesh=None, pos=None, **kwargs):
        """Init with mesh and some positions that are converted into ids."""
        super().__init__(**kwargs)
        self.pos = pos
        if mesh is not None:
            self.setMesh(mesh)

    def setMesh(self, mesh):
        """Set mesh, save index vector and compute Jacobian."""
        super().setMesh(mesh)
        self.ind = np.zeros(len(self.pos), dtype=np.int32)
        for i, po in enumerate(self.pos):
            cell = mesh.findCell(po)
            if cell is None:
                raise IndexError(f"Could not find cell at position {po}!")
            else:
                self.ind[i] = cell.id()

        # self.ind = np.array([mesh.findCell(po).id() for po in self.pos])
        self.J = pg.SparseMapMatrix()
        self.J.resize(len(self.ind), mesh.cellCount())
        for i, ii in enumerate(self.ind):
            self.J.setVal(i, ii, 1.0)

        self.setJacobian(self.J)

    def response(self, model):
        """Return values at the indexed cells."""
        return model[self.ind]

    def createJacobian(self, model):
        """Do nothing (linear)."""
        pass

    def createRefinedFwdMesh(self, mesh):
        """Create refined forward mesh: do nothing here to prevent this."""
        return mesh
