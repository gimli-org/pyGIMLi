"""ERT modelling operator classes.

* abstract base class providing a plotting function
* standard BERT modelling class using pygimli.core (C++) functionality
* 2.5D non-optimized totalfield forward operator for ERT (reference)
"""

import numpy as np

import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli import pf
from .visualization import showERTData


class ERTModellingBase(MeshModelling):
    """Modelling base class for ERT modelling."""

    def __init__(self, **kwargs):  # actually useless def
        super().__init__(**kwargs)

    def drawData(self, ax, data=None, **kwargs):
        """Draw data on given axes.

        Args
        ----
        ax: MPL Axe
            Axes instance to draw into.
        data: iterable | :py:func:pygimli.pyhsics.ert.DataContainer`
            Datavalues to draw

            * data is array, taken as values and draw on self.data structure
            * data is Datacontainer, draw data['rhoa'] on data structure
            * data is None, draw self.data['rhoa'] on self.data structure

        Keyword Args
        ------------
        vals: iterable
            Draw vals on self.data structure if no data are given.

        Remaining kwargs are forwarded to
        :py:func:pygimli.pyhsics.ert.showERTData`
        """
        kwargs['label'] = kwargs.pop('label', pg.unit('res'))
        kwargs['cMap'] = kwargs.pop('cMap', pg.utils.cMap('res'))

        vals = None
        if hasattr(data, '__iter__'):
            # show data values from array
            vals = data
            data = self.data
        elif isinstance(data, pg.DataContainerERT):
            # data is given DataContainer
            pass
        else:
            # data is self.DataContainer
            data = self.data

        if vals is None:
            vals = kwargs.pop('vals', data['rhoa'])

        return showERTData(data, vals=vals, ax=ax, **kwargs)

    def drawModel(self, ax, model, **kwargs):
        """Draw the para domain with option model values."""
        kwargs.setdefault('label', pg.unit('res'))
        kwargs.setdefault('cMap', pg.utils.cMap('res'))
        kwargs.setdefault('logScale', True)

        return super().drawModel(ax=ax, model=model, **kwargs)


class ERTModelling(ERTModellingBase):
    """Forward operator for Electrical Resistivity Tomography.

    Note
    ----
    Convention for complex resistiviy inversion:
    We want to use logarithm transformation for the imaginary part of model
    so we need the startmodel to have positive imaginary parts.
    The sign is flipped back to physical correct assumption before we call
    the response function.
    The Jacobian is calculated with negative imaginary parts and will
    be a conjugated complex block matrix for further calulations.
    """

    def __init__(self, sr=True, verbose=False):
        super().__init__()

        # don't use DC*fop or its regionmanager directly
        #
        self._core = None
        if sr:
            self._core = pg.core.DCSRMultiElectrodeModelling(verbose=verbose)
        else:
            self._core = pg.core.DCMultiElectrodeModelling(verbose=verbose)

        # Its good that the core knows about the actual RM
        self._core.setRegionManager(self.regionManagerRef())

        self._core.initJacobian()
        self.setJacobian(self._core.jacobian())

        # called from the ERTManager .. needed?
        self.solution = self._core.solution
        self.setComplex = self._core.setComplex
        self.complex = self._core.complex
        self.calculate = self._core.calculate
        self.calcGeometricFactor = self._core.calcGeometricFactor
        self.mapERTModel = self._core.mapERTModel

        self._conjImag = False  # the imaginary parts are flipped for log trans


    def setVerbose(self, v):
        """ Set verbosity.
        """
        super().setVerbose(v)
        self._core.setVerbose(v)


    def setDefaultBackground(self):
        """ Set the default background behavior.
        """
        # if self.complex(): # deactivated, do it by hand
        #     self.regionManager().addRegion(3, self._baseMesh, 2)
        super().setDefaultBackground()
        # regionIds = self.regionManager().regionIdxs()
        # pg.info("Found {} regions.".format(len(regionIds)))

        # if len(regionIds) > 1:
        #     bk = pg.sort(regionIds)[0]
        #     pg.info(f"(ERTModelling) Region with smallest marker ({bk}) "
        #             "set to background.")
        #     self.setRegionProperties(bk, background=True)


    @property
    def parameterCount(self):
        """ Return number of parameters.
        """
        return self.regionManager().parameterCount() * (1 + self.complex())


    def createConstraints(self, C=None):
        """ Create constraint matrix (special type for this).
        """
        super().createConstraints()  # standard C
        if self.complex():
            if C is not None:
                self.C1 = C
            elif isinstance(self.constraints(), pg.SparseMapMatrix):
                self.C1 = pg.SparseMapMatrix(self.constraintsRef())
                # make a copy because it will be overwritten
            else:
                self.C1 = self.constraints()

            self.C = pg.matrix.RepeatDMatrix(self.C1, 2)
            self.setConstraints(self.C)


    def createStartModel(self, dataVals):
        """ Create Starting model for ERT inversion.
        """
        if self.complex():
            dataC = pg.utils.toComplex(dataVals)
            nModel = self.parameterCount // 2
            smRe = np.ones(nModel) * np.median(np.median(dataC.real))
            smIm = np.ones(nModel) * np.median(np.median(dataC.imag))

            if min(smIm) < 0:
                # we want positive phase model
                sm = smRe - 1j * smIm
                pg.info("Model imaginary part being flipped to positive.")
                self._conjImag = True
            else:
                sm = smRe + 1j * smIm

            return pg.utils.squeezeComplex(sm)  # complex impedance
        else:
            return super().createStartModel(dataVals)

    def flipImagPart(self, v):
        """Flip imaginary port (convention)."""
        z = pg.utils.toComplex(v)
        pg.warn('pre min/max={0} / {1} im: {2} / {3}'.format(
            pf(min(z.real)), pf(max(z.real)),
            pf(min(z.imag)), pf(max(z.imag))))

        v = pg.utils.squeezeComplex(pg.utils.toComplex(v),
                                    conj=self._conjImag)

        z = pg.utils.toComplex(v)
        pg.warn('pos min/max={0} / {1} im: {2} / {3}'.format(
            pf(min(z.real)), pf(max(z.real)),
            pf(min(z.imag)), pf(max(z.imag))))
        return v

    def response(self, mod):
        """Forward response (apparent resistivity)."""
        # ensure the mesh is initialized
        self.mesh()

        if self.complex() and self._conjImag:
            pg.warn('flip imaginary part for response calc')
            mod = self.flipImagPart(mod)

        resp = self._core.response(mod)

        if self.complex() and self._conjImag:
            pg.warn('backflip imaginary part after response calc')
            resp = self.flipImagPart(resp)

        return resp

    def createJacobian(self, mod):
        """Compute Jacobian matrix and store but not return."""
        # ensure the mesh is initialized
        self.mesh()
        if self.complex():
            if self._conjImag:
                pg.warn("Flipping imaginary part for jacobian calc")
                mod = self.flipImagPart(mod)

            self._core.createJacobian(mod)
            self._J = pg.utils.squeezeComplex(self._core.jacobian(),
                                              conj=self._conjImag
                                              )
            self.setJacobian(self._J)
            # pg._r("create Jacobian", self, self._J)
            return self._J

        return self._core.createJacobian(mod)

    def setDataPost(self, data):
        """Set data (at a later stage)."""
        self._core.setData(data)

    def setMeshPost(self, mesh):
        """Set mesh (at a later stage)."""
        self._core.setMesh(mesh, ignoreRegionManager=True)


class ERTModellingReference(ERTModellingBase):
    """Reference implementation for 2.5D Electrical Resistivity Tomography."""

    def __init__(self, **kwargs):
        super().__init__()

        self.subPotentials = None
        self.lastResponse = None

        # only for mixed boundary hack since this need to know resistivies.
        self.resistivity = None

        # abscissa k and weight for 2.5 inverse cos-transform
        self.k = None
        self.w = None

    def response(self, model):
        """Solve forward task and return apparent resistivity for self.mesh."""
        # NOTE TODO can't be MT until mixed boundary condition depends on
        # self.resistivity
        pg.tic()
        if not self.data.allNonZero('k'):
            pg.error('Need valid geometric factors: "k".')
            pg.warn('Fallback "k" values to -sign("rhoa")')
            self.data.set('k', -pg.math.sign(self.data['rhoa']))

        mesh = self.mesh()

        nDof = mesh.nodeCount()
        elecs = self.data.sensorPositions()

        nEle = len(elecs)
        nData = self.data.size()

        self.resistivity = res = self.createMappedModel(model, -1.0)

        if self.verbose:
            print("Calculate response for model:", min(res), max(res))

        rMin = elecs[0].dist(elecs[1]) / 2.0
        rMax = elecs[0].dist(elecs[-1]) * 2.0

        k, w = self.getIntegrationWeights(rMin, rMax)

        self.k = k
        self.w = w

        # pg.show(mesh, res, label='res')
        # pg.wait()

        rhs = self.createRHS(mesh, elecs)

        # store all potential fields
        u = np.zeros((nEle, nDof))
        self.subPotentials = [pg.Matrix(nEle, nDof) for i in range(len(k))]

        for i, ki in enumerate(k):
            uE = pg.solve(mesh, a=1./res, b=-(ki * ki)/res, f=rhs,
                          bc={'Robin': ['*', self.mixedBC]},
                          userData={'sourcePos': elecs, 'k': ki},
                          verbose=False, stats=0, debug=False)
            self.subPotentials[i] = uE
            u += w[i] * uE

        # collect potential matrix,
        # i.e., potential for all electrodes and all injections
        pM = np.zeros((nEle, nEle))

        for i in range(nEle):
            pM[i] = pg.interpolate(mesh, u[i, :], destPos=elecs)

        # collect resistivity values for all 4 pole measurements
        r = np.zeros(nData)

        for i in range(nData):
            iA = int(self.data['a'][i])
            iB = int(self.data['b'][i])
            iM = int(self.data['m'][i])
            iN = int(self.data['n'][i])

            uAB = pM[iA] - pM[iB]
            r[i] = uAB[iM] - uAB[iN]

        self.lastResponse = r * self.data['k']

        if self.verbose:
            print("Resp min/max: {0} {1} {2}s".format(min(self.lastResponse),
                                                      max(self.lastResponse),
                                                      pg.dur()))

        return self.lastResponse

    def createJacobian(self, model):
        """Create Jacobian matrix for model and store it in self.jacobian()."""
        if self.subPotentials is None:
            self.response(model)

        J = self.jacobian()
        J.resize(self.data.size(), self.parameterCount)

        cells = self.mesh().findCellByMarker(0, -1)
        Si = pg.matrix.ElementMatrix()
        St = pg.matrix.ElementMatrix()

        u = self.subPotentials

        pg.tic()
        if self.verbose:
            print("Calculate sensitivity matrix for model: ",
                  min(model), max(model))

        Jt = pg.Matrix(self.data.size(), self.parameterCount)

        for kIdx, w in enumerate(self.w):
            k = self.k[kIdx]
            w = self.w[kIdx]

            Jt *= 0.
            A = pg.matrix.ElementMatrixMap()

            for i, c in enumerate(cells):
                modelIdx = c.marker()

                # 2.5D
                Si.u2(c)
                Si *= k * k
                Si += St.ux2uy2uz2(c)

                # 3D
                # Si.ux2uy2uz2(c); w = w* 2

                A.add(modelIdx, Si)

            for dataIdx in range(self.data.size()):

                a = int(self.data['a'][dataIdx])
                b = int(self.data['b'][dataIdx])
                m = int(self.data['m'][dataIdx])
                n = int(self.data['n'][dataIdx])
                Jt[dataIdx] = A.mult(u[kIdx][a] - u[kIdx][b],
                                     u[kIdx][m] - u[kIdx][n])

            J += w * Jt

        m2 = model*model
        k = self.data['k']

        for i in range(J.rows()):
            J[i] /= (m2 / k[i])

        if self.verbose:
            sumsens = np.zeros(J.rows())
            for i in range(J.rows()):
                sumsens[i] = pg.sum(J[i])
            print("sens sum: median = ", pg.math.median(sumsens),
                  " min = ", pg.min(sumsens),
                  " max = ", pg.max(sumsens))

    def calcGeometricFactor(self, data):
        """Calculate geometry factors for a given dataset."""
        if pg.y(data.sensorPositions()) == pg.z(data.sensorPositions()):
            k = np.zeros(data.size())
            for i in range(data.size()):
                a = data.sensorPosition(data['a'][i])
                b = data.sensorPosition(data['b'][i])
                m = data.sensorPosition(data['m'][i])
                n = data.sensorPosition(data['n'][i])
                k[i] = 1./(2.*np.pi) * (1./a.dist(m) - 1./a.dist(n) -
                                        1./b.dist(m) + 1./b.dist(n))
            return k
        else:
            raise BaseException("Please use BERT for non-standard "
                                "data sets" + str(data))

    def uAnalytical(self, p, sourcePos, k):
        """
        Calculate analytical potential for homogeneous halfspace.

            For sigma = 1 [S m]
        """
        r1A = (p - sourcePos).abs()
        # Mirror on surface at depth=0
        r2A = (p - pg.RVector3(1.0, -1.0, 1.0) * sourcePos).abs()

        if r1A > 1e-12 and r2A > 1e-12:
            return (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k)) / \
                (2.0 * np.pi)
        else:
            return 0.

    def getIntegrationWeights(self, rMin, rMax):
        """TODO WRITEME."""
        nGauLegendre = max(int((6.0 * np.log10(rMax / rMin))), 4)
        nGauLaguerre = 4

        k = pg.Vector()
        w = pg.Vector()

        k0 = 1.0 / (2.0 * rMin)
        pg.GaussLegendre(0.0, 1.0, nGauLegendre, k, w)
        kLeg = k0 * k * k
        wLeg = 2.0 * k0 * k * w / np.pi

        pg.GaussLaguerre(nGauLaguerre, k, w)
        kLag = k0 * (k + 1.0)
        wLag = k0 * np.exp(k) * w / np.pi

        return pg.cat(kLeg, kLag), pg.cat(wLeg, wLag)

    def mixedBC(self, boundary, userData):
        """Apply mixed boundary conditions."""
        if boundary.marker() != pg.core.MARKER_BOUND_MIXED:
            return 0

        sourcePos = pg.center(userData['sourcePos'])

        k = userData['k']
        r1 = boundary.center() - sourcePos

        # Mirror on surface at depth=0
        r2 = boundary.center() - pg.RVector3(1.0, -1.0, 1.0) * sourcePos
        r1A = r1.abs()
        r2A = r2.abs()

        rho = 1.

        if self.resistivity is not None:
            rho = self.resistivity[boundary.leftCell().id()]

        n = boundary.norm()

        if r1A > 1e-12 and r2A > 1e-12:
            # see mod-dc-2d example for robin like BC and the negative sign
            if (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k)) > 1e-12:

                return k / rho * (r1.dot(n) / r1A * pg.math.besselK1(r1A * k) +
                                  r2.dot(n) / r2A * pg.math.besselK1(r2A * k))\
                    / (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k))
            else:
                return 0.
        else:
            return 0.

    def pointSource(self, cell, f, userData):
        r"""
        Define function for the current source term.

        :math:`\delta(x-pos), \int f(x) \delta(x-pos)=f(pos)=N(pos)`
            Right hand side entries will be shape functions(pos)
        """
        i = userData['i']
        sourcePos = userData['sourcePos'][i]

        if cell.shape().isInside(sourcePos):
            f.setVal(cell.N(cell.shape().rst(sourcePos)), cell.ids())

    def createRHS(self, mesh, elecs):
        """Create right-hand-side vector."""
        rhs = np.zeros((len(elecs), mesh.nodeCount()))
        for i, e in enumerate(elecs):
            c = mesh.findCell(e)
            rhs[i][c.ids()] = c.N(c.shape().rst(e))
        return rhs


if __name__ == "__main__":
    pass
