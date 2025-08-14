# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

Basic inversion frameworks that usually needs a forward operator to run.
"""
from math import sqrt
import numpy as np
import pygimli as pg
from pygimli.solver.leastsquares import lsqr as lssolver
from pygimli.core.trans import str2Trans
from pygimli.utils import prettyFloat as pf
from pygimli.core.matrix import asDense
from .linesearch import lineSearch


class InversionBase(object):
    """ Inversion base class for all inversions.

    New inversion frame not using RInversion anymore.
    """
    def __init__(self, fop=None, **kwargs):
        """ Initialize inversion base class.
        """
        self._debug = kwargs.pop('debug', False)
        self._verbose = kwargs.pop('verbose', False)
        self._preStep = None
        self._postStep = None
        self._fop = fop
        self._model = None
        self._response = None
        self.lam = 20      # lambda regularization
        self.lambdaFactor = 1.0
        self.minLambda = 0
        self.iter = 0
        self.debug = False
        self.robustData = False
        self.blockyModel = False
        self.chi2History = []
        self.modelHistory = []
        self.deltaPhiPercent = 1
        self.minDPhi = 1
        self.stopAtChi1 = True
        self.localRegularization = False
        self.dataTrans = kwargs.pop('dataTrans', pg.trans.TransLin())
        self.cWeight = 1
        self.axs = None  # for showProgress only
        self.LSiter = 100
        self.maxIter = kwargs.pop('maxIter', 20)
        self.G = None
        self._jacobianOutdated = False
        self.lineSearchMethod = None  # auto inter-quad
        # self.minTau/maxTau


    @property
    def fop(self):
        """ Forward operator.
        """
        return self._fop


    @fop.setter
    def fop(self, f):
        """ Set forward operator.
        """
        self.setForwardOperator(f)


    def setForwardOperator(self, fop):
        """ Set forward operator.
        """
        self._fop = fop
        # we need to initialize the regionmanager by calling it once
        self._fop.regionManager()


    @property
    def verbose(self):
        """ Verbosity level.
        """
        return self._verbose


    @verbose.setter
    def verbose(self, v):
        """ Set verbosity level for both forward operator and inversion.
        """
        self._verbose = v
        self.fop.setVerbose(self._verbose)


    @property
    def dataTrans(self):
        """ Data transformation.
        """
        return self._dataTrans


    @dataTrans.setter
    def dataTrans(self, dt):
        """ Set data transformation.
        """
        if isinstance(dt, str):
            dt = str2Trans(dt)

        self._dataTrans = dt


    @property
    def modelTrans(self):
        """ Model transformation.
        """
        return self.fop.modelTrans


    @modelTrans.setter
    def modelTrans(self, mt):
        """ Set model transformation.
        """
        if isinstance(mt, str):
            mt = str2Trans(mt)

        self.fop.modelTrans = mt  # self._modelTrans # ????


    def convertStartModel(self, model):
        """ Convert scalar or array into startmodel vector.

        Use valid range or self.fop.parameterCount, if possible.

        Arguments
        ---------
        model: float|int|array|None
            starting model value or array
        """
        if model is None:
            return None
        elif isinstance(model, float) or isinstance(model, int):
            pg.debug("Homogeneous starting model set to:", float(model))
            return np.full(self.fop.parameterCount, float(model))
        elif hasattr(model, '__iter__'):
            if len(model) == self.fop.parameterCount:
                pg.debug("Starting model set from given array.", model)
                return model
            else:
                pg.error("Starting model size invalid {0} != {1}.".
                         format(len(model), self.fop.parameterCount))
        return None


    @property
    def model(self):
        """ The last active, i.e., current model.
        """
        return self._model


    @model.setter
    def model(self, m):
        """ Set model.

        Arguments
        ---------
        m: [float,]
            Model used as current model.
        """
        if self._model is not None:
            assert len(self._model) == len(m), "Model size mismatch!"

        self._model = m
        if np.any(m - self._model):
            self._jacobianOutdated = True
            self._response = None  # not known


    @property  # not sure if we need it
    def response(self):
        """ Return last forward response.
        """
        if self._response is None:
            self._response = self.fop.response(self.model)

        return self._response


    @response.setter
    def response(self, v):
        """ Set response vector from outside (e.g. postprocessing).
        """
        # does this even make sense? visualization is Manager's job
        assert len(self.dataVals) == len(v), "Response size not matching."
        self._response = v


    @property
    def dataVals(self):
        """ Data vector (deprecated).
        """
        # why deprecated? whats the alternative?
        if self._dataVals is None:
            pg.critical("No data. Inversion framework needs data values to run")

        return self._dataVals


    @dataVals.setter
    def dataVals(self, d):
        """ Set mandatory data values.

        Values == 0.0 will be set to tolerance
        """
        # also deprecated? whats the alternative?
        if d is None:
            pg.critical("Data values can't be set to None")
        self._dataVals = d


    @property
    def errorVals(self):
        """ Errors vector (deprecated).
        """
        # why deprecated? whats the alternative?
        return self._errorVals


    @errorVals.setter
    def errorVals(self, d):
        """ Set mandatory error values.

        Values == 0.0 will be set to tolerance.
        """
        # also deprecated? whats the alternative?
        self._errorVals = d

        if self._errorVals is None:
            pg._y(d)
            pg.critical("Inversion framework needs error values to run")

        if min(abs(self._errorVals)) < 1e-12:
            print(self._errorVals)
            pg.warn(
                "Found zero error values. Setting them to fallback value of 1")
            pg.core.fixZero(self._errorVals, 1)


    def echoStatus(self):
        """ Echo inversion status (model, response, rms, chi^2, phi).
        """
        # implementation?


    def setPostStep(self, p:callable):
        """ Set a function to be called after each iteration.

        Arguments
        ---------
        p: callable
            Function to be called after each iteration.
            called with `p(IterationNumber:int, self:InversionBase)`

        """
        self._postStep = p


    def setPreStep(self, p:callable):
        """ Set a function to be called before each iteration.

        Arguments
        ---------
        p: callable
            Function to be called before each iteration.
            called with `p(IterationNumber:int, self:InversionBase)`
        """
        self._preStep = p



    def setData(self, data):
        """ Set data.
        """
        # QUESTION_ISNEEDED
        if isinstance(data, pg.DataContainer):
            raise Exception("should not be here .. its Managers job")
            self.fop.setData(data)
        else:
            self.dataVals = data


    def chi2(self, response=None):
        """ Chi-squared misfit (mean of squared error-weighted misfit).
        """
        return self.phiData(response) / len(self.dataVals)


    def phiData(self, response=None):
        """ Data objective function (sum of squared error-weighted misfit).
        """
        if response is None:
            response = self.response

        dT = self.dataTrans
        dData = (dT.trans(self.dataVals) - dT.trans(response)) / \
            dT.error(self.dataVals, self.errorVals)

        return pg.math.dot(dData, dData)


    def roughness(self, model=None, weighted=True):
        """ Return (weighted) roughness vector.
        """
        if model is None:
            model or self.model

        modelVector = self.modelTrans(model)
        pureRoughness = self.fop.constraints().mult(modelVector)
        if weighted:
            return self.cWeight * pureRoughness
        else:
            return pureRoughness


    def phiModel(self, model=None):
        """ Model objective function (norm of regularization term).
        """
        if model is None:
            model = self.model

        rough = self.roughness(model, weighted=True)
        return pg.math.dot(rough, rough)


    def phi(self, model=None, response=None):
        """ Total objective function (phiD + lambda * phiM).
        """
        if response is None:
            response = self.response

        phiD = self.phiData(response)
        if self.localRegularization:
            return phiD
        else:
            return phiD + self.phiModel(model) * self.lam


    def relrms(self):
        """ Relative root-mean-square misfit of the last run.
        """
        return pg.math.rrms(self.data, self.response)


    def absrms(self):
        """ Absolute root-mean-square misfit of the last run.
        """
        return pg.math.rms(self.data, self.response)

    def setRegularization(self, *args, **kwargs):
        """ Set regularization properties for the inverse problem.

        This can be for specific regions (args) or all regions (no args).

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
        if len(args) == 0:
            args = ('*',)

        if "operator" in kwargs:
            self.fop.setCustomConstraints(kwargs.pop("operator"))
        if "C" in kwargs:
            self.fop.setCustomConstraints(kwargs.pop("C"))

        if len(kwargs) > 0:
            self.fop.setRegionProperties(*args, **kwargs)


    def setInterRegionConstraint(self, region1, region2, strength):
        """ Set constraints between neighboring regions.

        Parameters
        ----------
        region1, region2 : int|'*'
            Region IDs
        strength : float
            weighting factor for roughness across regions
        """
        self.fop.regionManager().setInterRegionConstraint(
            region1, region2, strength)


    def setInterfaceConstraint(self, marker, strength):
        """ Set regularization strength on specific interface.

        Parameters
        ----------
        marker : int
            Boundary marker of the interface
        strength : float
            weighting factor for roughness across boundary
        """
        self.fop.regionManager().setInterfaceConstraint(
            marker, strength)


    def setConstraintWeights(self, cWeight):
        """ Set weighting factors for the individual rows of the C matrix.
        """
        self.cWeight = cWeight


    def reset(self):
        """ Reset function currently called at beginning of every inversion.
        """
        # FW: Note that this is called at the beginning of run. I therefore
        # removed the startingModel here to allow explicitly set starting
        # models by the user.
        # if self._keepStartModel is False:
        #     self._startModel = None

        self._model = None
        self._dataVals = None
        self._errorVals = None


    def oneStep(self):
        """ Carry out one iteration step (e.g. good for coupling etc.).
        """
        dModel = self.modelUpdate()
        if self.verbose is True:
            print("dM: ", dModel)
        tau, responseLS = lineSearch(self, dModel)

        if self.verbose is True:
            print("tau: ", tau)

        pg.debug(f"tau={tau}")
        if tau >= 0.95:  # practically 1
            tau = 1

        self.model = self.modelTrans.update(self.model, dModel*tau)
        if tau == 1.0:
            self.response = responseLS
        else:  # compute new response
            self.response = self.fop.response(self.model)


    def start(self, dataVals, model=None, errorVals=None, response=None, **kwargs):
        """Initialize inversion run by setting data, model and response.

        Parameters
        ----------
        dataVals : iterable
            Data values 
        model : iterable
            Model values, if not given, the fops start model will be used.
        errorVals : iterable
            Relative error values. dv / v
            Can be omitted if absoluteError and/or relativeError kwargs given   
        response : iterable
            Response vector, if not given, it is computed from the model.

        Keyword Arguments
        -----------------
        absoluteError : float | iterable
            absolute error in units of dataVals
        relativeError : float | iterable
            relative error related to dataVals
        """
        self.dataVals = dataVals

        if errorVals is None:  # use absoluteError and/or relativeError instead
            absErr = kwargs.pop("absoluteError", 0)
            relErr = kwargs.pop("relativeError", 0)
            if np.any(np.isclose(absErr + relErr, 0, atol=0)):
                raise Exception("Zero error occurred, check abs/relErr")

            errorVals = pg.abs(absErr / np.asarray(dataVals)) + relErr

        if isinstance(errorVals, (float, int)):
            errorVals = np.ones_like(dataVals) * errorVals

        # Triggers update of fop properties, any property to be set before.
        self.errorVals = errorVals

        if self.fop is None:
            raise Exception("Need valid forward operator for inversion run.")

        if model is None:
            model = self.fop.startModel()
        
        self.model = model
        if response is None:
            self.response = self.fop.response(model)

    def run(self, dataVals, errorVals=None, **kwargs):
        """ Run inversion.

        The inversion will always start from the starting model taken from
        the forward operator.
        If you want to run the inversion from a specified prior model,
        e.g., from a other run, set this model as starting model to the FOP
        (fop.setStartModel).

        Parameters
        ----------
        dataVals : iterable
            Data values
        errorVals : iterable
            Relative error values. dv / v
            Can be omitted if absoluteError and/or relativeError kwargs given

        Keyword Arguments
        -----------------
        absoluteError : float | iterable
            absolute error in units of dataVals
        relativeError : float | iterable
            relative error related to dataVals
        maxIter : int
            Overwrite class settings for maximal iterations number.
        dPhi : float [1]
            Overwrite class settings for delta data phi aborting criteria.
            Default is 1%
        cType: int [1]
            Temporary global constraint type for all regions.
        startModel: array
            Temporary starting model for the current inversion run.
        lam: float
            Temporary regularization parameter lambda.
        lambdaFactor : float [1]
            Factor to change lam with every iteration
        robustData : bool
            Robust (L1 norm mimicking) data reweighting
        blockyModel : bool
            Robust (L1 norm mimicking) model roughness reweighting
        isReference : bool [False]
            Starting model is also a reference to constrain against
        showProgress : bool
            Show progress in form of updating models
        verbose : bool
            Verbose output on the console
        debug : bool
            Even more verbose console and file output
        """
        self.reset()
        startModel = self.convertStartModel(kwargs.pop('startModel', None))
        self.start(dataVals, errorVals=errorVals, model=startModel, **kwargs)
        self.fop.setVerbose(False)  # gets rid of CHOLMOD messages
        maxIter = kwargs.pop('maxIter', self.maxIter)
        minDPhi = kwargs.pop('dPhi', self.minDPhi)
        showProgress = kwargs.pop('showProgress', False)
        if 'blockyModel' in kwargs:
            self.blockyModel = kwargs['blockyModel']

        self.verbose = kwargs.pop('verbose', self.verbose)
        self.debug = kwargs.pop('debug', self.debug)
        self.robustData = kwargs.pop('robustData', False)
        self.stopAtChi1 = kwargs.pop("stopAtChi1", True)

        self.lam = kwargs.pop('lam', self.lam)
        self.lambdaFactor = kwargs.pop('lambdaFactor', 1.0)

        # catch a few regularization options that don't go into run
        for opt in ["cType", "limits", "correlationLengths", "C"]:
            if opt in kwargs:
                di = {opt: kwargs.pop(opt)}
                pg.verbose("Set regularization", di)
                self.setRegularization(**di)

        # temporary set max iter to one for the initial run call
        maxIterTmp = self.maxIter
        self.maxIter = 1

        if self.verbose:
            pg.info('Starting inversion.')
            print("fop:", self.fop)
            print("Data transformation:", self.dataTrans)
            print("Model transformation:", self.modelTrans)

            print("min/max (data): {0}/{1}".format(pf(min(self._dataVals)),
                                                   pf(max(self._dataVals))))
            print("min/max (error): {0}%/{1}%".format(
                pf(100 * min(self._errorVals)),
                pf(100 * max(self._errorVals))))
            print("min/max (start model): {0}/{1}".format(
                pf(min(startModel)), pf(max(startModel))))

        # To ensure reproduceability of the run() call, inv.start() will
        # reset self.inv.model() to fop.startModel().
        self.fop.setStartModel(startModel)
        if kwargs.pop("isReference", False):
            self.referenceModel = startModel
            pg.info("Setting starting model as reference!")

        self.model = pg.Vector(startModel)
        if self.verbose:
            print("-" * 80)
        if self._preStep and callable(self._preStep):
            self._preStep(0, self)

        self.maxIter = 0
        # self.start() # what's done here?
        self.maxIter = maxIterTmp
        if self.verbose:
            print("inv.iter 0 ... chi² = {:7.2f}".format(self.chi2()))
            # print("inv.iter 0 ... chi² = {0}".format(round(self.chi2(), 2)))

        if self._postStep and callable(self._postStep):
            self._postStep(0, self)

        if showProgress:
            self.showProgress(showProgress)

        # self.fop.checkConstraints()
        lastPhi = self.phi()
        self.chi2History = [self.chi2()]
        self.modelHistory = [startModel]

        for i in range(1, maxIter+1):
            if self._preStep and callable(self._preStep):
                self._preStep(i, self)

            if self.verbose:
                print("-" * 80)
                print("inv.iter", i, "... ", end='')

            self.oneStep()
            if np.isnan(self.model).any():
                pg.info(self.model)
                pg.critical('invalid model')

            chi2 = self.chi2()

            self.chi2History.append(chi2)
            self.modelHistory.append(self.model)

            if showProgress:
                self.showProgress(showProgress)

            if self._postStep and callable(self._postStep):
                self._postStep(i, self)

            if self.robustData:
                pass
                # self.inv.robustWeighting()

            if self.blockyModel:
                pass
                # self.inv.constrainBlocky()

            phi = self.phi()
            dPhi = phi / lastPhi

            if self.verbose:
                print("chi² = {:7.2f} (dPhi = {:.2f}%) lam: {:.1f}".format(
                    chi2, (1 - dPhi) * 100, self.lam))

            if chi2 <= 1 and self.stopAtChi1:
                print("\n")
                if self.verbose:
                    pg.boxprint(
                        "Abort criterion reached: chi² <= 1 (%.2f)" % chi2)
                break

            # if dPhi < -minDPhi:
            if (dPhi > (1.0 - minDPhi / 100.0)) and i > 2:  # should be minIter
                if self.verbose:
                    pg.boxprint(
                        "Abort criterion reached: dPhi = {0} (< {1}%)".format(
                            round((1 - dPhi) * 100, 2), minDPhi))
                break

            lastPhi = phi
            self.lam = max(self.lam*self.lambdaFactor, self.minLambda)

        # will never work as expected until we unpack kwargs .. any idea for
        # better strategy?
        # if len(kwargs.keys()) > 0:
        #     print("Warning! unused keyword arguments", kwargs)

        return self.model


    def jacobianMatrix(self, error_weighted=False, numpy_matrix=False):
        """ Jacobian matrix of the inverse (data/model-transformed) problem.

        Whereas the forward operator holds the Jacobian matrix of the forward,
        i.e. the intrinsic (untransformed) problem, this function returns the
        Jacobian of the (transformed) inverse problem, i.e. taking model and
        data transformations into account by using their inner derivatives.

        Parameters
        ----------
        self : pg.Inversion
            inversion instance with model, response and fop.jacobian
        error_weighted : bool
            add error weighting according to data transform
        numpy_matrix : bool
            return numpy matrix instead of MultLeftRightMatrix
        """
        tData = self.dataTrans.deriv(self.response)
        tModel = 1 / self.modelTrans.deriv(self.model)
        if error_weighted:
            tData /= self.dataTrans.error(self.response, self.errorVals)
        if numpy_matrix:
            return np.reshape(tData, [-1, 1]) * \
                pg.utils.gmat2numpy(self.fop.jacobian()) * \
                np.reshape(tModel, [1, -1])
        else:
            return pg.matrix.MultLeftRightMatrix(self.fop.jacobian(),
                                                tData, tModel)


    def residual(self):
        """ Residual vector (data-reponse)/error using data transform.
        """
        return (self.dataTrans.fwd(self.dataVals) -
                self.dataTrans.fwd(self.response)) / \
            self.dataTrans.error(self.response, self.errorVals)

    def dataGradientFormal(self):  # formal but restricted to existent J
        """ Data gradient from Jacobian and residual, i.e. J^T * dData.
        """
        return -self.jacobianMatrix(error_weighted=True).transMult(
            self.residual())


    def dataGradient(self, error_weighted=True):  # also works for fop.STy
        """ Data gradient from jacobian and residual, i.e. J^T * dData.
        """
        tData = self.dataTrans.deriv(self.response)
        if error_weighted:
            tData /= self.dataTrans.error(self.response, self.errorVals)

        return self.fop.STy(-self.residual()*tData) / \
            self.modelTrans.deriv(self.model)


    def modelGradient(self):
        """ Model gradient, i.e. C^T * C * (m - m0).
        """
        # self.inv.checkConstraints() # not necessary?
        if isinstance(self.cWeight, (float, int)):
            C = pg.matrix.ScaledMatrix(self.fop.constraints(), self.cWeight)
        else:
            C = pg.matrix.MultLeftMatrix(self.fop.constraints(), self.cWeight)

        return C.transMult(C.mult(self.modelTrans(self.model)))


    def gradient(self):
        """ Gradient of the objective function.
        """
        return self.dataGradient() + self.modelGradient() * self.lam


class GaussNewtonInversion(InversionBase):
    """ Gauss-Newton based inversion.
    """
    def __init__(self, fop=None, **kwargs):
        """ Initialize Gauss-Newton inversion.

        Arguments
        ---------
        fop : ForwardOperator
            Forward operator to be used for the inversion.
        """
        super().__init__(fop=fop, **kwargs)


    def modelUpdate(self):
        """ Compute (full) model update from inverse.
        """
        pg.verbose("Running LSQR inversion step!")
        model = self.model
        if len(self.response) != len(self.dataVals):
            self.setResponse(self.fop.response(model))

        self.fop.createJacobian(model)
        # self.checkTransFunctions()
        tD = self.dataTrans
        tM = self.modelTrans
        nData = len(self.dataVals)
        self.A = pg.BlockMatrix()  # to be filled with scaled J and C matrices

        ### part 1: data part
        self.JJ = self.jacobianMatrix(error_weighted=True)
        self.mat1 = self.A.addMatrix(self.JJ)
        self.A.addMatrixEntry(self.mat1, 0, 0)

        ### part 2: normal constraints
        # self.checkConstraints()
        self.C = self.fop.constraints()
        self.leftC = pg.Vector(self.C.rows(), 1.0)
        self.rightC = pg.Vector(self.C.cols(), 1.0)
        self.CC = pg.matrix.MultLeftRightMatrix(self.C,
                                                self.leftC, self.rightC)
        self.mat2 = self.A.addMatrix(self.CC)
        self.A.addMatrixEntry(self.mat2, nData, 0, sqrt(self.lam))

        ### part 3: parameter constraints
        if self.G is not None:
            self.rightG = 1.0 / tM.deriv(model)
            self.GG = pg.matrix.MultRightMatrix(self.G, self.rightG)
            self.mat3 = self.A.addMatrix(self.GG)
            nConst = self.C.rows()
            self.A.addMatrixEntry(self.mat3, nData+nConst, 0, sqrt(self.my))

        self.A.recalcMatrixSize()
        # right-hand side vector
        # deltaD = (tD.fwd(self.dataVals)-tD.fwd(self.response)) * self.dScale
        deltaD = self.residual()
        deltaC = -(self.CC * tM.fwd(model) * sqrt(self.lam))
        deltaC *= 1.0 - self.localRegularization  # oper. on DeltaM only
        rhs = pg.cat(deltaD, deltaC)
        if self.G is not None:
            deltaG = (self.c - self.G * model) * sqrt(self.my)
            rhs = pg.cat(rhs, deltaG)

        dM = lssolver(self.A, rhs, maxiter=self.LSiter, verbose=self.verbose)
        return dM


class DescentInversion(InversionBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def modelUpdate(self):
        """The negative gradient of  objective function as search direction."""
        if self.fop.STy.__doc__ is pg.Modelling.STy.__doc__:  # original
            if len(self.model) != self.fop.jacobian().cols():
                self._jacobianOutdated = True

            if self._jacobianOutdated:
                self.fop.createJacobian(self.model)

        return -self.gradient()


class NLCGInversion(InversionBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def modelUpdate(self):
        return None


class LBFGSInversion(InversionBase):
    """Limited-memory BFGS minimization."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def modelUpdate(self):
        return None


# Note that there is a lot of redundancy but this class is to be removed upon
# thorough testing of the purely Python Gauss-NewtonInversion class.
class ClassicInversion(object):
    """Basic Gauss-Newton based inversion framework (to be removed).

    Changes to prior Versions (remove me)

        * holds the starting model itself, forward operator only provides a
        method to create it: fop.createStartModel(dataValues)

    Attributes
    ----------
    verbose : bool
        Give verbose output
    debug : bool
        Give debug output
    startModel : float|array|None
        Current starting model that can be set in the init or as property.
        If not set explicitly, it will be estimated from the forward operator
        methods. This property will be recalulated for every run call if not
        set explicitly with self.startModel = float|array, or None to reforce
        autogeneration. Note that the run call accepts a temporary startModel
        (for the current calculation only).
    model : array
        Holds the last active model
    maxIter : int [20]
        Maximal interation number.
    stopAtChi1 : bool [True]
        Stop iteration when chi² is one. If set to False the iteration stops
        after maxIter or convergence reached (self.inv.deltaPhiAbortPercent())
    """
    def __init__(self, fop=None, inv=None, **kwargs):
        self._debug = kwargs.pop('debug', False)
        self._verbose = kwargs.pop('verbose', False)

        # If this class or its derived is a Framework the _inv holds another
        # Inversion which allows us (remove me)........
        # this will be probably removed in the future
        self.isFrameWork = False  # check if needed
        self._stopAtChi1 = True

        self._preStep = None
        self._postStep = None

        self._inv = None
        self._fop = None
        self._lam = 20      # lambda regularization
        self.chi2History = []

        # cache: keep startmodel if set explicitly or calculated from FOP, will
        # be recalulated for every run if not set explicitly
        self._startModel = None
        # flag to keep startModel if set manually by init or self.startModel
        # unless self.startModel = None
        self._keepStartModel = False

        self.reset()

        if inv is not None:
            self._inv = inv
            self.isFrameWork = True
        else:
            self._inv = pg.core.RInversion(self._verbose, self._debug)

        self.dataTrans = kwargs.pop('dataTrans', pg.trans.TransLin())
        self.axs = None  # for showProgress only
        self.maxIter = kwargs.pop('maxIter', 20)

        if fop is not None:
            self.setForwardOperator(fop)

        if "startModel" in kwargs:
            self.startModel = kwargs["startModel"]

    def reset(self):
        """Reset function currently called at beginning of every inversion."""
        # FW: Note that this is called at the beginning of run. I therefore
        # removed the startingModel here to allow explicitly set starting
        # models by the user.
        if self._keepStartModel is False:
            self._startModel = None

        self._model = None
        self._dataVals = None
        self._errorVals = None

    @property
    def iter(self):
        """Number of iterations."""
        return self._inv.iter()

    @property
    def inv(self):
        """Return (core) inversion object."""
        return self._inv

    @property
    def fop(self):
        """Forward operator."""
        return self._fop

    @fop.setter
    def fop(self, f):
        """Set forward operator."""
        self.setForwardOperator(f)

    def setForwardOperator(self, fop):
        """Set forward operator."""
        self._fop = fop
        # we need to initialize the regionmanager by calling it once
        self._fop.regionManager()
        self._inv.setForwardOperator(fop)

    @property
    def verbose(self):
        """Verbosity level."""
        return self._verbose

    @verbose.setter
    def verbose(self, v):
        """Set verbosity level for both forward operator and inversion."""
        self._verbose = v
        self.inv.setVerbose(self._verbose)
        self.fop.setVerbose(self._verbose)

    @property
    def debug(self):
        """Debug level."""
        return self._debug

    @debug.setter
    def debug(self, v):
        """Set debug (output and files) level for both inversion."""
        self._debug = v
        self.inv.setDoSave(self._debug)

    @property
    def dataTrans(self):
        """Data transformation."""
        return self._dataTrans

    @dataTrans.setter
    def dataTrans(self, dt):
        """Set data transformation."""
        if isinstance(dt, str):
            dt = str2Trans(dt)

        self._dataTrans = dt
        self.inv.setTransData(self._dataTrans)

    @property
    def modelTrans(self):
        """Model transformation."""
        return self.fop.modelTrans

    @modelTrans.setter
    def modelTrans(self, mt):
        """Set model transformation."""
        if isinstance(mt, str):
            mt = str2Trans(mt)

        self.fop.modelTrans = mt  # self._modelTrans # ????

    @property
    def startModel(self):
        """Return current default starting model.

        Returns the current default starting model or
        calls `fop.createStartmodel()` if none is defined.
        """
        if self._startModel is None:
            sm = self.fop.regionManager().createStartModel()
            if len(sm) > 0 and max(abs(np.atleast_1d(sm))) > 0.0:
                self._startModel = sm
                if self.verbose:
                    pg.info("Created startmodel from region infos:", sm)
            else:
                if self.verbose:
                    pg.verbose("No region infos for startmodel")

        if self._startModel is None:
            sm = self.fop.createStartModel(self.dataVals)
            # pg.info("Created startmodel from forward operator:", sm)
            if self.verbose:
                pg.info("Created startmodel from forward operator:" +
                        "{:d}, min/max={:f}/{:f}".format(
                            len(sm), min(sm), max(sm)))
            self._startModel = sm
        return self._startModel

    @startModel.setter
    def startModel(self, model):
        """Set starting model.

        model: [float] | float
            Model used as starting model.
            Float value is used as constant model.
        """
        sm = self.convertStartModel(model)
        if sm is None:
            self._keepStartModel = False
        else:
            self._keepStartModel = True
        self._startModel = sm

    def convertStartModel(self, model):
        """Convert scalar or array into startmodel vector.

        Use valid range or self.fop.parameterCount, if possible.

        Attributes
        ----------
        model: float|int|array|None
            starting model value or array
        """
        if model is None:
            return None
        elif isinstance(model, float) or isinstance(model, int):
            pg.debug("Homogeneous starting model set to:", float(model))
            return np.full(self.fop.parameterCount, float(model))
        elif hasattr(model, '__iter__'):
            if len(model) == self.fop.parameterCount:
                pg.debug("Starting model set from given array.", model)
                return model
            else:
                pg.error("Starting model size invalid {0} != {1}.".
                         format(len(model), self.fop.parameterCount))
        return None

    @property
    def model(self):
        """The last active (i.e. current) model."""
        if self._model is None:
            if hasattr(self.inv, 'model'):
                # inv is RInversion()
                if len(self.inv.model()) > 0:
                    return self.inv.model()
                else:
                    raise Exception(pg.critical(
                        "There was no inversion run so there is no last model"))
            else:
                return self.inv.model
        return self._model

    @model.setter
    def model(self, m):
        self._model = m

    @property
    def response(self):
        """Return last forward response."""
        if len(self.inv.response()) > 0:
            return self.inv.response()
        else:
            raise Exception(
                "There was no inversion run so there is no response yet")

    @response.setter
    def response(self, v):
        """Set response vector in the inside RInversion."""
        self.inv.setResponse(v)

    # backward compatibility
    @property
    def dataErrs(self):
        """Data errors (deprecated)."""
        pg.warn('do not use')
        return self._errorVals

    @dataErrs.setter
    def dataErrs(self, v):
        """Set data errors (deprecated)."""
        pg.warn('do not use')
        self._errorVals = v

    @property
    def dataVals(self):
        """Data vector (deprecated)."""
        return self._dataVals

    @dataVals.setter
    def dataVals(self, d):
        """Set mandatory data values.

        Values == 0.0 will be set to tolerance
        """
        self._dataVals = d

        if self._dataVals is None:
            pg._y(d)
            pg.critical("Inversion framework needs data values to run")

        # zero can be a valid data value
        #
        # if min(abs(self._dataVals)) < 1e-12:
        #     print(self._dataVals)
        #     pg.warn("Found zero data values. \
        #             Setting them to a TOLERANCE value of 1e-12")
        #     pg.core.fixZero(self._dataVals, 1e-12)

    @property
    def errorVals(self):
        """Errors vector (deprecated)."""
        return self._errorVals

    @errorVals.setter
    def errorVals(self, d):
        """Set mandatory error values.

        Values == 0.0. Will be set to Tolerance
        """
        self._errorVals = d

        if self._errorVals is None:
            pg._y(d)
            pg.critical("Inversion framework needs error values to run")

        if min(abs(self._errorVals)) < 1e-12:
            print(self._errorVals)
            pg.warn(
                "Found zero error values. Setting them to fallback value of 1")
            pg.core.fixZero(self._errorVals, 1)

    @property
    def robustData(self):
        """Return robust data reweighting (IRLS L1 scheme) bool."""
        return self.inv.robustData()

    @robustData.setter
    def robustData(self, v):
        """Set robust data reweighting (IRLS L1 scheme) True or False."""
        if self.inv is not None:
            self.inv.setRobustData(v)

    @property
    def blockyModel(self):
        """Return blocky model roughness reweighting (IRLS L1 scheme) bool."""
        return self.inv.blockyModel()

    @blockyModel.setter
    def blockyModel(self, v):
        """Set blocky model roughness reweighting (IRLS L1 scheme) bool."""
        if self.inv is not None:
            self.inv.setBlockyModel(v)

    @property
    def maxIter(self):
        """Maximum iterations."""
        return self.inv.maxIter()

    @maxIter.setter
    def maxIter(self, v):
        """Set maximum iterations."""
        if self.inv is not None:
            self.inv.setMaxIter(v)

    @property
    def stopAtChi1(self):
        """Stop at chi^2=1 behaviour (bool)."""
        return self._stopAtChi1

    @stopAtChi1.setter
    def stopAtChi1(self, b):
        """Define whether to stop at chi^2=1."""
        self._stopAtChi1 = b

    @property
    def minDPhi(self):
        """Minimum data objective function decrease."""
        return self.inv.deltaPhiAbortPercent()

    @minDPhi.setter
    def minDPhi(self, dPhi):
        """Set minimum data objective function decrease."""
        return self.setDeltaChiStop(dPhi)

    @property
    def lam(self):
        """ Return regularization strength.
        """
        return self._lam

    @lam.setter
    def lam(self, lam):
        """ Set regularization strength.
        """
        self._lam = lam
        self.inv.setLambda(lam)

    def setDeltaPhiStop(self, it):
        """ Define minimum relative decrease in objective function to stop.
        """
        self.inv.setDeltaPhiAbortPercent(it)

    @pg.renamed(setDeltaPhiStop)
    def setDeltaChiStop(self, it):
        """ Set data fit change level (deprecated).
        """
        self.setDeltaPhiStop(it)

    def echoStatus(self):
        """ Echo inversion status (model, response, rms, chi^2, phi).
        """
        self.inv.echoStatus()

    def setPostStep(self, p):
        """ Set a function to be called after each iteration.

        The function is called with p(IterationNumber, self).
        """
        self._postStep = p

    def setPreStep(self, p):
        """ Set a function to be called before each iteration.
        """
        self._preStep = p

    def setData(self, data):
        """ Set data.
        """
        # QUESTION_ISNEEDED
        if isinstance(data, pg.DataContainer):
            raise Exception("should not be here .. its Managers job")
            self.fop.setData(data)
        else:
            self.dataVals = data

    def chi2(self, response=None):
        """ Chi-squared misfit (mean of squared error-weighted misfit).
        """
        return self.phiData(response) / len(self.dataVals)

    def phiData(self, response=None):
        """ Data objective function (sum of suqred error-weighted misfit).
        """
        if response is None:
            response = self.response

        dT = self.dataTrans
        dData = (dT.trans(self.dataVals) - dT.trans(response)) / \
            dT.error(self.dataVals, self.errorVals)

        return pg.math.dot(dData, dData)

    def phiModel(self, model=None):
        """ Model objective function (norm of regularization term).
        """
        if model is None:
            model = self.model

        rough = self.inv.roughness(model)
        return pg.math.dot(rough, rough)

    def phi(self, model=None, response=None):
        """ Total objective function (phiD + lambda * phiM).
        """
        phiD = self.phiData(response)
        if self.inv.localRegularization():
            return phiD
        else:
            return phiD + self.phiModel(model) * self.inv.getLambda()

    def relrms(self):
        """ Relative root-mean-square misfit of the last run.
        """
        return self.inv.relrms()

    def absrms(self):
        """ Absolute root-mean-square misfit of the last run.
        """
        return self.inv.absrms()

    def setRegularization(self, *args, **kwargs):
        """ Set regularization properties for the inverse problem.

        This can be for specific regions (args) or all regions (no args).

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
        if len(args) == 0:
            args = ('*',)

        if "operator" in kwargs:
            self.fop.setCustomConstraints(kwargs.pop("operator"))
        if "C" in kwargs:
            self.fop.setCustomConstraints(kwargs.pop("C"))

        if len(kwargs) > 0:
            self.fop.setRegionProperties(*args, **kwargs)

    def setInterRegionConstraint(self, region1, region2, strength):
        """ Set constraints between neighboring regions.

        Parameters
        ----------
        region1, region2 : int|'*'
            Region IDs
        strength : float
            weighting factor for roughness across regions
        """
        self.fop.regionManager().setInterRegionConstraint(
            region1, region2, strength)

    def setInterfaceConstraint(self, marker, strength):
        """ Set regularization strength on specific interface.

        Parameters
        ----------
        marker : int
            Boundary marker of the interface
        strength : float
            weighting factor for roughness across boundary
        """
        self.fop.regionManager().setInterfaceConstraint(
            marker, strength)

    def setConstraintWeights(self, cWeight):
        """ Set weighting factors for the invidual rows of the C matrix.
        """
        self.inv.setCWeight(cWeight)

    def run(self, dataVals, errorVals=None, **kwargs):
        """ Run inversion.

        The inversion will always start from the starting model taken from
        the forward operator.
        If you want to run the inversion from a specified prior model,
        e.g., from a other run, set this model as starting model to the FOP
        (fop.setStartModel).
        Any self.inv.setModel() settings will be overwritten.

        Parameters
        ----------
        dataVals : iterable
            Data values
        errorVals : iterable
            Relative error values. dv / v
            Can be omitted if absoluteError and/or relativeError kwargs given

        Keyword Arguments
        -----------------
        absoluteError : float | iterable
            absolute error in units of dataVals
        relativeError : float | iterable
            relative error related to dataVals
        maxIter : int
            Overwrite class settings for maximal iterations number.
        dPhi : float [1]
            Overwrite class settings for delta data phi aborting criteria.
            Default is 1%
        cType: int [1]
            Temporary global constraint type for all regions.
        startModel: array
            Temporary starting model for the current inversion run.
        lam: float
            Temporary regularization parameter lambda.
        lambdaFactor : float [1]
            Factor to change lam with every iteration
        robustData : bool
            Robust (L1 norm mimicking) data reweighting
        blockyModel : bool
            Robust (L1 norm mimicking) model roughness reweighting
        isReference : bool [False]
            Starting model is also a reference to constrain against
        showProgress : bool
            Show progress in form of updating models
        verbose : bool
            Verbose output on the console
        debug : bool
            Even more verbose console and file output
        """
        with pg.tictoc('run'):
            with pg.tictoc('prep'):
                with pg.tictoc('reset'):
                    self.reset()

                if errorVals is None:  # use absoluteError and/or relativeError instead
                    absErr = kwargs.pop("absoluteError", 0)
                    relErr = kwargs.pop("relativeError",
                                        0.01 if np.allclose(absErr, 0) else 0)
                    errorVals = pg.abs(absErr / np.asarray(dataVals)) + relErr

                if isinstance(errorVals, (float, int)):
                    errorVals = np.ones_like(dataVals) * errorVals

                # if self.isFrameWork: # removed 240317
                #     pg.critical('in use?')
                #     return self._inv.run(dataVals, errorVals, **kwargs)

                if self.fop is None:
                    pg.critical("Need a valid forward operator for the inversion run.")

                self.fop.setVerbose(False)  # gets rid of CHOLMOD messages
                maxIter = kwargs.pop('maxIter', self.maxIter)
                minDPhi = kwargs.pop('dPhi', self.minDPhi)
                showProgress = kwargs.pop('showProgress', False)

                if 'blockyModel' in kwargs:
                    self.blockyModel = kwargs['blockyModel']

                self.verbose = kwargs.pop('verbose', self.verbose)
                self.debug = kwargs.pop('debug', self.debug)
                self.robustData = kwargs.pop('robustData', False)

                if "stopAtChi1" in kwargs:
                    self._stopAtChi1 = kwargs["stopAtChi1"]

                lam = kwargs.pop('lam', self.lam)
                self.inv.setLambda(lam)
                self.inv.setLambdaFactor(kwargs.pop('lambdaFactor', 1.0))

                # catch a few regularization options that don't go into run
                for opt in ["cType", "limits", "correlationLengths", "C"]:
                    if opt in kwargs:
                        di = {opt: kwargs.pop(opt)}
                        pg.verbose("Set regularization", di)
                        self.setRegularization(**di)

                # Triggers update of fop properties, any property to be set before.
                self.inv.setTransModel(self.fop.modelTrans)  # why from fop??
                self.dataVals = dataVals
                self.errorVals = errorVals

                self.inv.setData(self._dataVals)
                self.inv.setRelativeError(self._errorVals)

                # temporary set max iter to one for the initial run call
                maxIterTmp = self.maxIter
                self.maxIter = 1

                startModel = self.convertStartModel(kwargs.pop('startModel', None))
                # we cannot add the following into kwargs.pop, since someone may call
                # with explicit startModel=None
                if startModel is None:
                    startModel = self.startModel

                if self.verbose:
                    pg.info('Starting inversion.')
                    print("fop:", self.inv.fop())
                    if isinstance(self.inv.transData(), pg.trans.TransCumulative):
                        print("Data transformation (cumulative):")
                        for i in range(self.inv.transData().size()):
                            print("\t", i, self.inv.transData().at(i))
                    else:
                        print("Data transformation:", self.inv.transData())

                    if isinstance(self.inv.transModel(), pg.trans.TransCumulative):
                        print("Model transformation (cumulative):")
                        for i in range(self.inv.transModel().size()):
                            if i < 10:
                                print("\t", i, self.inv.transModel().at(i))
                            else:
                                print(".", end='')
                    else:
                        print("Model transformation:", self.inv.transModel())

                    print("min/max (data): {0}/{1}".format(pf(min(self._dataVals)),
                                                        pf(max(self._dataVals))))
                    print("min/max (error): {0}%/{1}%".format(
                        pf(100 * min(self._errorVals)),
                        pf(100 * max(self._errorVals))))
                    print("min/max (start model): {0}/{1}".format(
                        pf(min(startModel)), pf(max(startModel))))

                # To ensure reproducibility of the run() call, inv.start() will
                # reset self.inv.model() to fop.startModel().
                self.fop.setStartModel(startModel)

                if kwargs.pop("isReference", False):
                    self.inv.setReferenceModel(startModel)
                    pg.info("Setting starting model as reference!")

                if self.verbose:
                    print("-" * 80)
                if self._preStep and callable(self._preStep):
                    self._preStep(0, self)

        if self._postStep and callable(self._postStep):
            self._postStep(0, self)

        if showProgress:
            self.showProgress(showProgress)

        lastPhi = self.phi()
        self.chi2History = [self.chi2()]
        self.modelHistory = [startModel]

        for i in range(1, maxIter+1):
            if self._preStep and callable(self._preStep):
                self._preStep(i, self)

            self.maxIter = maxIterTmp
            if self.verbose:
                print("inv.iter 0 ... chi² = {:7.2f}".format(self.chi2()))
                # print("inv.iter 0 ... chi² = {0}".format(round(self.chi2(), 2)))

            if self._postStep and callable(self._postStep):
                self._postStep(0, self)

            if showProgress:
                self.showProgress(showProgress)

            if self._postStep and callable(self._postStep):
                self._postStep(i, self)

            for i in range(1, maxIter+1):
                with pg.tictoc('iter'):
                    if self._preStep and callable(self._preStep):
                        self._preStep(i + 1, self)

                    if self.verbose:
                        print("-" * 80)
                        print("inv.iter", i, "... ", end='')

                    self.oneStep()

                    if np.isnan(self.model).any():
                        print(self.model)
                        pg.critical('invalid model')

                    # resp = self.inv.response()  # NOT USED
                    chi2 = self.inv.chi2()

                    self.chi2History.append(chi2)
                    self.modelHistory.append(self.model)

                    if showProgress:
                        self.showProgress(showProgress)

                    if self._postStep and callable(self._postStep):
                        self._postStep(i + 1, self)

                    if self.robustData:
                        self.inv.robustWeighting()

                    if self.blockyModel:
                        self.inv.constrainBlocky()

                    phi = self.phi()
                    dPhi = phi / lastPhi

                    if self.verbose:
                        print("chi² = {:7.2f} (dPhi = {:.2f}%) lam: {:.1f}".format(
                            chi2, (1 - dPhi) * 100, lam))

                    if chi2 <= 1 and self.stopAtChi1:
                        print("\n")
                        if self.verbose:
                            pg.boxprint(
                                "Abort criterion reached: chi² <= 1 (%.2f)" % chi2)
                        break

                    # if dPhi < -minDPhi:
                    if (dPhi > (1.0 - minDPhi / 100.0)) and i > 2:  # should be minIter
                        if self.verbose:
                            pg.boxprint(
                                "Abort criterion reached: dPhi = {0} (< {1}%)".format(
                                    round((1 - dPhi) * 100, 2), minDPhi))
                        break

                    lastPhi = phi

                    lam *= self.inv.lambdaFactor()
                    self.inv.setLambda(lam)

            # will never work as expected until we unpack kwargs .. any idea for
            # better strategy?
            # if len(kwargs.keys()) > 0:
            #     print("Warning! unused keyword arguments", kwargs)

            self.model = self.inv.model()
            return self.model


    def oneStep(self):
        r"""Run one single iteration step.
        """
        with pg.tictoc('onestep'):
            try:
                self.inv.oneStep()
            except RuntimeError as e:
                print(e)
                pg.error('One step failed. '
                        'Aborting and going back to last model')


    def showProgress(self, style='all'):
        r""" Show inversion progress after every iteration.

        Can show models if `drawModel` method exists. The default fallback is
        plotting the :math:`\chi^2` fit as a function of iterations. Called if
        `showProgress=True` is set for the inversion run.
        """
        if self.fop.drawModel is None:
            style = 'convergence'

        if self.axs is None:
            axs = None
            if style == 'all' or style is True:
                fig, axs = pg.plt.subplots(1, 2)
            else:
                fig, axs = pg.plt.subplots(1, 1)
            self.axs = axs
        ax = self.axs

        if style == 'model':
            for other_ax in ax.figure.axes:
                # pg._y(type(other_ax).mro())
                if type(other_ax).mro()[0] == type(ax):
                    # only clear axes not the colorbars
                    other_ax.clear()

            self.fop.drawModel(ax, self.inv.model())
        elif style == 'all':
            # for other_ax in ax[0].figure.axes:
            #     other_ax.clear()
            for _ax in self.axs:
                _ax.clear()
                try:
                    pg.viewer.mpl.twin(_ax).clear()
                except Exception:
                    pass

            self.fop.drawModel(ax[0], self.inv.model(),
                               label='Model')
            self.fop.drawData(ax[1], self._dataVals, self._errorVals,
                              label='Data')
            self.fop.drawData(ax[1], self.inv.response(),
                              label='Response')

            ax[1].text(
                0.99, 0.005, r"Iter: {0}, rrms: {1}%, $\chi^2$: {2}".format(
                    self.inv.iter(), pf(self.inv.relrms()),
                    pf(self.inv.chi2())),
                transform=ax[1].transAxes,
                horizontalalignment='right',
                verticalalignment='bottom',
                fontsize=8)

            ax[1].figure.tight_layout()

        elif style == 'convergence':
            ax.semilogy(self.inv.iter(), self.inv.chi2(), "ro")
            if self.inv.iter() == 1:
                ax.set_xlabel("Iteration")
                ax.set_ylabel(r"$\chi^2$")
                ax.autoscale(tight=True)
                ax.axhline(y=1, ls="--")

        pg.plt.pause(0.05)

    def jacobianMatrix(self, error_weighted=False, numpy_matrix=False):
        """ Jacobian matrix of the inverse (data/model-transformed) problem.

        Whereas the forward operator holds the Jacobian matrix of the forward,
        i.e. the intrinsic (untransformed) problem, this function returns the
        Jacobian of the (transformed) inverse problem, i.e. taking model and
        data transformations into account by using their inner derivatives.

        Parameters
        ----------
        self : pg.Inversion
            inversion instance with model, response and fop.jacobian
        error_weighted : bool
            add error weighting according to data transform
        numpy_matrix : bool
            return numpy matrix instead of MultLeftRightMatrix
        """
        tData = self.dataTrans.deriv(self.response)
        tModel = 1 / self.modelTrans.deriv(self.model)
        if error_weighted:
            tData *= self.dataTrans.error(self.response, self.errorVals)
        if numpy_matrix:
            J = self.fop.jacobian()
            if isinstance(J, pg.SparseMapMatrix):
                J = asDense(self.fop.jacobian())

            return np.reshape(tData, [-1, 1]) * \
                pg.utils.gmat2numpy(J) * \
                np.reshape(tModel, [1, -1])
        else:
            return pg.matrix.MultLeftRightMatrix(self.fop.jacobian(),
                                                tData, tModel)

    def residual(self):
        """ Residual vector (data-reponse)/error using data transform.
        """
        return (self.dataTrans(self.dataVals) - self.dataTrans(self.response)) / \
            self.dataTrans.error(self.response, self.errorVals)

    def dataGradient(self):
        """ Data gradient from jacobian and residual, i.e. J^T * dData.
        """
        return -self.jacobianMatrix(error_weighted=True).transMult(self.residual())

    def modelGradient(self):
        """ Model gradient, i.e. C^T * C * (m - m0).
        """
        # self.inv.checkConstraints() # not necessary?
        C = pg.matrix.MultLeftMatrix(self.fop.constraints(),
                                     self.inv.cWeight())
        return C.transMult(C.mult(self.modelTrans(self.model)))

    def gradient(self):
        """ Gradient of the objective function.
        """
        return self.dataGradient() + self.modelGradient() * self.lam

# END OF REMOVAL upon pg 1.6

Inversion = ClassicInversion  # pg<1.6
# Inversion = GaussNewtonInversion  # pg>=1.6
# maybe even an inversion chooser


class MarquardtInversion(Inversion):
    """ Marquardt scheme, i.e. local damping with decreasing strength.
    """

    def __init__(self, fop=None, **kwargs):
        super().__init__(fop, **kwargs)
        self.stopAtChi1 = False
        self.inv.setLocalRegularization(True)
        self.inv.setLambdaFactor(0.8)
        self.inv.setDeltaPhiAbortPercent(0.5)

    def run(self, dataVals, errorVals=None, **kwargs):
        r"""Run inversion with given data and error vectors.

        Parameters
        ----------
        dataVals : iterable
            data vector
        errorVals : iterable
            error vector (relative errors), can also be computed from
        absoluteError : float | iterable
            absolute error in units of dataVals
        relativeError : float | iterable
            relative error related to dataVals
        **kwargs:
            Forwarded to the parent class.
            See: :py:mod:`pygimli.modelling.Inversion`
        """
        if errorVals is None:  # use absoluteError and/or relativeError instead
            if not "absoluteError" in kwargs and \
                not "relativeError" in kwargs:
                errorVals = 0.0
            else:
                absErr = kwargs.pop("absoluteError", 0)
                relErr = kwargs.pop("relativeError",
                                    0.01 if np.allclose(absErr, 0, atol=0) else 0)
                errorVals = pg.abs(absErr / dataVals) + relErr

        self.fop.regionManager().setConstraintType(0)
        self.fop.setRegionProperties('*', cType=0)

        self.model = super().run(dataVals, errorVals, **kwargs)
        return self.model


class Block1DInversion(MarquardtInversion):
    """Inversion of layered models (including layer thickness).

    Attributes
    ----------
    nLayers : int
    """

    def __init__(self, fop=None, **kwargs):
        # pg.warn("move this to the manager")
        super(Block1DInversion, self).__init__(fop=fop, **kwargs)

    def setForwardOperator(self, fop):
        """Set forward operator."""
        if not isinstance(fop, pg.frameworks.Block1DModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.Block1DModelling but is of type:', fop)

        return super(Block1DInversion, self).setForwardOperator(fop)

    def fixLayers(self, fixLayers):
        """Fix layer thicknesses.

        Parameters
        ----------
        fixLayers : bool | [float]
            Fix all layers to the last value or set the fix layer
            thickness for all layers
        """
        if fixLayers is False:
            self.fop.setRegionProperties(0, modelControl=1.0)
        elif fixLayers is not None:
            # how do we fix values without modelControl?
            # maybe set the region to be fixed here
            self.fop.setRegionProperties(0, modelControl=1e6)
            if hasattr(fixLayers, '__iter__'):
                if len(fixLayers) != self.fop.nLayers:
                    print("fixLayers:", fixLayers)
                    pg.error("fixlayers needs to have a length of nLayers-1=" +
                             str(self.fop.nLayers - 1))
                self.fop.setRegionProperties(0, startModel=fixLayers)

    def setLayerLimits(self, limits):
        """Set min and max layer thickness.

        Parameters
        ----------
        limits : False | [min, max]
        """
        if limits is False:
            self.fop.setRegionProperties(0, limits=[0.0, 0.0], trans='log')
        else:
            self.fop.setRegionProperties(0, limits=limits, trans='log')

    def setParaLimits(self, limits):
        """Set the limits for each parameter region."""
        for i in range(1, 1 + self.fop.nPara):
            if self.fop.nPara == 1:
                self.fop.setRegionProperties(i, limits=limits, trans='log')
            else:
                self.fop.setRegionProperties(i, limits=limits[i - 1],
                                             trans='log')

    def run(self, dataVals, errorVals,
            nLayers=None, fixLayers=None, layerLimits=None, paraLimits=None,
            **kwargs):
        r"""Run inversion with given data and error vectors.

        Parameters
        ----------
        nLayers : int [4]
            Number of layers.
        fixLayers : bool | [thicknesses]
            See: :py:mod:`pygimli.modelling.Block1DInversion.fixLayers`
            For fixLayers=None, preset or defaults are uses.
        layerLimits : [min, max]
            Limits the thickness off all layers.
            For layerLimits=None, preset or defaults are uses.
        paraLimits : [min, max] | [[min, max],...]
            Limits the range of the model parameter. If you have multiple
            parameters you can set them with a list of limits.

        **kwargs:
            Forwarded to the parent class.
            See: :py:mod:`pygimli.modelling.MarquardtInversion`
        """
        if nLayers is None and "startModel" in kwargs:
            nLayers = (len(kwargs["startModel"]) + 1) // (self.fop.nPara + 1)
        if nLayers is not None:
            self.fop.nLayers = nLayers

        if layerLimits is not None:
            self.setLayerLimits(layerLimits)

        if fixLayers is not None:
            self.fixLayers(fixLayers)

        if paraLimits is not None:
            self.setParaLimits(paraLimits)

        self.model = super(Block1DInversion, self).run(dataVals,
                                                       errorVals, **kwargs)
        return self.model


class MeshInversion(Inversion):
    """Mesh-based inversion.

    ** UNUSED ** TO BE REMOVED or reactivated?
    Attributes
    ----------

    zWeight : float
        relative vertical weight
    """

    def __init__(self, fop=None, **kwargs):
        pg.critical('Obsolete .. to be removed.')
        super(MeshInversion, self).__init__(fop=fop, **kwargs)
        self._zWeight = 1.0

    def setForwardOperator(self, fop):
        """Set forward operator."""
        if not isinstance(fop, pg.frameworks.MeshModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.MeshModelling but is of type:', fop)

        return super(MeshInversion, self).setForwardOperator(fop)

    def run(self, dataVals, errorVals, mesh=None, zWeight=None, **kwargs):
        """Run inversion with given data and error values."""
        if mesh is not None:
            self.fop.setMesh(mesh)

        pg._r('ami here')
        # maybe move this to the fop
        if zWeight is None:
            zWeight = self._zWeight
        self.fop.setRegionProperties('*', zWeight=zWeight)
        # maybe move this to the fop

        pg.debug('run with: ', self.fop.regionProperties())
        # more mesh-related inversion attributes to set?

        # ensure the mesh is generated
        self.fop.mesh()  # not a nice way to ensure something

        self.model = super(MeshInversion, self).run(dataVals,
                                                    errorVals, **kwargs)

        return self.model


class PetroInversion(Inversion):
    """Inversion of a petrophysically related property instead of intrinsic."""

    def __init__(self, petro, fop=None, **kwargs):
        """Initialize inversion.

        Parameters
        ----------
        petro : transformation object (pg.trans)
            petrophysical transformation
        fop : forward operator (pg.Modelling)
            underlying forward operator to be combined with petrophysics
        """
        pg.critical('Obsolete .. to be removed.')
        if fop is not None:
            if not isinstance(fop, pg.frameworks.PetroModelling):
                fop = pg.frameworks.PetroModelling(fop, petro)

        super(PetroInversion, self).__init__(fop=fop, **kwargs)

    def setForwardOperator(self, fop):
        """Set forward operator."""
        if not isinstance(fop, pg.frameworks.PetroModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.PetroModelling but is of type:', fop)

        return super(PetroInversion, self).setForwardOperator(fop)

    def run(self, dataVals, errorVals, **kwargs):
        """Run inversion with given data and error vectors."""
        if 'limits' in kwargs:
            limits = kwargs.pop('limits')

            if len(self.fop.regionManager().regionIdxs()) > 1:
                pg.critical('implement')
            else:
                self.fop.setRegionProperties('*', limits=limits)

        # ensure the mesh is there
        self.fop.mesh()
        return super(PetroInversion, self).run(dataVals, errorVals, **kwargs)


class LCInversion(Inversion):
    """Quasi-2D Laterally constrained inversion (LCI) framework."""

    def __init__(self, fop=None, **kwargs):

        if fop is not None:
            f = pg.frameworks.LCModelling(fop, **kwargs)

        super(LCInversion, self).__init__(f, **kwargs)
        self.dataTrans = pg.trans.TransLog()
        # self.setDeltaChiStop(0.1)

    def prepare(self, dataVals, errorVals, nLayers=4, **kwargs):
        """Prepare inversion with given data and error vectors."""
        dataVec = pg.RVector()
        for d in dataVals:
            dataVec = pg.cat(dataVec, d)

        errVec = pg.RVector()
        for e in errorVals:
            errVec = pg.cat(errVec, e)

        self.fop.initJacobian(dataVals=dataVals, nLayers=nLayers,
                              nPar=kwargs.pop('nPar', None))

        # self.fop.initJacobian resets prior set startmodels
        if self._startModel is not None:
            self.fop.setStartModel(self._startModel)

        rC = self.fop.regionManager().regionCount()

        if kwargs.pop('disableLCI', False):
            self.inv.setMarquardtScheme(0.7)
            # self.inv.setLocalRegularization(True)
            for r in self.fop.regionManager().regionIdxs():
                self.fop.setRegionProperties(r, cType=0)
        else:
            # self.inv.stopAtChi1(False)
            cType = kwargs.pop('cType', None)
            if cType is None:
                cType = [1] * rC

            zWeight = kwargs.pop('zWeight', None)
            if zWeight is None:
                zWeight = [0.0] * rC

            self.fop.setRegionProperties('*',
                                         cType=cType,
                                         zWeight=zWeight,
                                         **kwargs)
            self.inv.setReferenceModel(self.fop.startModel())

        return dataVec, errVec

    def run(self, dataVals, errorVals, nLayers=4, **kwargs):
        """Run inversion with given data and error vectors."""
        lam = kwargs.pop('lam', 20)
        dataVec, errVec = self.prepare(dataVals, errorVals, nLayers, **kwargs)
        print('#'*50)
        print(kwargs)
        print('#'*50)
        return super(LCInversion, self).run(dataVec, errVec, lam=lam, **kwargs)
