"""Framwork for structural coupled cooperative inversion (SCCI).

Cooperative of several geophysical methods by structural coupling of
model roughnesses. The methodology goes back to the following papers:

* Günther & Rücker (2006): https://doi.org/10.4133/1.2923578
* Günther et al. (2010): https://doi.org/10.4133/1.3445447
* Hellman et al. (2017): http://dx.doi.org/10.1016/j.jappgeo.2017.06.008
* Ronczka et al. (2017): http://dx.doi.org/10.5194/se-8-671-2017

and was reimplemented in the projec
COMET - COupled Magnetic resonance Electrical resistivity Tomography
by Nico Skibbe et al. (2020), and applied in the papers

* Skibbe et al. (2018): http://dx.doi.org/10.1190/geo2018-0046.1
* Skibbe et al. (2020): http://dx.doi.org/10.1190/geo2019-0484.1
* Skibbe et al. (2021): http://dx.doi.org/10.1190/geo2020-0593.1

It uses managers holding data and inversion instances.
"""
import os
import numpy as np
import pygimli as pg


class SCCI(object):
    """Structurally coupled cooperative inversion class."""

    def __init__(self, managers=[], **kwargs):
        """Structurally coupled cooperative inversion framework."""
        # scci parameters
        self.a = None
        self.b = None
        self.c = None
        self.cmin = None
        self.cmax = None
        self.setCouplingPara()

        # init managers list
        self.managers = managers
        # I think we don't need the managers unless we plot etc.

        self.resultdir = None  # "scci_results"

        self.names = kwargs.pop("names",
                                ["M{0}".format(i+1) for i in range(10)])
        # init constraint weights

    def setCouplingPara(self, a=0.1, b=0.1, c=1.0, cmin=None, cmax=1.0):
        """Set coupling parameter."""
        self.a = a
        self.b = b
        self.c = c
        self.cmin = cmin or b
        self.cmax = cmax

    @property
    def invs(self):
        """Return inversions."""
        return self._gather('inv', raise_error=True)

    @property
    def fops(self):
        """Return forward operators."""
        return self._gather('fop', raise_error=True)

    @property
    def npara(self):
        """Number of total model parameters."""
        npara = 0
        paras = self._gather('_numPara', default=1)
        npara += np.sum(paras)
        return npara

    @property
    def nparas(self):
        """Return model parameters for the individual methods."""
        paras = self._gather('_numPara', default=1)
        return paras

    @property
    def roughnesses(self):
        """Return roughness vectors from individual models."""
        output = []
        # do for all managers == inversion instances
        for im, man in enumerate(self.managers):
            model = man.inv.model

            if man.fop.constraints().cols() != model.size():
                man.fop.createConstraints()

            roughness = man.inv.inv.pureRoughness(model)

            seg = int(len(roughness)/self.nparas[im])
            # do for all parameters
            for i in range(self.nparas[im]):
                output.append(np.array(roughness[i * seg:i * seg + seg]))

        return output

    def singleCWeights(self):
        """Constraint weight vectors of the individual models."""
        cweights = []
        for rough in self.roughnesses:
            cweight = roughness2CWeight(
                rough, a=self.a, b=self.b, c=self.c, Min=0,  # self.cmin,
                Max=self.cmax)
            pg.debug([np.min(cweight), np.max(cweight)])
            cweights.append(cweight)

        return cweights

    def updateConstraintWeights(self):
        """Set new constraint weights based on the model roughnesses."""
        pg.debug('SCCI: updateConstraintWeights()')
        # write down single c weights
        single_weights = np.array(self.singleCWeights())
        pg.debug('single c weights before updateConstraintWeights: {}'
                 .format(single_weights))

        new_cweight = np.ones_like(single_weights)

        # each new c weight consists of the combination of cweights of all
        # other parameters
        for ipar in range(self.npara):
            all_others = np.arange(self.npara) != ipar
            pg.debug('SCCI: all others: {}'.format(all_others))
            new_cweight[all_others] *= single_weights[ipar]
            # new_cweight *= single_weights[ipar]

        # cut extreme values according to cmin, cmax
        new_cweight = np.minimum(new_cweight, self.cmax)
        new_cweight = np.maximum(new_cweight, self.cmin)

        pg.debug('min/max new c weight: {}/{}'
                 .format(np.min(new_cweight), np.max(new_cweight)))

        # set c weights
        total = 0
        matrices = []
        for iinv, inv in enumerate(self.invs):
            weight = []
            for j in range(self.nparas[iinv]):
                weight.extend(new_cweight[total])
                total += 1

            pg.debug('SCCI: manager {}, weights {}'.format(
                iinv, np.shape(weight)))

            inv.inv.setCWeight(weight)

            matrices.append(inv.fop.constraints())

        return new_cweight, matrices

    def run(self, save=False, **kwargs):
        """Run coupled inversion."""
        maxIter = kwargs.pop("maxIter", 8)
        if self.resultdir is None:
            self.resultdir = "result_a{}_b{}".format(self.a, self.b).replace(
                ".", "_")

        if save and not os.path.exists(self.resultdir):
            os.makedirs(self.resultdir)

        for i in range(maxIter):
            print("Coupled inversion {0}".format(i+1))
            self.updateConstraintWeights()
            for i, inv in enumerate(self.invs):
                basename = self.resultdir + "/" + self.names[i]
                if save:
                    np.save(basename+'_cWeight_{}.npy'.format(i + 1),
                            inv.inv.cWeight())

                inv.inv.oneStep()
                if save:
                    np.save(basename + '_response_{}.npy'.format(i + 1),
                            inv.response)
                    np.save(basename + '_model_{}.npy'.format(i + 1),
                            inv.model)

    def _gather(self, attribute, raise_error=False, default=None):
        """Gather variables from underlaying managers for convenience."""
        ret = []
        for mi, manager in enumerate(self.managers):
            if hasattr(manager, attribute):
                ret.append(getattr(manager, attribute))
            else:
                if raise_error:
                    raise AttributeError(
                        'Manager {} of type {} has no attribute "{}"'
                        .format(mi, type(manager), attribute))
                else:
                    ret.append(default)
        return ret

    def _gather_from_cpp(self, call, raise_error=False, default=None):
        """Gather variables from given callable from underlaying managers."""
        ret = []
        for mi, manager in enumerate(self.managers):
            if hasattr(manager, call):
                to_call = getattr(manager, call)
                if not callable(to_call):
                    if raise_error:
                        raise TypeError('{} object of manager {} of type {} '
                                        'is not callable'.format(
                                            type(to_call), mi, type(manager)))
                    else:
                        ret.append(default)
                else:
                    ret.append(to_call())
            else:
                if raise_error:
                    raise AttributeError(
                        'Manager {} of type {} has no function "{}"'
                        .format(mi, type(manager), call))
                else:
                    ret.append(default)
        return ret


def roughness2CWeight(vec, a=0.1, b=0.1, c=1.0, Max=1.0, Min=0.2):
    """Structural constraint weight as function of roughness.

    See Günther et al. (2010, SAGEEP extended abstract) (case a>0) for details.
    """
    avec = np.absolute(vec)
    cfun = (a / (avec + a) + b)**c
    # confine between min and max
    # cfun = (cfun - min(cfun)) / (max(cfun) - min(cfun)) * (Max-Min) + Min
    cfun = np.minimum(cfun, Max)
    cfun = np.maximum(cfun, Min)
    return cfun


def getNormalizedRoughness(mesh, model, trans='log', tmin=None, tmax=None):
    """Return normalized roughness of the model given mesh and transformation.

    The roughness is defined on each cell boundary with connection to another
    cell (so not for cells at the boundary of the mesh), so the number of
    values is not number of nodes nor number of cells or boundaries and
    strongly depends on the mesh.
    """
    # Step 1: init of mesh, fop and inv
    # we just need the methods of the inv and fop instances, so no input needed
    mesh.createNeighbourInfos()

    fop = pg.core.ModellingBase()
    fop.setMesh(mesh)

    # pseudo data, fop, verbose, debug
    inv = pg.core.Inversion([1, 2, 3], fop, True, False)

    # Step 2: take care of transformation
    if trans.lower() == 'log':
        if tmin is not None and tmax is not None:
            transmodel = pg.core.TransLogLU(tmin, tmax)
        else:
            transmodel = pg.core.TransLogLU()

    elif trans.lower() == 'cot':
        if tmin is None:
            tmin = 0

        if tmax is None:
            tmax = 0.7

        transmodel = pg.core.TransCotLU(tmin, tmax)

    else:
        raise Exception(
            'Transformation not implemented, expect "log" or "cot", not "{}".'
            .format(trans))

    # Step 3: region manager initialization (enables use of createConstraints)
    inv.setTransModel(transmodel)
    fop.regionManager()
    fop.createConstraints()

    # Step 4: with active constraints for the mesh, we simply get the roughness
    roughness = inv.pureRoughness(model)

    # return as numpy array
    return np.asarray(roughness)

