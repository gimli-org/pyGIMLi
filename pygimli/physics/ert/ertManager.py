"""Method Manager for Electrical Resistivity Tomography (ERT)."""
from pathlib import Path
import numpy as np
import pygimli as pg
from pygimli.frameworks import MeshMethodManager
from .ertModelling import ERTModelling, ERTModellingReference
from .ert import createInversionMesh, estimateError, simulate
from pygimli.utils import getSavePath


class ERTManager(MeshMethodManager):
    """ERT Manager.

    Method Manager for Electrical Resistivity Tomography (ERT)

    Todo
    ----
        * 3d
        * 3dtopo
        * complex on/off
        * closed geometry
        * transdim
        * singularity removal
        * ERT specific inversion options:
            * ...
    """

    def __init__(self, data=None, **kwargs):
        """Create ERT Manager instance.

        Parameters
        ----------
        data: :gimliapi:`GIMLI::DataContainerERT` | str
            You can initialize the Manager with data or give them a dataset
            when calling the inversion.

        Other Parameters
        ----------------
        * useBert: bool [True]
            Use Bert forward operator instead of the reference implementation.
        * sr: bool [True]
            Calculate with singularity removal technique.
            Recommended but needs the primary potential.
            For flat earth cases the primary potential will be calculated
            analytical. For domains with topography the primary potential
            will be calculated numerical using a p2 refined mesh or
            you provide primary potentials with setPrimPot.
        """
        self.useBert = kwargs.pop('useBert', True)
        self.sr = kwargs.pop('sr', True)

        super().__init__(data=data, **kwargs)
        self.inv.dataTrans = pg.trans.TransLogLU()

    def setSingularityRemoval(self, sr=True):
        """Turn singularity removal on or off."""
        self.reinitForwardOperator(sr=sr)

    def createForwardOperator(self, **kwargs):
        """Create and choose forward operator."""
        verbose = kwargs.pop('verbose', False)
        self.useBert = kwargs.pop('useBert', self.useBert)
        self.sr = kwargs.pop('sr', self.sr)
        if self.useBert:
            pg.verbose('Create ERTModelling FOP')
            fop = ERTModelling(sr=self.sr, verbose=verbose)
        else:
            pg.verbose('Create ERTModellingReference FOP')
            fop = ERTModellingReference(**kwargs)

        return fop

    def load(self, fileName):
        """Load ERT data.

        Forwarded to :py:mod:`pygimli.physics.ert.load`

        Parameters
        ----------
        fileName: str
            Filename for the data.

        Returns
        -------
        data: :gimliapi:`GIMLI::DataContainerERT`
        """
        self.data = pg.physics.ert.load(fileName)
        return self.data

    def createMesh(self, data=None, **kwargs):
        """Create default inversion mesh.

        Forwarded to :py:mod:`pygimli.physics.ert.createInversionMesh`
        """
        data = data or self.data

        if data is None:
            pg.critical('Please provide a data file for mesh generation')

        mesh = createInversionMesh(data, **kwargs)
        self.setMesh(mesh)
        return mesh

    # def simulate(self, mesh, scheme, res, **kwargs):
    def simulate(self, *args, **kwargs):
        """Simulate an ERT measurement.

        Perform the forward task for a given mesh, resistivity distribution &
        measuring scheme and return data (apparent resistivity) or potentials.

        For complex resistivity, the apparent resistivities is complex as well.

        The forward operator itself only calculates potential values for the
        electrodes in the given data scheme.
        To calculate apparent resistivities, geometric factors (k) are needed.
        If there are no values k in the DataContainerERT scheme, the function
        tries to calculate them, either analytically or numerically by using a
        p2-refined version of the given mesh.

        TODO
        ----
        * 2D + Complex + SR

        Args
        ----
        mesh : :gimliapi:`GIMLI::Mesh`
            2D or 3D Mesh to calculate for.

        res : float, array(mesh.cellCount()) | array(N, mesh.cellCount()) |
              list
            Resistivity distribution for the given mesh cells can be:
            - float for homogeneous resistivity (e.g. 1.0)
            - single array of length mesh.cellCount()
            - matrix of N resistivity distributions of length mesh.cellCount()
            - resistivity map [[regionMarker0, res0], [regionMarker0, res1],...]

        scheme : :gimliapi:`GIMLI::DataContainerERT`
            Data measurement scheme.

        Keyword Args
        ------------
        verbose: bool[False]
            Be verbose. Will override class settings.
        calcOnly: bool [False]
            Use fop.calculate instead of fop.response. Useful if you want
            to force the calculation of impedances for homogeneous models.
            No noise handling. Solution is put as token 'u' in the returned
            DataContainerERT.
        noiseLevel: float [0.0]
            add normally distributed noise based on
            scheme['err'] or on noiseLevel if error>0 is not contained
        noiseAbs: float [0.0]
            Absolute voltage error in V
        returnArray: bool [False]
            Returns an array of apparent resistivities instead of
            a DataContainerERT
        returnFields: bool [False]
            Returns a matrix of all potential values (per mesh nodes)
            for each injection electrodes.

        Returns
        -------
        DataContainerERT | array(data.size()) | array(N, data.size()) |
        array(N, mesh.nodeCount()):
            Data container with resulting apparent resistivity data and
            errors (if noiseLevel or noiseAbs is set).
            Optional returns a Matrix of rhoa values
            (for returnArray==True forces noiseLevel=0).
            In case of a complex valued resistivity model, phase values are
            returned in the DataContainerERT (see example below), or as an
            additionally returned array.

        Examples
        --------
        # >>> from pygimli.physics import ert
        # >>> import pygimli as pg
        # >>> import pygimli.meshtools as mt
        # >>> world = mt.createWorld(start=[-50, 0], end=[50, -50],
        # ...                        layers=[-1, -5], worldMarker=True)
        # >>> scheme = ert.createData(
        # ...     elecs=pg.utils.grange(start=-10, end=10, n=21),
        # ...     schemeName='dd')
        # >>> for pos in scheme.sensorPositions():
        # ...     _= world.createNode(pos)
        # ...     _= world.createNode(pos + [0.0, -0.1])
        # >>> mesh = mt.createMesh(world, quality=34)
        # >>> rhomap = [
        # ...     [1, 100. + 0j],
        # ...     [2, 50. + 0j],
        # ...     [3, 10.+ 0j],
        # ... ]
        # >>> data = ert.simulate(mesh, res=rhomap, scheme=scheme, verbose=1)
        # >>> rhoa = data.get('rhoa').array()
        # >>> phia = data.get('phia').array()
        """
        kwargs.setdefault("scheme", self.data)  # to use existing data
        kwargs.setdefault("sr", self.sr)
        return simulate(*args, **kwargs)

    def checkData(self, data=None):
        """Return data from container.

        THINKABOUT: Data will be changed, or should the manager keep a copy?
        """
        data = data or pg.DataContainerERT(self.data)
        # topo = min(pg.z(data)) != max(pg.z(data))  # not used!!

        if isinstance(data, pg.DataContainer):
            if not data.allNonZero('k'):
                pg.critical("Data contains no geometric factors data['k'].")
                # numeric = min(pg.z(data)) != max(pg.z(data))
                # data['k'] = createGeometricFactors(data,
                #                                 numerical=topo,
                #                                 #p2=True,
                #                                 verbose=self.verbose)

            if self.fop.complex():
                if not data.haveData('rhoa'):
                    pg.critical('Datacontainer have no "rhoa" values.')
                if not data.haveData('ip'):
                    pg.critical('Datacontainer have no "ip" values.')

                # pg.warn('check sign of phases')
                rhoa = data['rhoa']
                phia = -data['ip']/1000  # 'ip' is defined for neg mrad.
                # we should think about some 'phia' in rad

                return pg.utils.squeezeComplex(pg.utils.toComplex(rhoa, phia))

            else:
                if not data.haveData('rhoa'):

                    if data.allNonZero('r'):
                        pg.info("Creating apparent resistivies from "
                                "impedences rhoa = r * k")
                        data['rhoa'] = data['r'] * data['k']

                    elif data.allNonZero('u') and data.allNonZero('i'):

                        pg.info("Creating apparent resistivies from "
                                "voltage and currrent rhoa = u/i * k")
                        data['rhoa'] = data['u']/data['i'] * data['k']

                    else:
                        pg.critical("Datacontainer have neither: "
                                    "apparent resistivies 'rhoa', "
                                    "or impedances 'r', "
                                    "or voltage 'u' along with current 'i'.")

                if any(data['rhoa'] < 0) and \
                        isinstance(self.inv.dataTrans, pg.trans.TransLog):
                    # print(pg.find(data['rhoa'] < 0))
                    # print(data['rhoa'][data['rhoa'] < 0])
                    pg.warning("Found negative apparent resistivities. "
                               "These can't be processed with logarithmic "
                               "data transformation. You should consider to "
                               "filter them out using "
                               "data.remove(data['rhoa'] < 0).")

                return data['rhoa']

        return data

    def checkErrors(self, err, dataVals):
        """Check (estimate) and return relative error.

        By default we assume 'err' are relative values.
        """
        if isinstance(err, pg.DataContainer):
            rae = None

            if not err.allNonZero('err'):
                pg.warn("Datacontainer have no 'err' values. "
                        "Fallback of 1mV + 3% using "
                        "ERTManager.estimateError(...) ")
                rae = self.estimateError(err, absoluteError=0.001,
                                         relativeError=0.03)
            else:
                rae = err['err']

            if self.fop.complex():
                ipe = None

                if err.haveData('iperr'):
                    _, phi = pg.utils.toPolar(dataVals)
                    # assuming ipErr are absolute dPhi in mrad
                    ipe = err['iperr'] / abs(phi*1000)
                else:
                    pg.warn("Datacontainer have no 'iperr' values. "
                            "Fallback set to 0.01")
                    ipe = np.ones(err.size()) * 0.01

                return pg.cat(rae, ipe)

        return rae  # not set if err is no DataContainer (else missing)

    def estimateError(self, data=None, **kwargs):
        """Estimate error composed of an absolute and a relative part.

        Parameters
        ----------
        absoluteError : float [0.001]
            Absolute data error in Ohm m. Need 'rhoa' values in data.

        relativeError : float [0.03]
            relative error level in %/100

        absoluteUError : float [0.001]
            Absolute potential error in V. Need 'u' values in data. Or
            calculate them from 'rhoa', 'k' and absoluteCurrent if no 'i'
            is given

        absoluteCurrent : float [0.1]
            Current level in A for reconstruction for absolute potential V

        Returns
        -------
        error : Array
        """
        if data is None:  #
            error = estimateError(self.data, **kwargs)
            self.data["err"] = error
        else:  # the old way: better use ert.estimateError directly
            error = estimateError(data, **kwargs)

        return error

    def coverage(self):
        """Coverage vector considering the logarithmic transformation."""
        covTrans = pg.core.coverageDCtrans(self.fop.jacobian(),
                                           1.0 / self.inv.response,
                                           1.0 / self.inv.model)

        paramSizes = np.zeros(len(self.inv.model))
        for c in self.fop.paraDomain.cells():
            paramSizes[c.marker()] += c.size()

        wCovLog = np.log10(covTrans / paramSizes)
        return wCovLog[self.fop.paraDomain.cellMarkers()]

    def standardizedCoverage(self, threshold=0.01):
        """Return standardized coverage vector (0|1) using thresholding."""
        return 1.0*(self.coverage() > threshold)

    def showMisfit(self, errorWeighted=False, **kwargs):
        """Show relative or error-weighted data misfit."""
        if errorWeighted:
            misfit = (np.log(self.data["rhoa"]) - np.log(self.inv.response)) /\
                self.data["err"]
            kwargs.setdefault("label", "error-weighted misfit")
        else:
            misfit = - self.inv.response / self.data["rhoa"] * 100 + 100
            kwargs.setdefault("label", "relative misfit (%)")

        kwargs.setdefault("cMax", np.max(np.abs(misfit)))
        kwargs.setdefault("cMin", -kwargs["cMax"])
        kwargs.setdefault("cMap", "bwr")
        kwargs.setdefault("logScale", False)

        self.showData(vals=misfit, **kwargs)

    def showModel(self, model=None, ax=None, elecs=True, **kwargs):
        """Show the last inversion result.

        Parameters
        ----------
        model : iterable [None]
            Model vector to be drawn. Default is self.model from the last run
        ax : mpl axes
            Axes object to draw into. Create a new if its not given.
        elecs : bool
            Draw electrodes

        **kwargs : dict
            arguments passed to pg.show (cMin, cMax, cMap, logScale)

        Returns
        -------
        ax, cbar
        """
        if model is None:
            model = self.model

        if ax is None:
            _, ax = pg.plt.subplots()

        kwargs.setdefault("coverage", self.coverage())
        color = kwargs.pop("color", "magenta")
        ax, cBar = self.fop.drawModel(ax, model, **kwargs)
        if elecs:
            if isinstance(elecs, str):
                color = elecs

            pg.viewer.mpl.drawSensors(ax, self.fop.data.sensors(), color=color)

        return ax, cBar

    def saveResult(self, folder=None, size=(16, 10), **kwargs):
        """Save all results in the specified folder.

        Saved items are:
            Inverted profile
            Resistivity vector
            Coverage vector
            Standardized coverage vector
            Mesh (bms and vtk with results)
        """
        subfolder = self.__class__.__name__
        path = Path(getSavePath(folder, subfolder))
        pg.info(f'Saving inversion results to: {path}')
        np.savetxt(path / 'resistivity.vector', self.model)
        np.savetxt(path / 'resistivity-cov.vector', self.coverage())
        np.savetxt(path / 'resistivity-scov.vector',
                   self.standardizedCoverage())

        self.mesh.save(str(path / 'mesh.bms'))

        m = pg.Mesh(self.paraDomain)
        m['Resistivity'] = self.model
        m['Resistivity (log10)'] = np.log10(m['Resistivity'])
        m['Coverage'] = self.coverage()
        m['S_Coverage'] = self.standardizedCoverage()
        nM = m.cellCount()
        for k, v in kwargs.items():
            if hasattr(v, "__iter__") and len(v) == nM:
                m[k] = v

        m.exportVTK(str(path / 'resistivity.vtk'))
        m.saveBinaryV2(str(path / 'resistivity-pd.bms'))
        self.fop.mesh().save(str(path / 'resistivity-mesh'))

        np.savetxt(path / 'response.vector', self.inv.response)
        residual = self.inv.residual()  # includes error-(re)weighting
        np.savetxt(path / 'residual.vector', residual)
        data = self.fop.data.copy()
        data['response'] = self.inv.response
        data['residual'] = residual
        data.save(str(path / 'data.dat'),
                  'a b m n rhoa k err ip iperr response residual')

        if self.paraDomain.dim() == 2:
            fig, ax = pg.plt.subplots(figsize=size)
            self.showModel(ax=ax, **kwargs)
            fig.savefig(path / 'resistivity.pdf', bbox_inches="tight")
            return path, fig, ax

        return path


if __name__ == "__main__":
    pass
