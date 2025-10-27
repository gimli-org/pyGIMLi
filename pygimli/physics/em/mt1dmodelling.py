"""1D Magnetotelluric modelling classes."""
import numpy as np
import pygimli as pg
from pygimli.frameworks import Block1DModelling, MeshModelling


class MT1dBlockModelling(Block1DModelling):
    """MT 1d few-layer modelling with resistivity & thickness."""

    def __init__(self, T=None, nLayers=3, verbose=True, **kwargs):
        """Set up class with periods."""
        super().__init__(verbose=verbose, **kwargs)
        self.T = T
        self.nLayers = nLayers
        self.fwd = pg.core.MT1dModelling(self.T, nLayers, verbose=verbose)

    def response(self, par):
        """Model response."""
        return self.response_mt(par, 0)

    def response_mt(self, par, i=0):
        """Multi-threading model response."""
        nLayers = (len(par)+1) // 2
        fop = pg.core.MT1dModelling(self.T, int(nLayers), verbose=self.verbose())

        return fop.response(par)

    def createStartModel(self, rhoa):
        """Create starting model."""
        if self.nLayers == 0:
            pg.critical("Model space is not been initialized.")

        # layer thickness properties
        self.setRegionProperties(0, startModel=1000, trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa), trans='log')

        return super().createStartModel()

    def drawData(self, ax, data, **kwargs):
        """Draw MT sounding curves (app. resistivity & phase)."""
        rhoa, phase = np.split(data, 2)
        ax2 = pg.viewer.mpl.createTwinY(ax)
        kwargs.setdefault('marker', 'x')
        color = kwargs.pop('color', 'C0')
        kwargs.setdefault('label', 'rhoa')
        ax.loglog(rhoa, self.T, color=color, **kwargs)
        kwargs['label'] = 'phase'
        ax2.semilogy(np.rad2deg(phase), self.T, color="C1", **kwargs)
        ax.set_ylim(max(self.T), min(self.T))
        ax.set_ylabel(r"$T$ [s]")
        ax.set_xlabel(r"$\rho_a$ [$\Omega$m]")
        ax2.set_xlabel(r"$\phi_a$ [°]")
        ax.legend()
        ax2.legend()
        ax.grid()

    def drawModel(self, ax, model, **kwargs):
        """Draw model as 1D block model."""
        kwargs.setdefault('plot', 'loglog')
        pg.viewer.mpl.drawModel1D(ax=ax, model=model,
                                  xlabel=r'Resistivity ($\Omega$m)', **kwargs)
        ax.set_ylabel('Depth in (m)')
        return ax, None  # should return gci and not ax


class MT1dSmoothModelling(MeshModelling):
    """MT 1d few-layer modelling based on predefined thickness."""

    def __init__(self, T=None, thk=None, verbose=True, **kwargs):
        """Set up class with periods."""
        super().__init__(**kwargs)
        self.T = T
        self.thk = thk
        self.fwd = pg.core.MT1dRhoModelling(pg.Vector(self.T), pg.Vector(self.thk),
                                            verbose=verbose)
        self.mesh_ = pg.meshtools.createMesh1D(len(thk)+1)
        self.setMesh(self.mesh_)

    def response(self, model):
        """Return forward response of model."""
        return self.fwd.response(model)

    def createStartModel(self, dataVals=None):
        """Create starting model."""
        return pg.Vector(len(self.thk)+1, np.median(dataVals))

    def drawModel(self, ax, model, **kwargs):
        """Draw model as 1D multi-layered model."""
        kwargs.setdefault('plot', 'loglog')
        pg.viewer.mpl.drawModel1D(ax=ax, thickness=self.thk, values=model,
                                  xlabel=r'Resistivity ($\Omega$m)', **kwargs)
        ax.set_ylabel('Depth in (m)')
        return ax, None  # should return gci and not ax

MT1dSmoothModelling.drawData = MT1dBlockModelling.drawData
