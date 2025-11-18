"""Tools for (time domain) IP modelling using chargeability."""
import pygimli as pg


class IPSeigelModelling(pg.frameworks.MeshModelling):
    """DC/IP modelling class using an (FD-based) approach."""

    def __init__(self, f, mesh, rho, verbose=False, response=None):
        """Init class with DC forward operator and resistivity vector."""
        super().__init__(mesh=mesh, verbose=verbose)
        # self.setMesh(mesh)
        self.f = f
        self.rhoDC = rho  # DC resistivity
        self.rhoAC = rho * 1  # AC resistivity
        self.drhodm = -self.rhoDC
        if response is None:
            response = f.response(self.rhoDC)

        self.rhoaDC = response
        self.rhoaAC = pg.Vector(self.rhoaDC)  # a copy to have correct size.

        self.dmdrhoa = -1.0 / self.rhoaDC
        self.J = pg.matrix.MultLeftRightMatrix(f.jacobian(), self.dmdrhoa,
                                               self.drhodm)
        self.setJacobian(self.J)
        self.fullJacobian = False

    def response(self, m):
        """Return forward response as function of chargeability model."""
        self.rhoaAC = self.f.response(self.rhoDC / (1. - m))
        if self.fullJacobian:
            self.dmdrhoa = -1.0 / self.rhoaAC

        return pg.abs(1.0 - self.rhoaDC / self.rhoAC)

    # RhoAC = RhoDC * (1-m) => dRhoa / m = dRhoa/dRho * dRho/m = JDC * (-RhoDC)
    def createJacobian(self, model):
        """Create jacobian matrix using unchanged DC jacobian and m model."""
        if self.fullJacobian:
            self.rhoAC = self.rhoDC * (1 - model)
            self.J.right = -self.rhoAC
            self.f.createJacobian(self.rhoAC)


class DCIPSeigelModelling(pg.frameworks.MeshModelling):
    """DC/IP modelling class using an (FD-based) approach."""

    # def __init__(self, f=None, mesh=None, rho=None, resp=None, verbose=False)
    def __init__(self, ERT, verbose=False):
        """Init class with DC forward operator and resistivity vector."""
        super().__init__(verbose=verbose)
        self.setMesh(ERT.paraDomain)
        self.resp = ERT.inv.response
        self.res = ERT.inv.model
        self.J = pg.matrix.MultLeftRightMatrix(ERT.fop.jacobian(),
                                               1./self.resp, self.res)
        # G = pg.utils.gmat2numpy(ERT.fop.jacobian())
        # Glog = np.reshape(1./self.resp, (-1, 1)) * G * self.res
        # self.J = Glog / np.sum(Glog, 1)
        self.setJacobian(self.J)

    def response(self, m):
        """Return forward response as function of chargeability model."""
        return self.J.dot(m)

    def createJacobian(self, model):
        """Create jacobian matrix using unchanged DC jacobian and m model."""
        pass


class DCIPMModelling(pg.frameworks.MeshModelling):
    """DC/IP modelling class using an (FD-based) approach."""

    def __init__(self, f, mesh, rho, verbose=False, response=None):
        """Init class with DC forward operator and resistivity vector."""
        super().__init__(verbose=verbose)
        self.setMesh(mesh)
        self.f = f
        self.rhoDC = rho  # DC resistivity
        self.rhoAC = rho * 1  # AC resistivity
        self.drhodm = -self.rhoDC
        if response is None:
            self.rhoaDC = f.response(self.rhoDC)
        else:
            self.rhoaDC = response

        self.rhoaAC = self.rhoaDC * 1
        self.dmdrhoa = -1.0 / self.rhoaDC
        self.J = pg.matrix.MultLeftRightMatrix(f.jacobian(), self.dmdrhoa,
                                               self.drhodm)
        self.setJacobian(self.J)
        self.fullJacobian = False

    def response(self, m):
        """Return forward response as function of chargeability model."""
        self.rhoaAC = self.f.response(self.rhoDC * (1. - m))
        if self.fullJacobian:
            self.dmdrhoa = -1.0 / self.rhoaAC

        return pg.abs(1.0 - self.rhoaAC / self.rhoaDC)

    # RhoAC = RhoDC * (1-m) => dRhoa / m = dRhoa/dRho * dRho/m = JDC * (-RhoDC)
    def createJacobian(self, model):
        """Create jacobian matrix using unchanged DC jacobian and m model."""
        if self.fullJacobian:
            self.rhoAC = self.rhoDC * (1 - model)
            self.J.right = -self.rhoAC
            self.f.createJacobian(self.rhoAC)


if __name__ == "__main__":
    pass
