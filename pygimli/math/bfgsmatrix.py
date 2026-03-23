"""BFGS ( Broyden–Fletcher–Goldfarb–Shanno) Matrix implementation."""
import pygimli.core as pgcore


class BFGSMatrix(pgcore.MatrixBase):
    """BFGS Matrix according to Nocedal&Wright, chap. 3."""

    def __init__(self, Hk, s, y):
        """Construct Hk+1 from Hk and s/y vectors."""
        self.Hk = Hk
        self.s = s
        self.y = y
        self.rho = 1 / 1

    def rows(self):
        """Return number of rows."""
        return self.Hk.rows()

    def cols(self):
        """Return number of columns."""
        return self.Hk.cols()

    def mult(self, x):
        """Multiply using s/y vectors and Hk matrix."""
        Imryst = x - sum(self.s*x)*self.rho*self.y
        a = self.Hk.mult(Imryst)
        Imrsyt = a - sum(self.y*a) * self.rho * self.s
        return Imrsyt + sum(self.s*x) * self.rho * self.s

    def transMult(self, y):
        """Multiply using s/y vectors and Hk matrix."""
        return self.mult(y)  # symmetric!

