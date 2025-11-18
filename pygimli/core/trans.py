# -*- coding: utf-8 -*-
"""Some specialization to the trans functions."""

from .core import pgcore
import numpy as np

__TransCumulative_addForGC__ = pgcore.RTransCumulative.add


def __TransCumulative_addForGC_MP__(self, T, *args):
    """Don't use directly.

    Monkey patch to keep the GC happy until redesign.
    """
    if not hasattr(self, '__trans__'):
        self.__trans__ = []
    self.__trans__.append(T)

    if len(args) == 1:
        # be sure avoid auto conversion from int to IndexArray
        if isinstance(args[0], int):
            return __TransCumulative_addForGC__(self, T, size=args[0])
    return __TransCumulative_addForGC__(self, T, *args)


pgcore.RTransCumulative.add = __TransCumulative_addForGC_MP__

def __RTransCumulative_str(self):
    """String representation."""
    out = "Cumulative Data transformation:"
    for i in range(self.size()):
        itr = self.at(i).__repr__()
        out += f"\n{i}, {itr}"

    return out

pgcore.RTransCumulative.__repr__ = __RTransCumulative_str


def __RTransLog_str(self):
    out = "Logarithmic transform"
    lB = self.lowerBound()
    if lB != 0.0:
        out += f", lower bound {lB}"
    return out

pgcore.RTransLog.__repr__ = __RTransLog_str

def __RTransLogLU_str(self):
    """String representation."""
    out = "Logarithmic LU transform"
    out += f", lower bound {self.lowerBound()}"
    out += f", upper bound {self.upperBound()}"
    return out

pgcore.RTransLogLU.__repr__ = __RTransLogLU_str

def __RTransCotLU_str(self):
    """String representation."""
    out = "Cotangens LU transform"
    # bounds not available (change in C++)
    # out += f", lower bound {self.lowerBound()}"
    # out += f", upper bound {self.upperBound()}"
    return out

pgcore.RTransCotLU.__repr__ = __RTransCotLU_str

def __RTrans_str(self):
    """String representation."""
    return "Identity transform"

pgcore.RTrans.__repr__ = __RTrans_str


# Aliases
Trans = pgcore.RTrans
TransLinear = pgcore.RTransLinear
TransLin = pgcore.RTransLin
TransPower = pgcore.RTransPower
TransLog = pgcore.RTransLog
TransLogLU = pgcore.RTransLogLU
TransCotLU = pgcore.RTransCotLU
TransCumulative = pgcore.RTransCumulative


class TransSymLog(pgcore.RTrans):
    """Transformation using a bilogarithmic scaling."""

    def __init__(self, tol=1e-12):
        """Forward transformation."""
        super().__init__()
        self.tol = tol

    def __repr__(self):
        """String representation."""
        return f"Symlog transformation with threshold {self.tol}"

    def trans(self, x):
        """Forward transformation."""
        return pgcore.log(1 + np.abs(x) / self.tol) * np.sign(x)

    def invTrans(self, y):
        """Inverse transformation."""
        return (pgcore.exp(np.abs(y)) - 1.) * self.tol * np.sign(y)

    def deriv(self, x):
        """Derivative of the transformation."""
        return 1. / (np.abs(x) / self.tol + 1) / self.tol


def str2Trans(tr:str):
    """Convert string to transformation.

    Convention
    ----------
    lin : linear
    log : logarithmic
    logL : log with lower bound
    logL-U : log with lower and upper bound
    cotL-U : cotangent with lower and upper bound
    symlogT : symlog with T as lin-threshold
    """
    low = tr.lower()
    if low.startswith("lin"):
        return Trans()
    elif low.startswith("log"):
        if "-" in low:  # lower/upper
            sp = low[3:].split("-")
            return TransLogLU(float(sp[0]), float(sp[1]))
        elif len(low) > 3: # lower
            return TransLog(float(low[3:]))
        else: # just log
            return TransLog()
    elif low.startswith("cot") and "-" in low:
        sp = low[3:].split("-")
        return TransCotLU(float(sp[0]), float(sp[1]))
    elif low.startswith("symlog"):
        return TransSymLog(float(low[6:]))
    else:  # check for LU values, e.g. "Log1-1000" or "Cot0-1"
        raise KeyError("Transformation string unknown!")


