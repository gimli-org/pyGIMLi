"""Keep the documentation build headless.

Sphinx-gallery executes the gallery examples either inside the Sphinx process
or, with ``sphinx_gallery_conf["parallel"]``, inside joblib worker processes.
The workers are separate interpreters and only inherit the environment, so a
GUI backend that ends up in ``MPLBACKEND`` would open a window for every
figure -- and examples that end with :func:`pygimli.wait` then block the build
until the window is closed by hand.

:func:`forceHeadless` is used as a sphinx-gallery ``reset_modules`` hook, which
runs before every executed example, and therefore re-applies the headless
settings in whichever process actually runs the example. Set
``GIMLI_DOC_SHOW=1`` to keep the interactive behaviour for debugging an
example.
"""

import os
import sys

__all__ = ["forceHeadless"]


def isHeadless():
    """Return whether figures should be rendered off screen.

    Returns
    -------
    bool
        False if ``GIMLI_DOC_SHOW`` asks for interactive figure windows.
    """
    return os.environ.get("GIMLI_DOC_SHOW", "0") != "1"


def forceHeadless(galleryConf=None, fname=None):
    """Render figures off screen in the current process.

    Used as a sphinx-gallery ``reset_modules`` hook. The arguments are part of
    the hook signature and are not used.

    Parameters
    ----------
    galleryConf : dict
        Sphinx-gallery configuration, unused.
    fname : str
        Name of the example that is about to run, unused.
    """
    if not isHeadless():
        return

    os.environ["MPLBACKEND"] = "Agg"
    os.environ["PYVISTA_OFF_SCREEN"] = "true"

    import matplotlib

    # Switching the backend closes all open figures, so only do it when the
    # process is not headless yet.
    if matplotlib.get_backend().lower() != "agg":
        matplotlib.use("Agg", force=True)

    import matplotlib.pyplot as plt

    plt.ioff()

    # pyvista reads PYVISTA_OFF_SCREEN only at import time, so an already
    # imported pyvista has to be adjusted explicitly.
    pyvista = sys.modules.get("pyvista")
    if pyvista is not None:
        pyvista.OFF_SCREEN = True
