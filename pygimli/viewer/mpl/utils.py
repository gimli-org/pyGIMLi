#!/usr/bin/env python
"""Plotting utilities used throughout the viewer.mpl package."""
import sys
import os
import time
import numpy as np

import pygimli as pg
from pygimli.utils import prettyFloat

__holdAxes__ = True
__lastBackend__ = None


def updateFig(fig, force=True, sleep=.001):
    """For internal use."""
    if not globals()['__holdAxes__']:
        try:
            fig.canvas.draw_idle()
            if force:
                fig.canvas.flush_events()
                time.sleep(sleep)
        except BaseException as e:
            print(fig, e)
            pg.warn("Exception raised", e)


def updateAxes(ax, force=True):
    """For internal use."""
    updateFig(ax.figure, force=force)


def hold(val=True):
    """Set axis plotting state to hold (not done yet)."""
    globals()['__holdAxes__'] = val


def wait(**kwargs):
    """Show figure and wait until it is closed or key is pressed."""
    # plt.pause seems to be broken in mpl:2.1
    # ax.canvas.draw_onIdle()
    import matplotlib.pyplot as plt

    #if len(plt.get_fignums()) > 0:
    #for f in plt.get_fignums():
    #    updateFig(f)
    updateFig(plt.gca())
    kp = kwargs.pop('untilKeyPressed', False)
    if kp:
        plt.waitforbuttonpress(**kwargs)
    else:
        plt.show(**kwargs)


def noShow(on=True):
    """Toggle quiet mode to avoid popping figures.

    Arguments
    ---------
    on: bool[True]
        Set Matplotlib backend to 'agg' and restore old backend if set to False.
    """
    import matplotlib
    if on is True:
        globals()[__lastBackend__] = matplotlib.get_backend()
        matplotlib.use('agg')
    else:
        if globals()[__lastBackend__] is not None:
            matplotlib.use(globals()[__lastBackend__])


__registeredShowPendingFigsAtExit__ = False

def registerShowPendingFigsAtExit():
    """Register closing function to show all figures.

    If called, register a closing function that will ensure all pending MPL
    figures are shown. It is only set on demand by pg.show() since we only
    need it if matplotlib is used.
    """
    global __registeredShowPendingFigsAtExit__
    if not __registeredShowPendingFigsAtExit__:
        import atexit

        #pg._y('register wait on exit')
        ## first  call one empty show to initial QtManager before register
        # onExit to avoid RuntimeError: wrapped C/C++ object MainWindow been deleted
        if 'matplotlib.pyplot' in sys.modules:
            import matplotlib.pyplot as plt
            #pg._g('empty show')
            plt.show()

        @atexit.register
        def waitOnExit():
            """Call on script end to ensure to open all remaining mpl figures."""
            #pg._g('waitonexit')
            # only call it if mpl has been used

            if 'matplotlib.pyplot' in sys.modules:
                import matplotlib.pyplot as plt

                backend = sys.modules['matplotlib.pyplot'].get_backend()
                #backend = plt.get_backend()
                #print('backend', backend)

                if 'inline' not in backend \
                    and ('qt' in backend.lower()
                          or 'wx' in backend.lower()
                          or 'tk' in backend.lower()
                          or 'gtk' in backend.lower()):

                    if len(plt.get_fignums()) > 0:
                        pg.info(f'Showing pending widgets ({backend}) on exit. '
                                 'Close all figures or Ctrl-C to quit the ' \
                                 'program')
                        pg.wait()

    __registeredShowPendingFigsAtExit__ = True


def saveFigure(fig, filename, pdfTrim=False):
    """Save figure as pdf."""
    if '.pdf' in filename:
        filename = filename[0:filename.find('.pdf')]
    fig.savefig(filename + '.pdf', bbox_inches='tight')
    # pdfTrim=1
    if pdfTrim:
        try:
            print("trying pdf2pdfS ... ")
            os.system('pdf2pdfBB ' + filename + '.pdf')
            os.system('pdf2pdfS ' + filename + '.pdf')
        except BaseException as _:
            print("fail local convert. Should be no problem.")


def saveAxes(ax, filename, adjust=False):
    """Save axes as pdf."""
    if adjust:
        adjustWorldAxes(ax)

    updateAxes(ax, force=True)
    saveFigure(ax.figure, filename)


def insertUnitAtNextLastTick(ax, unit, xlabel=True, position=-2):
    """Replace the last-but-one tick label by unit symbol."""
    if xlabel:
        labels = ax.get_xticklabels()
        labels[position] = unit
        ax.set_xticklabels(labels)
    else:
        labels = ax.get_yticklabels()
        labels[-position] = unit
        ax.set_yticklabels(labels)


def adjustWorldAxes(ax, useDepth:bool=True, xl:str='$x$ in m', yl:str=None):
    """Set some common default properties for an axe."""
    if yl is None:
        if useDepth is True:
            ax.set_ylabel('Depth in m')
            renameDepthTicks(ax)
        else:
            ax.set_ylabel('$y$ in m')
    else:
        ax.set_ylabel(yl)

    ax.set_xlabel(xl)

    ax.figure.tight_layout()
    updateAxes(ax)


def renameDepthTicks(ax):
    """Switch signs of depth ticks to be positive."""
    from matplotlib import ticker

    @ticker.FuncFormatter
    def major_formatter(x, pos):
        return prettyFloat(-x) % x

    ax.yaxis.set_major_formatter(major_formatter)
    updateAxes(ax)


def setOutputStyle(dim='w', paperMargin=5, xScale=1.0, yScale=1.0, fontsize=9,
                   scale=1, usetex=True):
    """Set preferred output style."""
    dim = 0 if dim == 'w' else 1

    a4 = [21.0, 29.7]

    inches_per_cm = 1. / 2.54
    # inches_per_pt = 1.0 / 72.27  # pt/inch (latex)
    # goldenMean = (1.0 + np.sqrt(5.0)) / 2.0

    textwidth = (a4[0] - paperMargin) * inches_per_cm

    fig_width = textwidth * xScale  # fig width in inches
    fig_height = textwidth * yScale  # fig height in inches

    fig_size = [fig_width * scale, fig_height * scale]

    # print "figsize:", fig_size
    # fig.set_size_inches(fig_size)

    # from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # rc('font',**{'family':'serif','serif':['Palatino']})

    params = {
        'backend': 'ps',
        # 'font.weight'       : 'bold',
        'ax.labelsize': fontsize * scale,
        'font.size': fontsize * scale,
        'legend.fontsize': fontsize * scale,
        'xtick.labelsize': fontsize * scale,
        'ytick.labelsize': fontsize * scale,
        # font.sans-serif     : Bitstream Vera Sans, ...
        # 'font.cmb10'     : 'cmb10',
        # 'font.family'         : 'cursive',
        'font.family': 'sans-serif',
        # 'font.sans-serif'   : 'Helvetica',
        'text.usetex': usetex,
        'figure.figsize': fig_size,
        'xtick.major.pad': 4 * scale,
        'xtick.minor.pad': 4 * scale,
        'ytick.major.pad': 4 * scale,
        'ytick.minor.pad': 4 * scale,
        'xtick.major.size': 4 * scale,  # major tick size in points
        'xtick.minor.size': 2 * scale,  # minor tick size in points
        'ytick.major.size': 4 * scale,  # major tick size in points
        'ytick.minor.size': 2 * scale,  # minor tick size in points
        'lines.markersize': 6 * scale,
        'lines.linewidth': 0.6 * scale
    }
    import matplotlib
    matplotlib.rcParams.update(params)


def setPlotStuff(fontsize=7, dpi=None):
    """TODO merge with setOutputStyle.

    Change ugly name.
    """
    from matplotlib import rcParams

    # print(rcParams.keys())

    # rcParams['ax.labelsize'] = fontsize  # REMOVED IN MPL.1.5
    # rcParams['ax.titlesize'] = fontsize  # REMOVED IN MPL.1.5
    # rcParams['ax.linewidth'] = 0.3  # REMOVED IN MPL.1.5
    rcParams['font.size'] = fontsize
    rcParams['xtick.labelsize'] = fontsize
    rcParams['ytick.labelsize'] = fontsize
    rcParams['legend.fontsize'] = fontsize
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica']  # ['Times New Roman']
    rcParams['text.usetex'] = False

    #    rcParams['figure.figsize'] = 7.3, 4.2
    rcParams['xtick.major.size'] = 3
    rcParams['xtick.major.width'] = 0.3
    rcParams['xtick.minor.size'] = 1.5
    rcParams['xtick.minor.width'] = 0.3
    rcParams['ytick.major.size'] = rcParams['xtick.major.size']
    rcParams['ytick.major.width'] = rcParams['xtick.major.width']
    rcParams['ytick.minor.size'] = rcParams['xtick.minor.size']
    rcParams['ytick.minor.width'] = rcParams['xtick.minor.width']

    if dpi is not None:
        rcParams['figure.dpi'] = dpi
        rcParams['savefig.dpi'] = dpi


def createAnimation(fig, animate, nFrames, dpi, out):
    """Create animation for the content of a given matplotlib figure.

    Until I know a better place.
    """
    from matplotlib.animation import FuncAnimation
    anim = FuncAnimation(fig, animate, frames=nFrames,
                         interval=0.001, repeat=False)
    anim.save(out + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
              bitrate=24 * 1024, extra_args=None, metadata=None,
              extra_anim=None, savefig_kwargs=None)
    try:
        print("Create frames ... ")
        os.system('mkdir -p anim-' + out)
        os.system('ffmpeg -i ' + out + '.mp4 anim-' + out + '/movie%d.jpg')
    except BaseException as _:
        pass


def saveAnimation(mesh, data, out, vData=None, plc=None, label='', cMin=None,
                  cMax=None, logScale=False, cmap=None, **kwargs):
    """Create and save an animation for a given mesh with a set of field data.

    Until I know a better place.
    """
    import matplotlib.pyplot as plt

    dpi = 92
    scale = 1
    fig = pg.plt.figure(facecolor='white',
                     figsize=(scale * 800 / dpi, scale * 490 / dpi), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)

    gci = pg.viewer.mpl.drawModel(ax, mesh, data=data[0], cMin=cMin, cMax=cMax,
                                  cMap=cmap, logScale=logScale)

    pg.viewer.mpl.createColorbar(gci, label=label, pad=0.55)

    if plc:
        pg.show(plc, ax=ax)

    adjustWorldAxes(ax)

    def animate(i):
        """TODO WRITEME."""
        print(out + ": Frame:", i, "/", len(data))

        if vData is not None:
            ax.clear()
            pg.hold(True)
            pg.viewer.mpl.drawModel(ax, mesh, data=data[i], cMin=cMin,
                                   cMax=cMax, cMap=cmap, logScale=logScale)
            pg.viewer.mpl.drawStreams(ax, mesh, vData[i], **kwargs)
        else:
            print(min(data[i]), max(data[i]))
            pg.viewer.mpl.setMappableData(gci, data[i], cMin=cMin, cMax=cMax,
                                         logScale=logScale)

        plt.pause(0.001)

    createAnimation(fig, animate, int(len(data)), dpi, out)


def plotLines(ax, line_filename, linewidth=1.0, step=1):
    """Read lines from file and plot over model."""
    xz = np.loadtxt(line_filename)
    n_points = xz.shape[0]
    if step == 2:
        for i in range(0, n_points, step):
            x = xz[i:i + step, 0]
            z = xz[i:i + step, 1]
            ax.plot(x, z, 'k-', linewidth=linewidth)
    if step == 1:
        ax.plot(xz[:, 0], xz[:, 1], 'k-', linewidth=linewidth)


def twin(ax):
    """Return the twin of ax if exist."""
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            return other_ax
    return None


def createTwinX(ax):
    """Create (or return existing) twin x axes for ax."""
    return _createTwin(ax, 'twinx')


def createTwinY(ax):
    """Create (or return existing) twin x axes for ax."""
    return _createTwin(ax, 'twiny')


def _createTwin(ax, funct):
    """Create (or return existing) twin x axes for ax."""
    tax = None
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            tax = other_ax

    if tax is None:
        tax = getattr(ax, funct)()

    return tax

def isInteractive():
    """Return False if a non-interactive backend is used.

    e.g. for Jupyter Notebooks and sphinx builds.
    """
    backend = pg.plt.get_backend().lower()
    inlineBackend = "inline" in backend or backend == "agg"
    return not inlineBackend
