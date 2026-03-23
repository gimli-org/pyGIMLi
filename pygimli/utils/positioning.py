"""Several tools for position-based."""
import numpy as np


def distToLine(x, y, lx, ly):
    """Compute minimum distance to segmented Line.
    
    Parameters
    ----------
    x, y : iterable
        x and y positions
    lx, ly : iterable
        position of points creating a line
    """
    assert len(x) == len(y), "Lengths of x/y not matching"
    assert len(lx) == len(ly), "Lengths of lx/ly not matching"
    dist2 = np.ones_like(x) * 1e9
    for i in range(len(lx)-1):
        px = x - lx[i]
        py = y - ly[i]
        bx = lx[i + 1] - lx[i]
        by = ly[i + 1] - ly[i]
        t = (px * bx + by * py) / (bx**2 + by**2)
        t = np.maximum(np.minimum(t, 1), 0)
        dist2 = np.minimum(dist2, (px-t*bx)**2+(py-t*by)**2)

    return np.sqrt(dist2)


def pointInsidePolygon(x, y, polygon):
    """Determine whether a point is inside a closed polygon."""
    n = len(polygon)
    inside = False
    p1x, p1y, *_ = polygon[0]
    for i in range(n + 1):
        p2x, p2y, *_ = polygon[i % n]
        if y > min(p1y, p2y) and y <= max(p1y, p2y) and x <= max(p1x, p2x):
            if p1y != p2y:
                xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x or x <= xinters:
                    inside = not inside

        p1x, p1y = p2x, p2y

    return inside


def detectLinesAlongAxis(rx, ry, axis='x'):
    """Split data in lines based on x/y axis."""

    if axis == 'x':
        r = rx
    elif axis == 'y':
        r = ry
    else:
        print('Choose either *x* or *y* axis. Aborting this method ...')
        return

    dummy = np.zeros_like(rx, dtype=int)
    line = np.zeros_like(rx, dtype=int)
    li = 0
    last_sign = np.sign(r[1] - r[0])
    for ri in range(1, len(rx)):
        sign = np.sign(r[ri] - r[ri-1])
        dummy[ri-1] = li
        if sign != last_sign:
            li += 1
            last_sign *= -1
    dummy[-1] = li

    return sortLines(rx, ry, line, dummy, axis)


def detectLinesByDistance(rx, ry, minDist=200., axis='x'):
    """Split data in lines based on point distances."""

    dummy = np.zeros_like(rx, dtype=int)
    line = np.zeros_like(rx, dtype=int)
    li = 0
    for ri in range(1, len(rx)):
        dummy[ri-1] = li
        dist = np.sqrt((rx[ri]-rx[ri-1])**2 +
                       (ry[ri]-ry[ri-1])**2)
        if dist > minDist:
            li += 1
    
    dummy[-1] = li

    return sortLines(rx, ry, line, dummy, axis)


def detectLinesBySpacing(rx, ry, vec, axis='x'):
    """Detect line changes by jumps in point spacing."""
    if axis == 'x':
        r = rx
    elif axis == 'y':
        r = ry
    else:
        print('Choose either *x* or *y* axis. Aborting this method ...')
        return

    return np.argmin(np.abs(np.tile(r, (len(vec), 1)).T - vec), axis=1)


def detectLinesByDirection(rx, ry):
    """Split data in lines for line-wise processing."""
    dt = np.sqrt(np.diff(rx)**2 + np.diff(ry)**2)
    dtmin = np.median(dt) * 2
    dx = np.round(np.diff(rx) / dt * 2)
    dy = np.round(np.diff(ry) / dt * 2)
    sdx = np.hstack((0, np.diff(np.sign(dx)), 0))
    sdy = np.hstack((0, np.diff(np.sign(dy)), 0))
    line = np.zeros_like(rx, dtype=int)
    nLine = 1
    act = True
    for i, sdxi in enumerate(sdx):
        if sdxi != 0:
            act = not act
            if act:
                nLine += 1

        if sdy[i] != 0:
            act = not act
            if act:
                nLine += 1

        if i > 0 and dt[i-1] > dtmin:
            act = True
            nLine += 1

        if act:
            line[i] = nLine

    return line


def sortLines(rx, ry, line, dummy, axis):
    """Sort line elements by Rx or Ry coordinates."""
    means = []
    for li in np.unique(dummy):
        if axis == 'x':
            means.append(np.mean(ry[dummy==li], axis=0))
        elif axis == 'y':
            means.append(np.mean(rx[dummy==li], axis=0))

    lsorted = np.argsort(means)
    for li, lold in enumerate(lsorted):
        line[dummy == lold] = li + 1

    return line

def detectLines(x, y, mode=None, axis='x', show=False):
    """Split data in lines for line-wise processing.

    Several modes are available:
        'x'/'y': along coordinate axis
        spacing vector: by given spacing
        float: minimum distance
    """
    if isinstance(mode, (str)):
        line = detectLinesAlongAxis(x, y, axis=mode)
    elif hasattr(mode, "__iter__"):
        line = detectLinesBySpacing(x, y, mode, axis=axis)
    elif isinstance(mode, (int, float)):
        line = detectLinesByDistance(x, y, mode,
                                            axis=axis)
    else:
        line = detectLinesByDirection(x, y)

    return line
