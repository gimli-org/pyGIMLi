"""Kernel computation for gravity and magnetics."""
import numpy as np
import pygimli as pg
from pygimli.utils import ProgressBar


@pg.cache
def SolveGravMagHolstein(mesh, pnts, cmp, igrf=None):  # , foot=np.inf):
    """Solve gravity and/or magnetics problem after Holstein (1997).

    Parameters
    ----------
    mesh : pygimli:mesh
        tetrahedral or hexahedral mesh
    pnts : list|array of (x, y, z)
        measuring points
    cmp : list of str
        component list of type str, valid values are:
        gx, gy, gz, TFA, Bx, By, Bz, Bxx, Bxy, Bxz, Byy, Byz, Bzz
    igrf : list|array of size 3 or 7
        international geomagnetic reference field, either
        [D, I, H, X, Y, Z, F] - declination, inclination, horizontal field,
                               X/Y/Z components, total field OR
        [X, Y, Z] - X/Y/Z components

    Returns
    -------
    out : ndarray (nPoints x nComponents x nCells)
        kernel matrix to be multiplied with density or susceptibility
    """
    if pnts is None:
        pnts = [[0.0, 0.0]]

    doG = np.any([c[0] == "g" for c in cmp])
    doGT = np.any([c[0] == "g" and len(c) == 3 for c in cmp])
    doB = np.any([c[0] == "B" and len(c) == 2 for c in cmp]) or "TFA" in cmp
    doBT = np.any([c[0] == "B" and len(c) == 3 for c in cmp])
    B_tens = None

    kernel = np.zeros((mesh.cellCount(), len(pnts), len(cmp)))
    if np.any(igrf):
        if len(igrf) == 3:  # an X, Y, Z vector
            F = np.linalg.norm(igrf)
            fakt = F / (4*np.pi)
            B_dir = np.array(igrf) / F
        elif len(igrf) == 7:  # an IGRF vector (D, I, H, X, Y, Z, F)
            fakt = igrf[6] / (4*np.pi)
            myigrf = np.array(igrf[3:6])
            B_dir = myigrf / np.linalg.norm(myigrf)
        else:
            raise ValueError("Could not use IGRF vector. Len must be 3 or 7!")
    elif doB or doB:
        raise KeyError("Specify IGRF!")

    b_list, c_list = [], []
    for bd in mesh.boundaries():
        b_list.append([n.id() for n in bd.allNodes()])
        c_list.append([bd.leftCell(), bd.rightCell()])

    b_list = np.array(b_list)
    lb = b_list.shape

    n_list = np.array([n.pos() for n in mesh.nodes()])

    cl, cr = [], []
    for i, c in enumerate(c_list):
        if c[0]:
            cl.append([i, c[0].id()])
        if c[1]:
            cr.append([i, c[1].id()])

    cl = np.array(cl)
    cr = np.array(cr)

    rr = range(0, mesh.cellCount())
    rs = np.roll(range(0, lb[1]), -1)

    temp = np.zeros((len(pnts), lb[0], len(cmp)))
    pBar = ProgressBar(its=len(pnts), width=40, sign='+')
    nb = n_list[b_list]
    for i, p in enumerate(pnts):
        r1 = nb - p
        r2 = r1[:, rs, :]
        r0 = r2 - r1
        u = np.sum(np.cross(r1, r2), 1)
        u /= np.expand_dims(np.linalg.norm(u, axis=1)+1e-16, axis=1)
        ut = np.tile(u, lb[1]).reshape((lb[0], lb[1], 3))  # broadcasting!
        ll = np.linalg.norm(r0, axis=2)
        t = r0 / np.expand_dims(ll, axis=2)
        lm = (np.sum(r1*t, 2) + np.sum(r2*t, 2)) / 2
        h = np.cross(t, ut)
        hn = np.sum(h*r1, 2)
        v = np.sum(ut*r1, 2)
        r1n = np.linalg.norm(r1, axis=2)
        r2n = np.linalg.norm(r2, axis=2)
        rm = (r1n + r2n)/2
        lumbda = ll / (2*rm)
        jj = 0
        if doG: # gravitational field
            g = hn*np.arctanh(lumbda)-np.sign(v)*v*np.arctan2(
                hn*lumbda, (rm*(1-lumbda**2)+abs(v)))
            g_vec = 2 * u * np.expand_dims(np.sum(g, 1), axis=1)

            if 'g' in cmp:
                temp[i, :, jj] = g
                jj += 0

            if 'gx' in cmp:
                temp[i, :, jj] = g_vec[:, 0]
                jj += 1

            if 'gy' in cmp:
                temp[i, :, jj] = g_vec[:, 1]
                jj += 1

            if 'gz' in cmp:
                temp[i, :, jj] = g_vec[:, 2]
                jj += 1

            if doGT:
                raise NotImplementedError("Gravity tensor not yet supported!")
                G_tens = np.zeros([3, 3, 3])
                if 'gxx' in cmp:
                    temp[:, jj]=G_tens[:, 0, 0]
                if 'gxy' in cmp:
                    temp[:, jj]=G_tens[:, 0, 1]
                if 'gxz' in cmp:
                    temp[:, jj]=G_tens[:, 0, 2]
                if 'gyy' in cmp:
                    temp[:, jj]=G_tens[:, 1, 1]
                if 'gyz' in cmp:
                    temp[:, jj]=G_tens[:, 1, 2]
                if 'gzz' in cmp:
                    temp[:, jj]=G_tens[:, 2, 2]

        if doB or doBT:
            # magnetic field vector and gravity gradient tensor
            b = h*np.expand_dims(np.arctanh(lumbda), axis=2) - \
                ut*np.expand_dims(np.sign(v)*np.arctan2(
                    hn*lumbda, (rm*(1-lumbda**2)+abs(v))), axis=2)
            P = np.dot(u, B_dir)
            B_vec = np.expand_dims(P, 1) * np.sum(b, 1)
            B_vec = 2 * np.expand_dims(P, 1) * np.sum(b, 1)

            if 'TFA' in cmp:
                temp[i, :, jj] = fakt*B_vec.dot(B_dir)
                jj += 1

            if 'Bx' in cmp:
                temp[i, :, jj] = fakt*B_vec[:, 0]
                jj += 1

            if 'By' in cmp:
                temp[i, :, jj] = fakt*B_vec[:, 1]
                jj += 1

            if 'Bz' in cmp:
                temp[i, :, jj] = fakt*B_vec[:, 2]
                jj += 1

            if doBT:  # magnetic gradient tensor
                d = (-2*lumbda*hn) / (r1n*r2n*(1-lumbda**2))
                e = (-2*lumbda*lm) / (r1n*r2n)
                f = (-2*lumbda*v) / (r1n*r2n*(1-lumbda**2))

                h1 = np.expand_dims(h, axis=3)
                h2 = np.swapaxes(h1, 2, 3)
                t1 = np.expand_dims(t, axis=3)
                t2 = np.swapaxes(t1, 2, 3)
                u1 = np.expand_dims(ut, axis=3)
                u2 = np.swapaxes(u1, 2, 3)

                B = (h1*h2-u1*u2)*np.expand_dims(d, (2, 3)) + \
                    (t1*h2+h1*t2)*np.expand_dims(e, (2, 3))/2 + \
                    (h1*u2+u1*h2)*np.expand_dims(f, (2, 3))

                B_tens = np.expand_dims(P, (1, 2)) * np.sum(B, 1)

                if 'Bxx' in cmp:
                    temp[i, :, jj] = fakt*B_tens[:, 0, 0]
                    jj += 1

                if 'Bxy' in cmp:
                    temp[i, :, jj] = fakt*B_tens[:, 0, 1]
                    jj += 1

                if 'Bxz' in cmp:
                    temp[i, :, jj] = fakt*B_tens[:, 0, 2]
                    jj += 1

                if 'Byy' in cmp:
                    temp[i, :, jj] = fakt*B_tens[:, 1, 1]
                    jj += 1

                if 'Byz' in cmp:
                    temp[i, :, jj] = fakt*B_tens[:, 1, 2]
                    jj += 1

                if 'Bzz' in cmp:
                    temp[i, :, jj] = fakt * B_tens[:, 2, 2]
                    jj += 1

        pBar.update(i)

    kernel += np.array([np.sum(temp[:, cl[cl[:, 1] == j, 0]], 1) for j in rr])
    kernel -= np.array([np.sum(temp[:, cr[cr[:, 1] == j, 0]], 1) for j in rr])

    return kernel.transpose([1, 2, 0])
