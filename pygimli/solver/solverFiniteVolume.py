#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple Finite Volume Solver."""

import numpy as np
import pygimli as pg

from pygimli.solver import identity  # , parseArgToArray, parseArgToBoundaries


def cellDataToBoundaryData(mesh, v):
    """TODO Documentme."""
    if len(v) != mesh.cellCount():
        raise Exception("len(v) != mesh.cellCount():", len(v),
                        mesh.cellCount())
    return np.array([cellToFace(b, v) for b in mesh.boundaries()])

def boundaryNormals(mesh):
    """Collect all boundary outer normal vectors."""
    # implement with    return mesh.boundaryNorms()
    gB = np.zeros((mesh.boundaryCount(), 3))

    for b in mesh.boundaries():
        gB[b.id()] = b.norm()

    return gB

def _boundaryToCellDistances(mesh):
    """TODO Documentme."""
    d = [_boundaryToCellDistancesBound(b) for b in mesh.boundaries()]
    return np.array(d)

def _boundaryToCellDistancesBound(b):
    """TODO Documentme."""
    leftCell = b.leftCell()
    rightCell = b.rightCell()
    df1 = 0.
    df2 = 0.
    if leftCell:
        df1 = b.center().distance(leftCell.center())
    if rightCell:
        df2 = b.center().distance(rightCell.center())
    return df1 + df2


def cellDataToBoundaryGrad(mesh, v, vGrad):
    """TODO Documentme."""
    if len(v) != mesh.cellCount() or len(vGrad) != mesh.cellCount():
        raise BaseException("len(v) dismatch mesh.cellCount()")

    gB = mesh.cellDataToBoundaryGradient(v, vGrad)
    return np.vstack([pg.x(gB), pg.y(gB), pg.z(gB)]).T

    # gB = np.zeros((mesh.boundaryCount(), 3))
    #
    # for b in mesh.boundaries():
    #     leftCell = b.leftCell()
    #     rightCell = b.rightCell()
    #     gr = pg.RVector3(0.0, 0.0, 0.0)
    #     t = (b.node(1).pos() - b.node(0).pos()).norm()
    #
    #     if leftCell and rightCell:
    #         df1 = b.center().distance(leftCell.center())
    #         df2 = b.center().distance(rightCell.center())
    #
    #         gr = b.norm() * \
    #             (v[rightCell.id()] - v[leftCell.id()]) / (df1 + df2)
    #
    #         grL = t * t.dot(vGrad[leftCell.id()])
    #         grR = t * t.dot(vGrad[rightCell.id()])
    #
    #         gr += (grL + grR) * 0.5
    #
    #     elif leftCell:
    #         gr = t * t.dot(vGrad[leftCell.id()])
    #
    #     gB[b.id(), 0] = gr[0]
    #     gB[b.id(), 1] = gr[1]
    #     gB[b.id(), 2] = gr[2]
    # return gB


def cellToFace(boundary, vec, harmonic=False):
    """DEPRECATED.

    Interpolate cell to face values by weighted arithmetic/harmonic mean.
    """
    DEPRECATED
    leftCell = boundary.leftCell()
    rightCell = boundary.rightCell()
    df1 = 0.
    df2 = 0.
    u1 = 0.
    u2 = 0.
    if leftCell:
        df1 = boundary.center().distance(leftCell.center())
        u1 = vec[leftCell.id()]
    if rightCell:
        df2 = boundary.center().distance(rightCell.center())
        u2 = vec[rightCell.id()]
    uFace = 0.0
    d12 = (df1 + df2)

    if leftCell and rightCell:
        if harmonic:
            # harmonic mean
            uFace = (u1 * u2) / ((u2 - u1) * df2 / d12 + u1)
        else:
            # arithmetic mean
            # check left vs. right
            uFace = (u1 - u2) * df2 / d12 + u2

    elif leftCell:
        uFace = u1  # / df1
    elif rightCell:
        uFace = u2  # / df2
    return uFace


def cellToFaceArithmetic(boundary, AMM):
    """TODO Documentme."""
    leftCell = boundary.leftCell()
    rightCell = boundary.rightCell()
    df1 = 0.
    df2 = 0.
    harmonic = False
    if leftCell:
        df1 = boundary.center().distance(leftCell.center())
    if rightCell:
        df2 = boundary.center().distance(rightCell.center())
    d12 = (df1 + df2)

    if leftCell and rightCell:
        if harmonic:
            pass
            # harmonic mean
            # uFace = (u1 * u2) / ((u2-u1)*df2/d12 + u1)
        else:
            # arithmetic mean
            # check left vs. right
            AMM.addVal(boundary.id(), leftCell.id(), df2 / d12)
            AMM.addVal(boundary.id(), rightCell.id(), -df2 / d12 + 1.0)
            # uFace = (u1 - u2) * df2/d12 + u2

    elif leftCell:
        AMM.addVal(boundary.id(), leftCell.id(), 1.0)
    elif rightCell:
        AMM.addVal(boundary.id(), rightCell.id(), 1.0)


def cellDataToBoundaryDataMatrix(mesh):
    """TODO Documentme."""
    AMM = pg.matrix.SparseMapMatrix(mesh.boundaryCount(), mesh.cellCount())

    for b in mesh.boundaries():
        cellToFaceArithmetic(b, AMM)

    return AMM


def cellDataToCellGrad2(mesh, v):
    """TODO Documentme."""
    if len(v) != mesh.cellCount():
        print(len(v), mesh.cellCount())
        raise BaseException("len(v) dismatch mesh.cellCount()")

    vN = pg.meshtools.cellDataToNodeData(mesh, v)
    gC = np.zeros((mesh.cellCount(), 3))

    for c in mesh.cells():
        gr = c.grad(c.center(), vN)
        gC[c.id(), 0] = gr[0]
        gC[c.id(), 1] = gr[1]
        gC[c.id(), 2] = gr[2]
    return gC


def cellDataToCellGrad(mesh, v, CtB):
    """TODO Documentme."""
    if len(v) != mesh.cellCount():
        print(len(v), mesh.cellCount())
        raise BaseException("len of v missmatch mesh.cellCount()")
    div = mesh.boundaryDataToCellGradient(CtB * v)
    return np.vstack([pg.x(div), pg.y(div), pg.z(div)]).T

    # vF = cellDataToBoundaryData(mesh, v)
    # gC = np.zeros((mesh.cellCount(), 3))
    #
    # for b in mesh.boundaries():
    #
    #     leftCell = b.leftCell()
    #     rightCell = b.rightCell()
    #     vec = b.norm() * vF[b.id()] * b.size()
    #
    #     if leftCell:
    #         gC[leftCell.id(), 0] += vec[0]
    #         gC[leftCell.id(), 1] += vec[1]
    #         gC[leftCell.id(), 2] += vec[2]
    #     if rightCell:
    #         gC[rightCell.id(), 0] -= vec[0]
    #         gC[rightCell.id(), 1] -= vec[1]
    #         gC[rightCell.id(), 2] -= vec[2]
    #
    # gC[:, 0] /= mesh.cellSizes()
    # gC[:, 1] /= mesh.cellSizes()
    # gC[:, 2] /= mesh.cellSizes()
    #
    # return gC


def findVelocity(mesh, v, b, c, nc=None):
    """Find velocity for boundary b based on vector field v.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`

    v : array [(N,dim)]
        velocity field [[v_i,]_j,] with i=[1..3] for the mesh dimension
        and j = [0 .. N-1] per Cell or per Node so N is either
        mesh.cellCount() or mesh.nodeCount()

    b : :gimliapi:`GIMLI::Boundary`
        Boundary

    c : :gimliapi:`GIMLI::Cell`

        Associated Cell in flow direction

    nc : :gimliapi:`GIMLI::Cell`
        Associated neighbor cell .. if one given
        from flow direction
    """
    vel = [0.0, 0.0, 0.0]

    if hasattr(v, '__len__'):
        if len(v) == mesh.cellCount():
            if nc:
                # mean cell based vector-field v[x,y,z] for cell c and cell nc
                vel = (v[c.id()] + v[nc.id()]) / 2.0

                # vel1 = 1./ (1./v[c.id()] + 1./v[nc.id()])
                # print("findVel, check:", vel1, vel)

            else:
                # cell based vector-field v[x,y,z] for cell c
                vel = v[c.id()]
        elif len(v) == mesh.boundaryCount():
            vel = v[b.id()]
        else:
            # interpolate node based vector-field v[x,y,z] at point b.center()
            # print(len(vel), vel)
            vel = c.vec(b.center(), v)

    return vel


def findDiffusion(mesh, a, b, c, nc=None):
    """TODO Documentme.

    Parameters
    ----------
    a :
        Attribute for diffusion coefficient. Cell based or Boundary based

    b : :gimliapi:`GIMLI::Boundary`
        Boundary

    c : :gimliapi:`GIMLI::Cell`

        Associated cell in flow direction

    nc : :gimliapi:`GIMLI::Cell`
        associated neighbor cell .. if one given
        from flow direction
    """
    D = 0.

    if nc:
        if len(a) == mesh.boundaryCount():
            D = a[b.id()] / nc.center().distance(c.center()) * b.size()
        else:
            # Diffusion part
            # Interface harmonic median
            # D = (a[c.id()] / c.center().distance(b.center()) +
            # a[nc.id()] / nc.center().distance(b.center())) * 0.5 *
            # b.size()
            if a[c.id()] > 0 and a[nc.id()] > 0:
                D = 1. / (c.center().distance(b.center()) / a[c.id()] +
                        nc.center().distance(b.center()) / a[nc.id()]) * b.size()
            # print(D)
    else:
        if len(a) == mesh.boundaryCount():
            D = a[b.id()] / b.center().distance(c.center()) * b.size()
        else:
            D = a[c.id()] / b.center().distance(c.center()) * b.size()
    return D


def diffusionConvectionKernel(mesh, a=None, b=0.0,
                              uB=None, duB=None,
                              vel=0,
                              # u0=0,
                              fn=None,
                              scheme='CDS', sparse=False, time=0.0,
                              userData=None):
    """
    Generate system matrix for diffusion and convection in a velocity field.

    Particle concentration u inside a velocity field.

    Peclet Number - ratio between convection/diffusion = F/D
        F = velocity flow trough volume boundary,
        D = diffusion coefficient

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    a: value | array | callable(cell, userData)
        Diffusion coefficient per cell
    b: value | array | callable(cell, userData)
        TODO What is b
    fn: iterable(cell)
        TODO What is fn
    vel: ndarray (N,dim) | RMatrix(N,dim)
        velocity field [[v_i,]_j,] with i=[1..3] for the mesh dimension
        and j = [0 .. N-1] per Cell or per Node so N is either
        mesh.cellCount() or mesh.nodeCount()
    scheme: str [CDS]
        Finite volume scheme
        * CDS -- Central Difference Scheme.
            maybe irregular for Peclet no. |F/D| > 2
            Diffusion dominant. Error of order 2
        * UDS -- Upwind Scheme.
            Convection dominant. Error of order 1
        * HS -- Hybrid Scheme.
            Diffusion dominant for Peclet-number |(F/D)| < 2
            Convection dominant else.
        * PS -- Power Law Scheme.
            Identical to HS for Peclet-number |(F/D)| > 10 and near to ES else
            Convection dominant.
        * ES -- Exponential scheme
            Only stationary one-dimensional but exact solution
    Returns
    -------
    S: :gimliapi:`GIMLI::SparseMatrix` | numpy.ndarray(nCells, nCells)
        Kernel matrix, depends on vel, a, b, scheme, uB, duB .. if some of this
        has been changed you cannot cache these matrix
    rhsBoundaryScales: ndarray(nCells)
        RHS offset vector
    """
    if a is None:
        a = pg.Vector(mesh.boundaryCount(), 1.0)

    AScheme = None
    if scheme == 'CDS':
        AScheme = lambda peclet_: 1.0 - 0.5 * abs(peclet_)
    elif scheme == 'UDS':
        AScheme = lambda peclet_: 1.0
    elif scheme == 'HS':
        AScheme = lambda peclet_: max(0.0, 1.0 - 0.5 * abs(peclet_))
    elif scheme == 'PS':
        AScheme = lambda peclet_: max(0.0, (1.0 - 0.1 * abs(peclet_))**5.0)
    elif scheme == 'ES':
        AScheme = lambda peclet_: (peclet_) / (np.exp(abs(peclet_)) - 1.0) \
            if peclet_ != 0.0 else 1
    else:
        raise BaseException("Scheme unknwon:" + scheme)

    useHalfBoundaries = False

    dof = mesh.cellCount()

    if not uB:
        uB = []
    if not duB:
        duB = []

    if useHalfBoundaries:
        dof = mesh.cellCount() + len(uB)

    S = None
    if sparse:
        S = pg.matrix.SparseMapMatrix(dof, dof, stype=0) + identity(dof, scale=b)
    else:
        S = np.zeros((dof, dof))

    rhsBoundaryScales = np.zeros(dof)

    # we need this to fast identify uBoundary and value by boundary
    uBoundaryID = []
    uBoundaryVals = [None] * mesh.boundaryCount()

    for [b, val] in uB:

        if isinstance(b, pg.core.Boundary):
            uBoundaryID.append(b.id())
            uBoundaryVals[b.id()] = val
        elif isinstance(b, pg.core.Node):
            for _b in b.boundSet():
                if _b.rightCell() is None:
                    pg.warn('Dirichlet for one node considered for the nearest boundary.', _b.id())
                    uBoundaryID.append(_b.id())
                    uBoundaryVals[_b.id()] = val
                    break
        else:
            raise BaseException("Please give boundary, value list")


    duBoundaryID = []
    duBoundaryVals = [None] * mesh.boundaryCount()

    for [boundary, val] in duB:
        if not isinstance(boundary, pg.core.Boundary):
            raise BaseException("Please give boundary, value list")

        duBoundaryID.append(boundary.id())
        duBoundaryVals[boundary.id()] = val

    # iterate over all cells
    for cell in mesh.cells():
        cID = cell.id()
        for bi in range(cell.boundaryCount()):
            boundary = pg.core.findBoundary(cell.boundaryNodes(bi))

            ncell = boundary.leftCell()
            if ncell == cell:
                ncell = boundary.rightCell()

            v = findVelocity(mesh, vel, boundary, cell, ncell)

            # Convection part
            F = boundary.norm(cell).dot(v) * boundary.size()
            # print(F, boundary.size(), v, vel)
            # Diffusion part
            D = findDiffusion(mesh, a, boundary, cell, ncell)

            # print(F, D, F/D)
            # print((1.0 - 0.1 * abs(F/D))**5.0)
            if D > 0:
                aB = D * AScheme(F / D) + max(-F, 0.0)
            else:
                aB = max(-F, 0.0)

            aB /= cell.size()

            # print(cell.center(), boundary.center(), boundary.norm(cell), aB)
            if ncell:
                # no boundary
                if sparse:
                    S.addVal(cID, ncell.id(), -aB)
                    S.addVal(cID, cID, +aB)
                else:
                    S[cID, ncell.id()] -= aB
                    S[cID, cID] += aB

            elif not useHalfBoundaries:

                if boundary.id() in uBoundaryID:
                    val = pg.solver.generateBoundaryValue(boundary,
                                                uBoundaryVals[boundary.id()],
                                                time=time,
                                                userData=userData)

                    if sparse:
                        S.addVal(cID, cID, aB)
                    else:
                        S[cID, cID] += aB

                    #print(aB, uBoundaryVals[boundary.id()])
                    rhsBoundaryScales[cID] += aB * np.mean(val)

                if boundary.id() in duBoundaryID:
                    # Neumann boundary condition
                    val = pg.solver.generateBoundaryValue(
                        boundary,
                        duBoundaryVals[boundary.id()],
                        time=time,
                        userData=userData)

                    val = np.mean(val)

                    # amount of flow through the boundary .. maybe buggy
                    # fill be replaced by suitable FE solver
                    outflow = -val * boundary.size() / cell.size()

                    if sparse:
                        S.addVal(cID, cID, outflow)
                    else:
                        S[cID, cID] += outflow

        if fn is not None:
            if sparse:
                # * cell.shape().domainSize())
                S.addVal(cell.id(), cell.id(), -fn[cell.id()])
            else:
                # * cell.shape().domainSize()
                S[cell.id(), cell.id()] -= fn[cell.id()]

    return S, rhsBoundaryScales


def solveFiniteVolume(mesh, a=1.0, b=0.0, f=0.0, fn=0.0, vel=None, u0=0.0,
                      times=None,
                      uL=None, relax=1.0,
                      ws=None, scheme='CDS', **kwargs):
    r"""Solve partial differential equation with Finite Volumes.

    This function is a syntactic sugar proxy for using the Finite Volume
    functionality of the library core to solve elliptic and parabolic partial
    differential of the following type:

    .. math::

        \frac{\partial u}{\partial t} + \mathbf{v}\cdot\nabla u & = \nabla\cdot(a \nabla u) + b u + f(\mathbf{r},t)\\
        u(\mathbf{r}, t) & = u_B  \quad\mathbf{r}\in\Gamma_{\text{Dirichlet}}\\
        \frac{\partial u(\mathbf{r}, t)}{\partial \mathbf{n}} & = u_{\partial \text{B}}  \quad\mathbf{r}\in\Gamma_{\text{Neumann}}\\
        u(\mathbf{r}, t=0) & = u_0 \quad\text{with} \quad\mathbf{r}\in\Omega

    The Domain :math:`\Omega` and the Boundary :math:`\Gamma` are defined
    through the given mesh with appropriate boundary marker.

    The solution :math:`u(\mathbf{r}, t)` is given for each cell in the mesh.

    TODO:

     * Refactor with solver class and Runga-Kutte solver
     * non steady boundary conditions

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    a: value | array | callable(cell, userData)
        Stiffness weighting per cell values.
    b: value | array | callable(cell, userData)
        Scale for mass values b
    f: iterable(cell)
        Load vector
    fn: iterable(cell)
        TODO What is fn
    vel: ndarray (N,dim) | RMatrix(N,dim)
        Velocity field :math:`\mathbf{v}(\mathbf{r}, t=\text{const}) = \{[v_i]_j,\}`
        with :math:`i=[1\ldots 3]` for the mesh dimension
        and :math:`j = [0\ldots N-1]` with N either the amount of cells,
        nodes, or boundaries.
        Velocities per boundary are preferred and will be interpolated
        on demand.
    u0: value | array | callable(cell, userData)
        Starting field
    times: iterable
        Time steps to calculate for.
    ws Workspace
        This can be an empty class that will used as an Workspace to store and
        cache data.

        If ws is given: The system matrix is taken from ws or
        calculated once and stored in ws for further usage.

        The system matrix is cached in this Workspace as ws.S
        The LinearSolver with the factorized matrix is cached in
        this Workspace as ws.solver
        The rhs vector is only stored in this Workspace as ws.rhs
    scheme: str [CDS]
        Finite volume scheme:
        :py:mod:`pygimli.solver.diffusionConvectionKernel`
    **kwargs:
        * bc : Boundary Conditions dictionary, see pg.solver
        * uB : Dirichlet boundary conditions DEPRECATED
        * duB : Neumann boundary conditions DEPRECATED

    Returns
    -------
    u: ndarray(nTimes, nCells)
        Solution field for all time steps.
    """
    verbose = kwargs.pop('verbose', False)
    # The Workspace is to hold temporary data or preserve matrix rebuild
    # swatch = pg.Stopwatch(True)
    sparse = True

    workspace = pg.solver.WorkSpace()
    if ws:
        workspace = ws

    a = pg.solver.parseArgToArray(a, [mesh.cellCount(), mesh.boundaryCount()])
    b = pg.solver.parseArgToArray(b, mesh.cellCount())
    f = pg.solver.parseArgToArray(f, mesh.cellCount())
    fn = pg.solver.parseArgToArray(fn, mesh.cellCount())

    boundsDirichlet = None
    boundsNeumann = None

    # BEGIN check for Courant-Friedrichs-Lewy
    if vel is not None:

        if isinstance(vel, float):
            print("Warning! .. velocity is float and no vector field")

        # we need velocities for boundaries
        if len(vel) is not mesh.boundaryCount():
            if len(vel) == mesh.cellCount():
                vel = pg.meshtools.cellDataToNodeData(mesh, vel)
            elif len(vel) == mesh.nodeCount():
                vel = pg.meshtools.nodeDataToBoundaryData(mesh, vel)
            else:
                print("mesh:", mesh)
                print("vel:", vel.shape)

                raise Exception("Cannot determine data format for velocities")

        if times is not None:
            pg.solver.checkCFL(times, mesh, np.max(pg.abs(vel)))

    if not hasattr(workspace, 'S'):

        boundsDirichlet = []
        if 'bc' in kwargs:
            bct = dict(kwargs['bc'])
            if 'Dirichlet' in bct:
                boundsDirichlet += pg.solver.parseArgToBoundaries(bct.pop('Dirichlet'), mesh)

            if 'Node' in bct:
                n = bct.pop('Node')
                boundsDirichlet.append([mesh.node(n[0]), n[1]])

            if 'Neumann' in bct:
                boundsNeumann = pg.solver.parseArgToBoundaries(bct.pop('Neumann'), mesh)


        if 'uB' in kwargs:
            pg.deprecated('use new bc dictionary')
            boundsDirichlet = pg.solver.parseArgToBoundaries(kwargs['uB'],
                                                             mesh)

        if 'duB' in kwargs:
            pg.deprecated('use new bc dictionary')
            boundsNeumann = pg.solver.parseArgToBoundaries(kwargs['duB'], mesh)

        workspace.S, workspace.rhsBCScales = diffusionConvectionKernel(
            mesh=mesh,
            a=a,
            b=b,
            uB=boundsDirichlet,
            duB=boundsNeumann,
            # u0=u0,
            fn=fn,
            vel=vel,
            scheme=scheme,
            sparse=sparse,
            userData=kwargs.pop('userData', None))

        dof = len(workspace.rhsBCScales)
        workspace.ap = np.zeros(dof)

        # for nonlinears
        if uL is not None:
            for i in range(dof):
                val = 0.0
                if sparse:
                    val = workspace.S.getVal(i, i) / relax
                    workspace.S.setVal(i, i, val)
                else:
                    val = workspace.S[i, i] / relax
                    workspace.S[i, i] = val

                workspace.ap[i] = val

        # print('FVM kernel 2:', swatch.duration(True))
    # endif: not hasattr(workspace, 'S'):

    workspace.rhs = np.zeros(len(workspace.rhsBCScales))
    workspace.rhs[0:mesh.cellCount()] = f  # * mesh.cellSizes()

    workspace.rhs += workspace.rhsBCScales

    # for nonlinear: relax progress with scaled last result
    if uL is not None:
        workspace.rhs += (1. - relax) * workspace.ap * uL
    # print('FVM: Prep:', swatch.duration(True))

    if not hasattr(times, '__len__'):

        if sparse and not hasattr(workspace, 'solver'):
            Sm = pg.matrix.SparseMatrix(workspace.S)
            # hold Sm until we have reference counting,
            # loosing Sm here will kill LinSolver later
            workspace.Sm = Sm
            workspace.solver = pg.core.LinSolver(Sm, verbose=verbose)

        u = None
        if sparse:
            u = workspace.solver.solve(workspace.rhs)
        else:
            u = np.linalg.solve(workspace.S, workspace.rhs)
        # print('FVM solve:', swatch.duration(True))
        return u[0:mesh.cellCount():1]
    else:
        theta = kwargs.pop('theta', 0.5 + 1e-6)

        if sparse:
            I = pg.solver.identity(len(workspace.rhs))
        else:
            I = np.diag(np.ones(len(workspace.rhs)))

        progress = None
        if verbose:
            from pygimli.utils import ProgressBar
            progress = ProgressBar(its=len(times))

            print("Solve timesteps with Crank-Nicolson.")

        return pg.solver.crankNicolson(times, workspace.S, I,
                    f=workspace.rhs,
                    u0=pg.solver.cellValues(mesh, u0),
                    theta=theta,
                    progress=progress)


def createFVPostProzessMesh(mesh, u, uDirichlet):
    """Create node based post processing of cell centered FV solutions.

    Create a mesh suitable for node based post processing of cell
    centered Finite Volume solutions.
    This is something like cellDataToNodeData with extra Dirichlet points
    but without smoothing due to interpolation.

    IMPROVE DOC!!

    Parameters
    ----------

    """
    allBounds = pg.solver.parseArgToBoundaries(uDirichlet, mesh)
    bounds, vals = zip(*allBounds)
    uDirVals = pg.solver.generateBoundaryValue(bounds, vals)

    def isBoundary(b):
        """TODO Documentme."""
        return b.rightCell() is None and b.leftCell() is not None

    if bounds is None:
        bounds = []

        for b in mesh.boundaries():
            if isBoundary(b):
                bounds.append(b)

    boundsIdx = []
    for b in bounds:
        boundsIdx.append(b.id())

    poly2 = pg.Mesh(2)
    for p in mesh.cellCenters():
        poly2.createNode(p)

    boundSort = []
    boundSort.append(bounds[0])
    boundSort[-1].tag()
    boundSortIdx = [0]
    lastB = boundSort[-1]
    lastN = lastB.node(1)

    while len(boundSort) != len(bounds):
        newB = None
        newN = None
        for b in lastN.boundSet():
            if not b.tagged() and isBoundary(b):
                boundSort.append(b)
                b.tag()
                newB = b

                if newB.node(0) == lastN:
                    newN = newB.node(1)
                else:
                    newN = newB.node(0)

        lastN = newN
        lastB = newB

        if not newB:
            raise BaseException("implement me")

        boundSortIdx.append(boundsIdx.index(newB.id()))

    bNodes = [poly2.createNode(b.center()) for b in boundSort]

    for i, _ in enumerate(bNodes):
        poly2.createEdge(bNodes[i], bNodes[(i + 1) % len(bNodes)])

    mesh2 = createMesh(poly2, switches='-pezY')
    return mesh2, pg.cat(u, uDirVals[boundSortIdx])
    # def createFVPostProzessMesh(...)


def applyBoundaryValues(uB, mesh, uBBC):
    """TODO Documentme."""
    if isinstance(uBBC, dict):
        boundsDirichlet = pg.solver.parseArgToBoundaries(uBBC['Dirichlet'], mesh)
        for b, v in boundsDirichlet:
            uB[b.id()] = v
    else:
        for marker, val in uBBC:
            for b in mesh.findBoundaryByMarker(marker):
                uB[b.id()] = val


def __d(name, v, showAll=False):
    """TODO Documentme."""
    print(name, np.mean(v), min(v), max(v))
    if showAll:
        print(v)


def __solveStokes(mesh, viscosity, velBoundary=None, preBoundary=None,
                pre0=None, vel0=None,
                tol=1e-4, maxIter=1000,
                verbose=1, **kwargs):
    """Steady Navier-Stokes solver."""
    workspace = pg.solver.WorkSpace()
    wsux = pg.solver.WorkSpace()
    wsuy = pg.solver.WorkSpace()
    wsp = pg.solver.WorkSpace()

    # get cache values if given
    ws = kwargs.pop('ws', None)
    if ws:
        workspace = ws

        if hasattr(workspace, 'wsux'):
            wsux = workspace.wsux
        else:
            workspace.wsux = wsux

        if hasattr(workspace, 'wsuy'):
            wsuy = workspace.wsuy
        else:
            workspace.wsuy = wsuy

        if hasattr(workspace, 'wsp'):
            wsp = workspace.wsp
        else:
            workspace.wsp = wsp

    velocityRelaxation = kwargs.pop('vRelax', 0.5)
    pressureRelaxation = kwargs.pop('pRelax', 0.8)

    pressureCoeff = None
    preCNorm = []
    divVNorm = []

    # velBoundaryX = []
    # velBoundaryY = []

    # if velBoundary is not None:
    #     for marker, vel in velBoundary:
    #         if vel is not None:
    #             if not isinstance(vel[0], str):
    #                 velBoundaryX.append([marker, vel[0]])
    #             if not isinstance(vel[1], str):
    #                 velBoundaryY.append([marker, vel[1]])

    pressure = None
    if pre0 is None:
        pressure = np.zeros(mesh.cellCount())
    else:
        if len(pre0) == mesh.nodeCount():
            pressure = pg.meshtools.nodeDataToCellData(mesh, pre0)
        else:
            pressure = np.array(pre0)

    velocity = None
    if vel0 is None:
        velocity = np.zeros((mesh.cellCount(), mesh.dimension()))
    else:
        velocity = np.array(vel0)

    mesh.createNeighborInfos()

    CtB = mesh.cellToBoundaryInterpolation()
    controlVolumes = CtB * mesh.cellSizes()

    density = kwargs.pop('density', 1.0)
    force = kwargs.pop('f', [0.0, 0.0])
    density = pg.solver.parseArgToArray(density, nDof=mesh.cellCount())
    forceX = pg.solver.parseArgToArray(force[0], nDof=mesh.cellCount())
    forceY = pg.solver.parseArgToArray(force[1], nDof=mesh.cellCount())

    for i in range(maxIter):
        pressureGrad = cellDataToCellGrad(mesh, pressure, CtB)
        # __d('vx', pressureGrad[:,0])

        velocity[:, 0] = solveFiniteVolume(
                mesh, a=viscosity,
                f=-pressureGrad[:, 0] / density + forceX,
                bc=velBoundary[0],
                uL=velocity[:, 0],
                relax=velocityRelaxation,
                ws=wsux)
#        for s in wsux.S:
#            print(s)
#        __d('rhs', wsux.rhs, 1)
#        __d('ux', velocity[:,0])

        velocity[:, 1] = solveFiniteVolume(
                mesh, a=viscosity,
                f=-pressureGrad[:, 1] / density + forceY,
                bc=velBoundary[1],
                uL=velocity[:, 1],
                relax=velocityRelaxation,
                ws=wsuy)

        ap = np.array(wsux.ap * mesh.cellSizes())
        apF = CtB * ap
        uxF = CtB * velocity[:, 0]
        uyF = CtB * velocity[:, 1]

        applyBoundaryValues(uxF, mesh, velBoundary[0])
        applyBoundaryValues(uyF, mesh, velBoundary[1])

        pxF = CtB * pressureGrad[:, 0]
        pyF = CtB * pressureGrad[:, 1]

        pF2 = cellDataToBoundaryGrad(mesh, pressure, pressureGrad)

        velXF = uxF + controlVolumes / apF * (pxF - pF2[:, 0])
        velYF = uyF + controlVolumes / apF * (pyF - pF2[:, 1])

        applyBoundaryValues(velXF, mesh, velBoundary[0])
        applyBoundaryValues(velYF, mesh, velBoundary[1])

        if pressureCoeff is None:
            pressureCoeff = 1. / apF * mesh.boundarySizes() * \
                _boundaryToCellDistances(mesh)

        div = -mesh.divergence(np.vstack([velXF, velYF]).T)

        # boundsDirichlet = pg.solver.parseArgToBoundaries(preBoundary, mesh)
        # pB = CtB * pressure
        # for ix, [boundary, val] in enumerate(boundsDirichlet):
        # boundsDirichlet[ix][1] = -(-pg.solver.generateBoundaryValue(boundary,
        # val) + pB[boundary.id()])

        pressureCorrection = solveFiniteVolume(mesh,
                                               a=pressureCoeff,
                                               f=div,
                                               bc=preBoundary,
                                               # uB=boundsDirichlet,
                                               ws=wsp)

        # pg.show(mesh, pressureCorrection, label='pC')
        # print(i, pg.solver.generateBoundaryValue(mesh.boundary(58), val),
        # pB[58], boundsDirichlet[-1][1],
        # (CtB*pressureCorrection)[58])

        pressure += pressureCorrection * pressureRelaxation

        pressureCorrectionGrad = cellDataToCellGrad(mesh, pressureCorrection,
                                                    CtB)

        velocity[:, 0] -= pressureCorrectionGrad[:, 0] / ap * mesh.cellSizes()
        velocity[:, 1] -= pressureCorrectionGrad[:, 1] / ap * mesh.cellSizes()

        preCNorm.append(pg.core.norm(pressureCorrection))
        divVNorm.append(pg.core.norm(div))

        if workspace:
            workspace.div = div
        #        __d('div', div)
        #        if ( i == 1):
        #            sd

        convergenceTest = 100
        if i > 6:
            convergenceTest = (divVNorm[-1] - divVNorm[-2]) + \
                              (divVNorm[-2] - divVNorm[-3]) + \
                              (divVNorm[-3] - divVNorm[-4]) + \
                              (divVNorm[-4] - divVNorm[-5]) + \
                              (divVNorm[-5] - divVNorm[-6])
            convergenceTest /= 5

        if verbose:
            print("\r" + str(i) + " div V=" + str(divVNorm[-1]) +
                  " ddiv V=" + str(convergenceTest))
        if divVNorm[-1] > 1e6:
            print("\r" + str(i) + " div V=" + str(divVNorm[-1]) +
                  " ddiv V=" + str(convergenceTest))
            raise BaseException("Stokes solver seems diverging")

        if i > 2:
            if i == maxIter or divVNorm[-1] < tol or \
              abs(convergenceTest * divVNorm[-1]) < tol:
                break

        if verbose:
            print(str(i) + ": |pC|" + str(preCNorm[-1]))

    return velocity, pressure, preCNorm, divVNorm


def _test_ConvectionAdvection():
    """Test against a reference solution.
    taken from https://www.ctcms.nist.gov/fipy/examples/flow/generated/examples.flow.stokesCavity.html
    """
    N = 21  # 21 reference
    maxIter = 11  # 11 reference
    Nx = N
    Ny = N

    x = np.linspace(-1.0, 1.0, Nx + 1)
    y = np.linspace(-1.0, 1.0, Ny + 1)
    grid = pg.createGrid(x=x, y=y)

    a = pg.Vector(grid.cellCount(), 1.0)

    b7 = grid.findBoundaryByMarker(1)[0]

    for b in grid.findBoundaryByMarker(1):
        if b.center()[1] < b.center()[1]:
            b7 = b
    b7.setMarker(7)

    swatch = pg.Stopwatch(True)
    # velBoundary = [[1, [0.0, 0.0]],
    #                [2, [0.0, 0.0]],
    #                [3, [1.0, 0.0]],
    #                [4, [0.0, 0.0]],
    #                [7, [0.0, 0.0]]]

    velBoundaryX = {'Dirichlet':{'1,2,4,7': 0.0,
                                3: 1.0,}}
    velBoundaryY = {'Dirichlet':{'1,2,4,7': 0.0,
                                3: 0.0,}}

    # velBoundaryX = {'Dirichlet':{'1,2,7': 0.0,
    #                             3: 1.0,}}
    # velBoundaryY = {'Dirichlet':{'3,4,7': 0.0}}

    preBoundary = {'Dirichlet': {7: 0.0}}

    vel, pres, pCNorm, divVNorm = __solveStokes(grid, a,
                                                [velBoundaryX, velBoundaryY],
                                                preBoundary,
                                                maxIter=maxIter,
                                                verbose=1)

    print("time", len(pCNorm), swatch.duration(True))

    referenceSolutionDivV = 0.029187181920161752 #fipy
    print("divNorm: ", divVNorm[-1])
    print("to reference: ", divVNorm[-1] - referenceSolutionDivV)

    np.testing.assert_approx_equal(divVNorm[-1],
                                   referenceSolutionDivV, significant=8)
    fig = pg.plt.figure()
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    show(grid, data=pg.meshtools.cellDataToNodeData(grid, pres),
         label='Pressure', cMap='RdBu', ax=ax1)
    show(grid, data=pg.core.logTransDropTol(
         pg.meshtools.cellDataToNodeData(grid, vel[:, 0]), 1e-2),
         label='$v_x$', ax=ax2)
    show(grid, data=pg.core.logTransDropTol(
         pg.meshtools.cellDataToNodeData(grid, vel[:, 1]), 1e-2),
         label='$v_y$', ax=ax3)

    show(grid, data=vel, ax=ax1, showMesh=True)

    pg.plt.figure()
    pg.plt.semilogy(pCNorm, label='norm')
    pg.plt.legend()

    pg.wait()

if __name__ == '__main__':

    _test_ConvectionAdvection()
