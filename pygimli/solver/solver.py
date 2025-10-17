#!/usr/bin/env python3
"""Finite-element solver and utility functions."""
from copy import deepcopy

import numpy as np
import pygimli as pg


def parseDictKey_(key, markers):
    """Parse dictionary key of type str to marker list."""
    return parseMarkersDictKey(key, markers)


def parseMarkersDictKey(key, markers):
    """Parse dictionary key of type str to marker list.

    Utility function to parse a dictionary key string into a valid list of
    markers containing in a given markers list.

    Parameters
    ----------
    key: str | int | [int]
        Supported are
        - int: single markers
        - '*': all markers
        - 'm1': Single marker
        - 'm1,m2': Comma separated list
        - ':': Slice wildcard
        - 'start:stop:step': Slice like syntax
        - list of integers (not a dict key but nice for convenience)

    markers: [int]
        List of integers, e.g., cell or boundary markers

    Returns
    -------
    mas: [int]
        List of integers described by key
    """
    markers = pg.unique(markers)
    mas = None

    if isinstance(key, str):
        if key == '*':
            return markers

        if ',' in key:
            mas = [int(k) for k in key.split(',')]
        elif ':' in key:
            sse = key.split(':')

            start = markers[0]
            stop = markers[-1] + 1
            step = 1

            if len(sse) > 0:
                try:
                    start = int(sse[0])
                except BaseException:
                    pass
            if len(sse) > 1:
                try:
                    stop = int(sse[1])
                except BaseException:
                    pass
            if len(sse) > 2:
                try:
                    step = int(sse[2])
                except BaseException:
                    pass

            mas = list(range(start, stop, step))
        else:
            mas = [int(key)]
    elif isinstance(key, (int, float)):
        mas = [int(key)]
    else:
        mas = key

    return [m for m in mas if m in markers]


def boundaryIdsFromDictKey(mesh, key, outside=True):
    """Find all boundaries matching a dictionary key.

    Attributes
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`

    key: str|int
        Representation for boundary marker. Will be parsed by
        :py:mod:`pygimli.solver.solver.parseMarkersDictKey`
    outside: bool [True]
        Only select outside boundaries.

    Returns
    -------
    dict: {marker, [boundary.id()]}
    """
    mas = pg.solver.parseMarkersDictKey(key, mesh.boundaryMarkers())
    ret = {}
    for m in mas:
        for i in pg.find(mesh.boundaryMarkers() == m):
            if m not in ret:
                ret[m] = []
            if outside is True and not mesh.boundary(i).outside():
                continue
            ret[m].append(i)
    return ret


def cellValues(mesh, arg, **kwargs):
    """Get a value for each cell.

    Returns a array or vector of length mesh.cellCount() based on arg.
    The preferable arg is a dictionary for the cell marker and the appropriate
    cell value. The designated value can be calculated using a
    `callable(cell, **kwargs)`, which is called on demand.

    Attributes
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Used if arg is callable

    arg: float | int | complex | ndarray | iterable | callable | dict
        Argument to be parsed as cell data.
        If arg is a dictionary, its key will be interpreted as cell marker:

        Dictionary is {key: value, }. Value can be float, int, complex or
        ndarray. The latter for anistropic or elastic tensors.
        Key can be integer for the cell marker or a str, which will be
        interpreted as wildcard, splice or list. See examples or parseDictKey_.
        String keys will be interpreted before integer keys, the latter will
        overwite the first, i.e. you can fill all with wildcard '*' and
        specify others with integers.
        See examples or `py:mod:pygimli.solver.parseMarkersDictKey`.

        Iterable of length mesh.nodeCount() to be interpolated to cell centers.

    userData: class
        Used if arg is callable

    Returns
    -------
    ret: :gimliapi:`GIMLI::Vector` | ndarray(mesh.cellCount(), xx )
        Array of desired length filled with the appropriate values.

    Examples
    --------
    >>> import pygimli as pg
    >>> mesh = pg.createGrid(x=range(5))
    >>> mesh.setCellMarkers([1, 1, 2, 2])
    >>> print(mesh.cellCount())
    4
    >>> print(pg.solver.cellValues(mesh, [1, 2, 3, 4]))
    [1, 2, 3, 4]
    >>> print(pg.solver.cellValues(mesh, {1:1.0, 2:10}))
    [ 1.  1. 10. 10.]
    >>> print(pg.solver.cellValues(mesh, {':':2.0}))
    [2. 2. 2. 2.]
    >>> print(pg.solver.cellValues(mesh, {'0:2':3.0}))
    [3.0 3.0 None None]
    >>> print(pg.solver.cellValues(mesh, {'*':1.0, '2':2.0}))
    [1. 1. 2. 2.]
    >>> print(pg.solver.cellValues(mesh, {2:2.0, '*':1.0}))
    [1. 1. 2. 2.]
    >>> print(pg.solver.cellValues(mesh, np.ones(mesh.nodeCount())))
    4 [1.0, 1.0, 1.0, 1.0]
    >>> print(np.array(pg.solver.cellValues(mesh, {'1:3':np.diag([1.0, 2.0])})))
    [[[1. 0.]
      [0. 2.]]
    <BLANKLINE>
     [[1. 0.]
      [0. 2.]]
    <BLANKLINE>
     [[1. 0.]
      [0. 2.]]
    <BLANKLINE>
     [[1. 0.]
      [0. 2.]]]
    >>> print(np.array(pg.solver.cellValues(mesh, {':':pg.core.CMatrix(2, 2)})))
    [[[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]
    <BLANKLINE>
     [[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]
    <BLANKLINE>
     [[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]
    <BLANKLINE>
     [[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]]
    >>> print(pg.solver.cellValues(mesh, {'1,2':1 + 1j*2.0}))
    [1.+2.j 1.+2.j 1.+2.j 1.+2.j]
    >>> def cellVal(c, b=1):
    ...     return c.center()[0]*b
    >>> t = pg.solver.cellValues(mesh, {':' : cellVal})
    >>> print([t[c.id()](c) for c in mesh.cells()])
    [0.5, 1.5, 2.5, 3.5]
    """
    if isinstance(arg, dict):

        try:
            val = list(arg.values())[0]
        except BaseException:
            pg.error("Can't interpret empty dictionary:", arg)
            val = 1.0

        ret = [None] * mesh.cellCount()

        ### fill all with wildcard first
        for key, val in arg.items():
            if isinstance(key, str) and key == '*':
                mas = parseDictKey_(key, mesh.cellMarkers())

                for m in mas:
                    for i in pg.find(mesh.cellMarkers() == m):
                        ret[i] = val

        ### fill all string keys
        for key, val in arg.items():
            if isinstance(key, str) and not key == '*':
                mas = parseDictKey_(key, mesh.cellMarkers())

                for m in mas:
                    for i in pg.find(mesh.cellMarkers() == m):
                        ret[i] = val

        ### fill all integer keys
        for key, val in arg.items():
            if isinstance(key, int):
                for i in pg.find(mesh.cellMarkers() == key):
                    ret[i] = val

        return np.array(ret)

    # if arg have already the correct size
    if hasattr(arg, '__len__'):
        if len(arg) == mesh.cellCount():
            return arg
        if len(arg) == mesh.nodeCount():
            return pg.interpolate(mesh, arg, mesh.cellCenters())

    # if arg if scalar or global data type, ndarray or Matrix but not the right
    # size assume global tensor
    if isinstance(arg, (np.ndarray
                        | pg.core.RMatrix
                        | pg.core.CMatrix
                        | float
                        | int
                        | complex)):
        return [arg]*mesh.cellCount()

    return parseArgToArray(arg,
                           nDof=mesh.cellCount(),
                           mesh=mesh, **kwargs)


def parseArgToArray(arg, nDof, mesh=None, userData={}):
    """Parse array related arguments to create a valid value array.

    Parameters
    ----------
    arg: float | int | iterable | callable
        The target array value that will be converted to an array.

        If arg is a callable with it must fulfill:

        :: arg(cell|node|boundary, userData={})

        Where MeshEntity is one of
        :gimliapi:`GIMLI::Cell` ,
        :gimliapi:`GIMLI::Node` or
        :gimliapi:`GIMLI::Boundary`
        depending on nDof, where nDof is mesh.cellCount(),
        mesh.nodeCount() or mesh.boundaryCount(),
        respectively.

    nDof: int | [int]
        Desired array size.

    mesh: :gimliapi:`GIMLI::Mesh`
        Used if arg is callable

    userData: class
        Used if arg is callable

    Returns
    -------
    ret: :gimliapi:`GIMLI::Vector`
        Array of desired length filled with the appropriate values.
    """
    # pg.warn('check if obsolete: parseArgToArray')
    if not hasattr(nDof, '__len__'):
        nDof = [nDof]

    try:
        return pg.Vector(nDof[0], float(arg))
    except BaseException:
        pass

    if hasattr(arg, '__len__'):
        if isinstance(arg, np.ndarray):
            if len(arg) == nDof[0]:
                return arg
            else:
                raise ValueError('Given array does not have requested (' +
                                str(nDof) + ') size (' +
                                str(len(arg)) + ')')

        for n in nDof:
            if len(arg) == n:
                return arg

        try:
            # [marker, val] || [[marker, val]]
            return parseMapToCellArray(arg, mesh)
        except BaseException:
            raise ValueError("Array 'arg' has the wrong size: " +
                            str(len(arg)) + " != " + str(nDof))
    elif callable(arg):
        ret = pg.Vector(nDof[0], 0.0)

        if not mesh:
            raise ValueError("Please provide a mesh for the callable"
                            "argument to parse ")

        if nDof[0] == mesh.nodeCount():
            for n in mesh.nodes():
                if userData:
                    ret[n.id()] = arg(node=n, userData=userData)
                else:
                    ret[n.id()] = arg(node=n)
        elif nDof[0] == mesh.cellCount():
            for c in mesh.cells():
                if userData:
                    ret[c.id()] = arg(cell=c, userData=userData)
                else:
                    ret[c.id()] = arg(cell=c)
        elif nDof[0] == mesh.boundaryCount():
            for b in mesh.boundaries():
                if userData:
                    ret[b.id()] = arg(boundary=b, userData=userData)
                else:
                    ret[b.id()] = arg(boundary=b)
        else:
            raise Exception("Cannot parse callable argument " + str(nDof) +
                            " nodes: " + str(mesh.nodeCount()) +
                            " cells: " + str(mesh.cellCount()))

        return ret
    raise Exception("Cannot parse argument type " + str(type(arg)))


def generateBoundaryValue(boundary, arg, time=0.0, userData={},
                          expectList=False, nCoeff=1):
    """Generate a value for the given Boundary.

    TODO
    ----
        * support for complex vals

    Parameters
    ----------
    boundary: :gimliapi:`GIMLI::Boundary` or list of ..
        The related boundary.
    expectList: bool[False]
        Allow list values for Robin BC.
    arg: convertible | iterable | callable or list of ..
        - convertible into float
        - iterable of minimum length = boundary.id()
        - callable generator function

        If arg is a callable it must fulfill:

        :: arg(boundary=:gimliapi:`GIMLI::Boundary`, time=0.0, userData={})

        The callable function arg have to return appropriate values for all
        nodes of the boundary or one value for all nodes (scalar field only).
        Value can be scalar or vector field value, e.g., return force values
        for all nodes at a boundary to return an ndarray((nodes, dims)), e.g.
        'lambda _b: np.array([[forc_x, forc_y, forc_z] for n in _b.nodes()]).T'

    Returns
    -------
    val: [float]
        Value for all nodes of the boundary.
    """
    # pg._r('genBounds')
    val = 0.

    if callable(arg):
        kwargs = dict()
        # pg._r(arg)
        if time != 0.0 and time is not None:
            kwargs['time'] = time
        if userData is not None and userData.keys():
            kwargs['userData'] = userData
        try:
            # val(boundary, time=0, userData={})
            val = arg(boundary, **kwargs)

        except BaseException as e:
            print(arg, "(", kwargs, ")")
            pg.critical("Wrong arguments for callback function.", e)

    elif hasattr(arg, '__len__'):
        if callable(arg[0]):
            kwargs = arg[1]
            if time != 0.0 and time is not None:
                kwargs['time'] = time
            val = arg[0](boundary=boundary, **kwargs)
        else:
            val = arg
    else:
        try:
            val = float(arg)
        except ValueError:
            print(arg, val)
            pg.error("can't create boundary values.")

    # transform val into list of length nodeCount

    if expectList is True and np.array(val).ndim != 2:
        val = np.atleast_1d(val)

    if isinstance(boundary, pg.core.Node):
        return val

    if nCoeff == 1 and expectList is False:
        if isinstance(val, float):
            val = np.ones(boundary.nodeCount(), dtype=float) * val
        if len(val) != boundary.nodeCount():
            print(val)
            pg.critical("Boundary value cannot be generated for nCoeff=1 val:",
                        val)
    else:
        val = np.atleast_2d(val)
        # pg._y(val)
        if len(val) != boundary.nodeCount() or val.shape[1] != nCoeff:
            val = np.tile(val, (boundary.nodeCount(), 1))

    return val


def parseArgPairToBoundaryArray(pair, mesh):
    r"""Parse boundary-related argument pairs.

    Creates a list of [:gimliapi:`GIMLI::Boundary`, value|callable ].

    Parameters
    ----------
    pair: tuple
        - [marker, arg]
        - [marker, [callable, `**`kwargs]]
        - [marker, [arg_x, arg_y, arg_z]]
        - [boundary, arg]
        - [`*`, arg]
        - [node, arg]

        arg will be parsed by
        :py:mod:`pygimli.solver.solver.generateBoundaryValue`
        and distributed to each boundary.
        Callable functions will be executed at run time.
        wildcard `*` is interpreted as all boundary elements with one neighboring cell
    mesh: :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker.

    Returns
    -------
    bc: list()
        [:gimliapi:`GIMLI::Boundary`, value|callable]
    """
    bc = []
    bounds = []
    #  [marker, [callable, *kwargs]]
    if isinstance(pair[1], list) and callable(pair[1][0]):
        pair = [pair[0]] + pair[1]

    if pair[0] == '*':
        mesh.createNeighborInfos()
        for b in mesh.boundaries():
            if b.leftCell() is not None and b.rightCell() is None:
                bounds.append(b)

    elif isinstance(pair[0], int):
        bounds = mesh.findBoundaryByMarker(pair[0])
    elif isinstance(pair[0], pg.core.Node):
        bc.append(pair)
        return bc

    for b in bounds:
        val = None
        val = pair[1:] if len(pair) > 2 else pair[1]

        bc.append([b, val])

        # print('-'*50)
        # print(b, pair[1], callable(pair[1]))
        # print('+'*50)
        # if callable(pair[1]):
        #     # don't execute the callable here
        #     # we want to call them at runtime
        #     if len(pair) > 2:
        #         val = pair[1:]
        #     else:
        #         val = pair[1]
        # else:
        #     this will be executed
        #     val = generateBoundaryValue(b, pair[1])

    # print('#'*30)
    return bc


def parseArgToBoundaries(args, mesh):
    """Parse boundary related arguments to create a valid boundary value list.

    List is of type  `[ :gimliapi:`GIMLI::Boundary`, value|callable ]`

    TODO
    ----
    - callable dynamic at runtime

    Parameters
    ----------
    args : dict, float, callable
        Dictionary is preferred (key=value|callable).
        If args is just a callable or float every outer boundary is processed
        with args.

        List pairs will be removed or not correct parsed for vector valued
        problems. Callable will be evaluated at runtime. See examples.
        Else see :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray`

    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker

    Returns
    -------
    boundaries : list()
        [ :gimliapi:`GIMLI::Boundary`, value|callable ]

    Examples
    --------
    >>> # no need to import matplotlib. pygimli show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> plc = mt.createWorld([0, 0], [1, -1], worldMarker=0)
    >>> ax, _ = pg.show(plc, boundaryMarker=True)
    >>> mesh = mt.createMesh(plc)
    >>> # all four outer boundaries get value = 1.0
    >>> b = pg.solver.parseArgToBoundaries(1.0, mesh)
    >>> print(len(b))
    4
    >>> # all edges with marker 1 get value = 1.0
    >>> b = pg.solver.parseArgToBoundaries({1: 1.0}, mesh)
    >>> print(len(b))
    1
    >>> # same as above with marker 2 get value 2
    >>> b = pg.solver.parseArgToBoundaries({'1': 1.0, 2 : 2.0}, mesh)
    >>> print(len(b))
    2
    >>> # same as above with marker 3 get value 3
    >>> b = pg.solver.parseArgToBoundaries({1:1., 2:2., 3:3.}, mesh)
    >>> print(len(b))
    3
    >>> # Boundary values for vector valued problem
    >>> b = pg.solver.parseArgToBoundaries({1:[1.0, 1.0]}, mesh)
    >>> print(len(b), b[0][1])
    1 [1.0, 1.0]
    >>> # edges with marker 1 and 2 get value 1
    >>> b = pg.solver.parseArgToBoundaries({'1,2':1.0}, mesh)
    >>> print(len(b))
    2
    >>> b = pg.solver.parseArgToBoundaries({'1, 2, 3': 1.0}, mesh)
    >>> print(len(b))
    3
    >>> b = pg.solver.parseArgToBoundaries({'1:4':1.0, 4:4.0}, mesh)
    >>> print(len(b))
    4
    >>> b = pg.solver.parseArgToBoundaries({mesh.node(0):0.0}, mesh)
    >>> print(len(b))
    1
    >>> def bCall(boundary):
    ...     u = []
    ...     for i, n in enumerate(boundary.nodes()):
    ...         u.append(i)
    ...     return u
    >>> b = pg.solver.parseArgToBoundaries({1:bCall}, mesh)
    >>> print(len(b),b[0][1](b[0][0]))
    1 [0, 1]
    >>> def bCall(boundary, a1, a2):
    ...     return a1 + a2
    >>> b = pg.solver.parseArgToBoundaries({1: [bCall, {'a1':2, 'a2':3}]}, mesh)
    >>> print(len(b), b[0][1][0](b[0][0], **b[0][1][1]))
    1 5
    >>> b = pg.solver.parseArgToBoundaries({1: [bCall, {'a1':1, 'a2':2}],
    ...                                     2: [bCall, {'a1':3, 'a2':4}]}, mesh)
    >>> print(len(b), b[0][1][0](b[0][0], **b[0][1][1]))
    2 3
    >>> b = pg.solver.parseArgToBoundaries({'1,2': [bCall, {'a1':4, 'a2':5}]}, mesh)
    >>> print(len(b), b[1][1][0](b[1][0], **b[1][1][1]))
    2 9
    >>> pg.wait()
    """
    boundaries = list()

    if isinstance(args, dict):

        # try:
        #     val = list(args.values())[0]
        # except BaseException as _:
        #     return boundaries
        #     pg.error("Can't interpret empty dictionary:", args)

        for key, val in args.items():
            if isinstance(key, pg.core.stdVectorNodes) or \
               isinstance(key, list) and isinstance(key[0], pg.core.Node):
                for n in key:
                    boundaries += parseArgPairToBoundaryArray([n, val], mesh)

            elif isinstance(key, str) and key != '*':
                markers = parseDictKey_(key, mesh.boundaryMarkers())

                for m in markers:
                    boundaries += parseArgPairToBoundaryArray([m, val], mesh)
            else:
                boundaries += parseArgPairToBoundaryArray([key, val], mesh)

        return boundaries

    if callable(args) or isinstance(args, (float, int)):
        return parseArgToBoundaries({'*': args}, mesh)

    else:
        raise Exception('cannot interpret boundary token', args)

    return boundaries


def _bcIsForVectorValues(bc, mesh):
    """Guess if boundary condition is for vector-valued problems."""
    verbose = False

    def testForV3(t):
        if verbose:
            print("test for v3", t)
        try:
            test = t(mesh.boundary(0)) if callable(t) else t

            if verbose:
                print("test for v3 test", test)

            if hasattr(test, '__iter__'):
                test = np.array(test)
                # call(b): [v_i] in R with i==1..nodeCount() -> scalar values
                # call(b): [v_i] in R³ with i==1..nodeCount() -> value values
                if len(test) == mesh.boundary(0).nodeCount() and \
                    len(test[0]) == mesh.dim():
                        return True

                if test.ndim == 2 and len(test[0]) == mesh.dim():
                    return True
        except BaseException as e:
            if verbose:
                print(e)
        return False

    for key, _bVal in bc.items():
        if key == 'Robin':
            continue

        if verbose:
            print("test vector values for:", _bVal)

        if isinstance(_bVal, list):
            if isinstance(_bVal[0], list):
                # [[nodeID, [x, y, z]], [nodeID, [x, y, z]]]
                if testForV3(_bVal[0][1]):
                    return True
            else:
                # [nodeID, [x, y, z]]
                if testForV3(_bVal[1]):
                    return True

        elif isinstance(_bVal, dict):
            # {key, bc}
            for key, test in _bVal.items():
                if testForV3(test):
                    return True

        else:
            # [x, y, z] = call(boundary)
            if testForV3(_bVal):
                return True

    return False


def parseMapToCellArray(attributeMap, mesh, default=0.0):
    """Parse a value map to cell attributes.

    A map should consist of pairs of marker and value.
    A marker is an integer and corresponds to the cell.marker().

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        For each cell of mesh a value will be returned.

    attributeMap : list | dict
        List of pairs [marker, value] ] || [[marker, value]],
        or dictionary with marker keys

    default : float [0.0]
        Fill all unmapped attributes to this default.

    Returns
    -------
    att : array
        Array of length mesh.cellCount()
    """
    # pg.warn('check if obsolete: parseMapToCellArray')
    att = pg.Vector(mesh.cellCount(), default)

    if isinstance(attributeMap, dict):
        for marker, value in attributeMap.items():
            idx = pg.find(mesh.cellMarkers() == marker)
            if len(idx) == 0:
                pg.warn("parseMapToCellArray: cannot find marker " +
                        str(marker) + " within mesh.")
            else:
                if isinstance(value, complex):
                    if not isinstance(att, pg.CVector):
                        att = pg.math.toComplex(att)
                    att.setVal(val=value, ids=idx)
                else:
                    att.setVal(val=float(value), ids=idx)
    elif hasattr(attributeMap, '__len__'):
        if not hasattr(attributeMap[0], '__len__'):
            # assuming [marker, value]
            attributeMap = [attributeMap]

        for pair in attributeMap:
            if hasattr(pair, '__len__'):
                idx = pg.find(mesh.cellMarkers() == pair[0])
                if len(idx) == 0:
                    pg.warn("parseMapToCellArray: cannot find marker " +
                            str(pair[0]) + " within mesh.")
                else:
                    # print('---------------------')
                    # print(att, idx, pair[1], type(pair[1]), float(pair[1]))
                    if isinstance(pair[1], complex):
                        if not isinstance(att, pg.CVector):
                            att = pg.math.toComplex(att)
                        att.setVal(val=pair[1], ids=idx)
                    else:
                        att.setVal(val=float(pair[1]), ids=idx)
            else:
                raise Exception("Please provide a list of [int, value] pairs" +
                                str(pair))
    else:
        print("attributeMap: ", attributeMap)
        raise Exception("Cannot interpret attributeMap!")

    return att


def grad(mesh, u, r=None):
    r"""Return gradient for scalar field.

    Return the discrete interpolated gradient :math:`\mathbf{v}`
    for a given scalar field :math:`\mathbf{u}`.

    .. math::
        \mathbf{v}(\mathbf{r}_{\mathcal{C}})
        &=
        \nabla u(\mathbf{r}_{\mathcal{N}})
        \\
        (\mathbf{v_x}(\mathbf{r}_{\mathcal{C}}),
         \mathbf{v_y}(\mathbf{r}_{\mathcal{C}}),
         \mathbf{v_z}(\mathbf{r}_{\mathcal{C}}))^{\text{T}}
        &=
        \left(\frac{\partial u}{\partial x},
         \frac{\partial u}{\partial y},
         \frac{\partial u}{\partial z}\right)^{\text{T}}

    With :math:`\mathcal{N}=\cup_{i=0}^{N} \text{Node}_i`,
    :math:`\mathcal{C}=\cup_{j=0}^{M} \text{Cell}_j`,
    :math:`\mathbf{u}=\{u(\mathbf{r}_{i})\}\in I\!R` and
    :math:`\mathbf{r}_{i} = (x_i, y_i, z_i)^{\text{T}}`

    The discrete scalar field
    :math:`\mathbf{u}(\mathbf{r}_{\mathcal{N}})`
    need to be defined for each node position :math:`\mathbf{r}_{\mathcal{N}}`.
    The resulting vector field :math:`\mathbf{v}(\mathbf{r}_{\mathcal{C}})`
    is defined for each cell center position :math:`\mathbf{r}_{\mathcal{C}}`.
    If you need other positions than the cell center,
    provide an appropriate array of coordinates :math:`\mathbf{r}`.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Discretization base, interpolation will be performed via finite element
        base shape functions.
    u: array | callable
        Scalar field per mesh node position or an appropriate
        `callable([[x,y,z]])`
    r: ndarray((M, 3)) [mesh.cellCenter()]
        Alternative target coordinates :math:`\mathbf{r}` for the resulting
        gradient field. i.e., the positions where the vector field is defined.
        Default are all cell centers.

    Returns
    -------
    v: ndarray((M, 3))
        Resulting vector field defined on
        :math:`\mathbf{v}(\mathbf{r}_{\mathcal{C}})`.
        M is number of cells or length of given alternative coordinates r.

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> fig, ax = pg.plt.subplots()
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 20), y=np.linspace(0, 1, 20))
    >>> u = lambda p: pg.x(p)**2 * pg.y(p)
    >>> pg.show(mesh, u(mesh.positions()), ax=ax) # doctest: +ELLIPSIS
    (<Axes...
    >>> pg.show(mesh, [2.*pg.y(mesh.cellCenters())*pg.x(mesh.cellCenters()),
    ...              pg.x(mesh.cellCenters())**2], ax=ax)  # doctest: +ELLIPSIS
    (<Axes...
    >>> pg.show(mesh, pg.solver.grad(mesh, u), ax=ax, color='w',
    ...        linewidth=0.4) # doctest: +ELLIPSIS
    (<Axes...
    >>> pg.plt.show()
    """
    if r is None:
        r = mesh.cellCenters()

    uv = u
    if callable(u) and not isinstance(u, pg.Vector):
        uv = u(mesh.positions())

    if len(uv) == mesh.cellCount():
        uv = pg.meshtools.cellDataToNodeData(mesh, uv)

    v = np.ndarray((len(r), 3))

    for i, _ in enumerate(v):
        c = mesh.findCell(r[i])
        if c:
            v[i] = c.grad(r[i], uv)

    return v


def div(mesh, v):
    r"""Return the discrete interpolated divergence field.

        Return the discrete interpolated divergence field. :math:`\mathbf{u}`
        for each cell for a given vector field :math:`\mathbf{v}`.
        First order integration via boundary center.

    .. math::
        d(cells) & = \nabla\cdot\vec{v} \\
        d(c_i) & = \sum_{j=0}^{N_B}\vec{v}_{B_j} \cdot \vec{n}_{B_j}

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Discretization base, interpolation will be performed via finite element
        base shape functions.

    V : array(N,3) | R3Vector
        Vector field at cell centers or boundary centers

    Returns
    -------
    d : array(M)
        Array of divergence values for each cell in the given mesh.

    Examples
    --------
    >>> import pygimli as pg
    >>> import numpy as np
    >>> v = lambda p: p
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4))
    >>> print(pg.math.round(pg.solver.div(mesh, v(mesh.boundaryCenters())), 1e-5))
    3 [1.0, 1.0, 1.0]
    >>> print(pg.math.round(pg.solver.div(mesh, v(mesh.cellCenters())), 1e-5))
    3 [0.5, 1.0, 0.5]
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4),
    ...                      y=np.linspace(0, 1, 4))
    >>> print(pg.math.round(pg.solver.div(mesh, v(mesh.boundaryCenters())), 1e-5))
    9 [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    >>> divCells = pg.solver.div(mesh, v(mesh.cellCenters()))
    >>> # divergence from boundary values are exact where the divergence from
    >>> # interpolated cell center values wrong due to interpolation to boundary
    >>> print(np.round(sum(divCells),12))
    12.0
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4),
    ...                      y=np.linspace(0, 1, 4),
    ...                      z=np.linspace(0, 1, 4))
    >>> print(sum(pg.solver.div(mesh, v(mesh.boundaryCenters()))))
    81.0
    >>> divCells = pg.solver.div(mesh, v(mesh.cellCenters()))
    >>> print(np.round(sum(divCells),12))
    54.0
    """
    mesh.createNeighborInfos()
    d = None
    if hasattr(v, '__len__'):
        if len(v) == mesh.boundaryCount():
            d = mesh.divergence(v)
        elif len(v) == mesh.nodeCount():
            d = mesh.divergence(pg.meshtools.nodeDataToBoundaryData(mesh, v))
        elif len(v) == mesh.cellCount():
            CtB = mesh.cellToBoundaryInterpolation()
            d = mesh.divergence(np.array([CtB*pg.x(v),
                                          CtB*pg.y(v),
                                          CtB*pg.z(v)]).T)
        else:
            print(len(v), mesh)
            raise BaseException("implement me")
    elif callable(v):
        raise BaseException("implement me")

    return d


def divergence(mesh, func=None, normMap=None, order=1):
    """Divergence for callable function func((x,y,z)).

    MOVE THIS to a better place

    Return sum div over boundary.

    Parameters
    ----------
    func : function
        function to be called
    normMap : [[]]
        mapping

    Returns
    -------
    div : iterable
    """
    if func is None:
        func = lambda r: r

    diverg = 0
    directionCheck = False

    if mesh.cellCount() > 0:
        directionCheck = True

    bNorms = None
    if normMap is not None:
        bNorms = np.zeros((mesh.boundaryCount(), 2))
        for pair in normMap:
            bounds = mesh.findBoundaryByMarker(pair[0])
            for b in bounds:
                bN = [0.0, 0.0]
                if not isinstance(pair[1][0], str):
                    bN[0] = pair[1][0]
                if not isinstance(pair[1][1], str):
                    bN[1] = pair[1][1]

                bNorms[b.id()] = bN

    for b in mesh.boundaries():

        if directionCheck:
            if b.leftCell() is None and b.rightCell() is None:
                # print(b.id(), b.leftCell(), b.rightCell())
                sw = pg.Stopwatch(True)
                mesh.createNeighborInfos()
                print("NeighborInfos()", sw.duration(True))
                # return gauss(grid, F)

            # don't calc for inner boundaries
            if b.leftCell() is not None and b.rightCell() is not None:
                continue

        divS = 0
        shape = b.shape()

        if order == 1:
            if bNorms is not None:
                divS = shape.norm().dot(bNorms[b.id()]) * shape.domainSize()
            else:
                divS = shape.norm().dot(
                    func(shape.center())) * shape.domainSize()
        else:
            weights = pg.core.IntegrationRules.instance().weights(shape, order)
            abscissa = pg.core.IntegrationRules.instance().abscissa(shape,
                                                                    order)

            for i, p in enumerate(abscissa):
                rPos = shape.xyz(p)
                divS += shape.norm().dot(func(rPos)) * \
                    weights[i] * shape.domainSize()

        if directionCheck and b.leftCell() is None:
            divS *= -1
            # raise Exception("invalid mesh: left is None .. every
            # boundary need leftCell")

        diverg += divS

    return diverg

def triDiagToeplitz(dom, a, l, r, start=0, end=-1):
    """Create tri-diagonal Toeplitz matrix."""
    A = pg.matrix.SparseMapMatrix(dom, dom)

    if end == -1:
        end = dom

    for i in range(start, end):
        A.addVal(i, i, a)
        if i > start:
            A.addVal(i, i - 1, l)

        if i < end - 1:
            A.addVal(i, i + 1, r)
    return A


def identity(dom, start=0, end=-1, scale=1):
    """Create identity matrix."""
    A = pg.matrix.SparseMapMatrix(dom, dom)

    if end == -1:
        end = dom

    for i in range(start, end):
        if hasattr(scale, '__len__'):
            A.addVal(i, i, scale[i])
        else:
            A.addVal(i, i, scale)
    return A


def showSparseMatrix(mat, full=False):
    """Show the content of a sparse matrix.

    Parameters
    ----------
    mat: :gimliapi:`GIMLI::SparseMatrix` | :gimliapi:`GIMLI::SparseMapMatrix`
        Matrix to be shown.
    full: bool [False]
        Show as dense matrix.
    """
    if isinstance(mat, pg.matrix.SparseMapMatrix):
        m_ = pg.matrix.SparseMatrix(mat)
        return showSparseMatrix(m_, full)
    else:
        rows = mat.vecRowIdx()
        cols = mat.vecColPtr()
        vals = mat.vecVals()

        matD = None
        if full is True:
            matD = pg.Matrix(mat.rows(), mat.cols())

        for i in range(mat.rows()):
            for j in range(cols[i], cols[i + 1]):
                if full:
                    matD[i, rows[j]] = vals[j]
                else:
                    if vals[j] != 0:
                        print(i, rows[j], vals[j])

        if full is True:
            print(np.array(matD))


class LinSolver:
    """Proxy class for the solution of linear systems of equations."""
    def __init__(self, mat=None, solver=None, verbose=False, **kwargs):
        """Init the solver proxy class with a matrix and start factorization.

        Args
        ----
        solver: str [None]
            Name for the used solver (pg (umfpack or cholmod), scipy).
            If solver is none decide from matrix type.
        """
        self.verbose = verbose

        self._m = None ## hold local copy if we need to convert the matrix first
        self._x = None ## hold local copy for result
        self._solver = None ## (str) 'pg', 'scipy', ...
        self.factorTime = 0.0
        self.solvingTime = 0.0
        self.solverStr= ''
        self._factorized = False
        self._desiredArrayType = np.array

        # callback for solving
        self._solveCB = None

        if solver is None:
            if isinstance(mat, pg.matrix.MatrixBase):
                solver = 'PG'
            elif isinstance(mat, np.ndarray):
                solver = 'numpy'
                pg.critical("Not yet implemented!")
            else:
                from scipy.sparse import spmatrix
                if isinstance(mat, spmatrix):
                    solver = 'SciPy'
                else:
                    print(mat)
                    pg.critical("Data format for matrix not recognized.")

        if solver.lower() == 'pg':
            self.solverStr = 'PG'
        elif solver.lower() == 'scipy':
            self.solverStr = 'SciPy'
        else:
            pg.critical('inuse?')
            self.solverStr = solver

        self._solveCB = getattr(self, f'solve_{self.solverStr}')
        self._factorizeCB = getattr(self, f'factorize_{self.solverStr}')

        if self.verbose:
            pg.info(f"Solving with {self.solver}")

        if mat is not None:
            self.factorize(mat)

    def isFactorized(self):
        """Is factorized."""
        return self._factorized

    def factorize(self, mat):
        """Factorize matrix mat."""
        pg.tic(key='LinSolver.factorize')

        self._factorizeCB(mat)
        self._factorized = True

        self.factorTime = pg.dur(key='LinSolver.factorize', reset=True)

        if self.verbose:
            pg.info("Matrix factorization:", self.factorTime)

    def solve(self, b):
        """ """
        pg.tic(key='LinSolver.solve')

        x = self._solveCB(b)

        self.solverTime = pg.dur(key='LinSolver.solve', reset=True)
        if self.verbose:
            pg.info("Matrix solve:", self.solverTime)
        return x

    def factorize_PG(self, mat):
        """"""
        self._m = pg.matrix.toSparseMatrix(mat)
        self._desiredArrayType = pg.Vector
        self._solver = pg.core.LinSolver(self._m, verbose=self.verbose)

    def solve_PG(self, b):
        """"""
        self._x = pg.Vector()
        self._solver.solve(self._convertRHS(b), self._x)
        return self._x

    def factorize_SciPy(self, mat):
        """"""
        self._m = pg.matrix.asCSC(mat)

        scipy = pg.optImport('scipy', 'Used for sparse linear solver.')

        from scipy.sparse.linalg import factorized

        self._desiredArrayType = np.array
        self._solver = factorized(self._m)

    def solve_SciPy(self, b):
        """"""
        self._x = self._solver(self._convertRHS(b))
        return self._x

    def __call__(self, b):
        """Shortcut to self.solve(b)."""
        return self.solve(b)

    def _convertRHS(self, b):
        """Convert right hand side vector into the desired format."""
        if not isinstance(b, type(self._desiredArrayType(0))):
            return self._desiredArrayType(b)
        return b


def linSolve(mat, b, solver=None, verbose=False, **kwargs):
    r"""Direct linear solution after :math:`\textbf{x}` using core LinSolver.

    .. math::

        \textbf{A}\textbf{x} = \textbf{b}

    If :math:`\textbf{A}` is symmetric, sparse and positive definite.

    Parameters
    ----------
    mat: :gimliapi:`GIMLI::SparseMatrix`, :gimliapi:`GIMLI::SparseMapMatrix`,
        numpy.array

        System matrix. Need to be symmetric, sparse and positive definite.

    b: iterable array
        Right hand side of the equation.

    solver: str [None]
        Try to choose a solver, 'pg' for pygimli core cholmod or umfpack.
        'np' for numpy linalg or scipy.sparse.linalg.
        Automatic choosing if solver is None depending on matrixtype.

    verbose: bool [False]
        Be verbose.

    Returns
    -------
    x: :gimliapi:`GIMLI::Vector`
        Solution vector.
    """
    # TODO!! refactor with LinSolver
    swatch = pg.Stopwatch()
    reorder = kwargs.pop('reorder', False)
    # perm = None

    # determine the solver if none set
    if solver is None:
        if isinstance(mat, pg.matrix.MatrixBase):
            solver = 'pg'
        elif isinstance(mat, np.ndarray):
            solver = 'numpy'
        else:
            from scipy.sparse import spmatrix
            if isinstance(mat, spmatrix):
                solver = 'scipy'
            else:
                print(mat)
                pg.critical("Data format for matrix not recognized.")

    if solver.lower() == 'pg':
        # core proxy to cholmod and LDL for float and umfpack for complex
        if reorder is True:
            pg.warning('Matrix reordering for pg core solver not yet implemented')
        _m = pg.matrix.asSparseMatrix(mat)

        pg._b(verbose)
        solver = pg.core.LinSolver(_m, verbose=verbose)

        if verbose:
            pg.info(f"Solving with {solver.solverName()}")
            pg.info("Matrix factorization:", swatch.duration(restart=True))

        x = solver.solve(b)

        if verbose:
            pg.info("Matrix solution:", swatch.duration())

    elif solver.lower() == 'numpy':
        if verbose:
            pg.info("Solving with np.linalg.solve")

        x = np.linalg.solve(mat, b)

    elif solver.lower() == 'scipy':
        # pg._r(swatch.duration(restart=True))
        _m = pg.matrix.asCSR(mat)
        # pg._r('convert', swatch.duration(restart=True))

        # scipy is now a dependency
        # scipy = pg.optImport('scipy', 'Used for sparse linear solver.')
        # pg._r('import', swatch.duration(restart=True))

        from scipy.sparse.linalg import spsolve

        if verbose:
            pg.info("Solving with scipy.sparse.spsolve")

        if reorder is True and 0:  # deactivated
            def permCOO(M, perm):
                # M.indices = perm.take(M.indices)
                # M = M.tocsc()
                # M.indices = perm.take(M.indices)
                # return M.tocsr()

                return M[np.ix_(perm, perm)]

                # print(M.row.shape, M.col.shape)
                # rowP = perm[M.row]
                # colP = perm[M.col]
                # print(rowP.shape, colP.shape)
                # MP = scipy.sparse.coo_matrix((M.data, (rowP, colP)),
                #                              shape=M.shape)
                # return MP

            # perm = scipy.sparse.csgraph.reverse_cuthill_mckee(_m)
            # pg._r('reverse_cuthill_mckee', swatch.duration(restart=True))

            # ax, _ = pg.show(_m)
            # _m = permCOO(_m, perm)
            # _m = permCOO(_m.tocoo(), perm).tocsr()
            # pg.show(_m, ax=ax, color='green')

            # pg._r('perm matrix', swatch.duration(restart=True))

            # x = spsolve(_m, b.array())[perm]

            # x = spsolve(_m, b.array()[perm])#[perm]
            # x = x[perm]
        else:
            x = spsolve(_m, b)
            if verbose:
                print(f'spsolve: norm(b-Ax) = {np.linalg.norm(b-_m.dot(x))}')

        # pg._r(swatch.duration())
    return x


def applyDirichlet(mat, rhs, uDirIndex, uDirichlet):
    """Apply Dirichlet BC. This should be moved into the core."""
    # idx = np.argsort(uDirIndex)
    if mat is not None:
        if rhs is not None:
            uDir = pg.Vector(mat.rows(), 0.0)
            uDir.setVal(uDirichlet, uDirIndex)
            rhs -= mat * uDir

        for i in uDirIndex:
            mat.cleanRow(i)
            mat.cleanCol(i)
            mat.setVal(i, i, 1.0)

    if rhs is not None:
        rhs[uDirIndex] = uDirichlet
        # rhs.setVal(uDirichlet, uDirIndex)


def getDirichletMap(mat, boundaryPairs, time=0.0, userData={},
                    nodePairs=None, dofOffset=0, nCoeff=1, dofPerCoeff=None):
    r"""Get map of index: dirichlet value.

    Apply Dirichlet boundary condition to the system matrix S and rhs vector.
    The right hand side values for h can be given for each boundary
    element individually by setting proper boundary pair arguments.

    .. math::
        u(\textbf{r}, t) = h
        \quad\text{for}\quad\textbf{r}\quad\text{on}\quad\delta\Omega=
        \Gamma_{\text{Dirichlet}}

    Parameters
    ----------
    mat: :gimliapi:`GIMLI::SparseMatrix`
        System matrix of the system equation.
    boundaryPair: list()
        List of pairs [:gimliapi:`GIMLI::Boundary`, h].
        The value :math:`h` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`h` need to be a scalar value (float or int) or
        a value generator function that will be executed at runtime.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,
    nodePairs: list() | callable
        List of pairs [nodeID, uD].
        The value uD will assigned to the nodes given there ids.
        This node value settings will overwrite any prior settings due to
        boundaryPair.
    time: float
        Will be forwarded to value generator.
    userData: class
        Will be forwarded to value generator.
    dofOffset: int[0]
        Offset for matrix index.
    """
    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    # uDirNodes = []   # []
    uDirVal = dict()  # {nID: val}

    def _genVecUd(n, ud, dofOffset, nCoeff=1, dofPerCoeff=None):
        ret = {}
        if callable(ud):
            pg.error("callable node pairs need to be implemented.")

        if isinstance(n, pg.core.Node):
            idx = dofOffset + n.id()
        else:
            idx = dofOffset + n

        if hasattr(ud, '__iter__'):
            # vector valued problem
            if dofPerCoeff is None:
                if mat.shape[0] % len(ud) != 0:
                    print(mat)
                    print(mat.shape, len(ud))
                    pg.error("Matrix size missmatch for vector valued problem")
                else:
                    dofPerCoeff = mat.shape[0] // len(ud)

            if nCoeff == 1:
                nCoeff = len(ud)

            for i in range(nCoeff):
                if ud[i] is not None:
                    ret[idx + i * dofPerCoeff] = ud[i]
        else:
            if nCoeff > 1:
                print('nCoeff:', nCoeff, 'ud:', ud, 'idx:', idx)
                pg.error('number of coefficents > 1 but uDirichlet is scalar.')

            if ud is not None:
                ret[idx] = ud
        return ret

    for pair in boundaryPairs:
        ent = pair[0]
        val = pair[1]
        # print('**', ent, val)
        uD = generateBoundaryValue(ent, val, time, userData, nCoeff=nCoeff)
        # print('\t', uD)

        if uD is not None:

            if isinstance(ent, pg.core.Node):
                uDirVal.update(_genVecUd(ent, uD, dofOffset))
            else:
                if isinstance(uD, float):
                    pg.critical(uD)
                    uD = [uD] * ent.nodeCount()
                if len(uD) == ent.nodeCount():
                    # print('uD', uD, nCoeff, dofPerCoeff)
                    for i, n in enumerate(ent.nodes()):
                        uDirVal.update(_genVecUd(n, uD[i], dofOffset,
                                                 nCoeff=nCoeff,
                                                 dofPerCoeff=dofPerCoeff))
                else:
                    pg.error('Dirichlet values per boundary need to have '
                             'length of boundary.nodeCount()')

    if nodePairs is not None:
        # print("nodePairs", nodePairs)

        if len(nodePairs) == 2 and isinstance(nodePairs[0], int):
            # assume a single Node [NodeId, val]
            nodePairs = [nodePairs]

        for [n, val] in nodePairs:
            uDirVal.update(_genVecUd(n, val, dofOffset,
                           nCoeff=nCoeff, dofPerCoeff=dofPerCoeff))

    return uDirVal


def assembleDirichletBC(mat, boundaryPairs, rhs=None, time=0.0, userData={},
                        nodePairs=None,
                        dofOffset=0, nCoeff=1, dofPerCoeff=None):
    r"""Apply Dirichlet boundary condition.

    Args
    ----
    rhs: :py:mod:`Vector`
        Right hand side vector of the system equation will bet set to
        :math:`u_{\text{D}}`
    """
    uDirVal = getDirichletMap(mat, boundaryPairs, time=time,
                              userData=userData,
                              nodePairs=nodePairs,
                              dofOffset=dofOffset,
                              nCoeff=nCoeff,
                              dofPerCoeff=dofPerCoeff)

    # pg._g(list(uDirVal.keys()), list(uDirVal.values()))
    if not uDirVal.keys():
        return

    applyDirichlet(mat, rhs, list(uDirVal.keys()), list(uDirVal.values()))
    return uDirVal


def assembleNeumannBC(rhs, boundaryPairs, nDim=1, time=0.0, userData={},
                      dofOffset=0, nCoeff=1, dofPerCoeff=None):
    r"""Apply Neumann condition to the system matrix S.

    Apply Neumann condition to the system matrix S.
    The right hand side values for g can be given for each boundary
    element individually by setting proper boundary pair arguments.

    .. math::
        \frac{\partial u(\textbf{r}, t)}{\partial\textbf{n}}
        = \textbf{n}\nabla u(\textbf{r}, t) = g
        \quad\text{for}\quad\textbf{r}\quad\text{on}\quad\delta\Omega=\Gamma_{\text{Neumann}}

    Parameters
    ----------
    rhs: :py:mod:`Vector`
        Right hand side vector of length node count.
    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, g ].
        The value :math:`g` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`g` need to be a scalar value (float or int) or
        a value generator function that will be executed at run time.

        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,
    nDim: int [1]
        Number of dimensions for vector valued problems. The rhs array need to
        have the correct size, i.e., number of Nodes * mesh.dimension()
    time: float
        Will be forwarded to value generator.
    userData: class
        Will be forwarded to value generator.
    dofOffset: int[0]
        Offset for matrix index.
    """
    if rhs is None:
        raise BaseException("Neumann Boundary condition needs rhs vector.")
    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    Se = pg.matrix.ElementMatrix()

    dof = len(rhs) // nDim

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        # print('+++++', boundary)
        # print('\t', val)
        g = generateBoundaryValue(boundary, val, time, userData, nCoeff=nCoeff)
        # print('\t', g)

        # if a is not None:
        #     pg.warning('Scaling of Neumann values necessary? Check!')
        #     try:
        #         g *= a[boundary.leftCell().id()]
        #     except BaseException as e:
        #         print(boundary.leftCell())
        #         print(boundary.leftCell().id())
        #         print(len(a))
        #         pg.warn('Insufficient cell information.')

        if g is not None:
            Se.u(boundary)
            for dim in range(nDim):
                if nDim == 1:
                    gd = g
                else:
                    if isinstance(g, list) and len(g) == nDim:
                        gd = g[dim]
                    else:
                        gd = g.T[dim]

                idx = Se.ids() + dim*dof + dofOffset

                if isinstance(gd, float) and gd == 0:
                    continue
                if hasattr(gd, '__iter__') and not np.any(gd):
                    continue

                # print(nDim, g, gd)
                if isinstance(rhs, pg.Vector):
                    # print(Se)
                    # pg.info(sum(Se.row(0)))
                    # pg.info(Se.row(0), gd, idx)
                    rhs.addVal(Se.row_RM(0) * gd, idx)
                    #rhs.setVal(Se.row(0) * gd, idx)
                    # rhs.add(Se, g)
                else:
                    # check
                    #pg.error('check')
                    rhs[idx] += Se.row_RM(0) * gd

                    # for i, j in enumerate(Se.ids()):
                    #     rhs[j + dim*dof] += Se.row(0)[i] * gd


def assembleRobinBC(mat, boundaryPairs, rhs=None, time=0.0, userData={},
                    dofOffset=0, nCoeff=1, dofPerCoeff=None):
    r"""Apply Robin boundary condition.

    Apply Robin boundary condition to the system matrix and the rhs vector:

    .. math::
        \frac{\partial u(\textbf{r}, t)}{\partial\textbf{n}}
        & = \alpha(u_0-u) \quad\text{or} \\
        \beta\frac{\partial u(\textbf{r}, t)}{\partial\textbf{n}}
        + \alpha u
        & = \gamma \\
        & \quad\text{for}\quad\textbf{r}\quad\text{on}\quad\delta\Omega=
        \Gamma_{\text{Robin}}\\

    Arguments
    ---------
    mat: :gimliapi:`GIMLI::SparseMatrix`
        System matrix of the system equation.

    boundaryPair: list
        List of pairs [:gimliapi:`GIMLI::Boundary`,
        :math:`a, u_0` | :math:`\alpha, \beta, \gamma`].
        The values will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        Values can be a single value for :math:`\alpha` or :math:`a`,
        two values will be interpreted as :math:`a, u_0`,
        and three values will be :math:`\alpha, \beta, \gamma`.
        Also generator (callable) is possible which will be executed at runtime
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        :ref:`tut:modelling_bc` or testing/test_FEM.py for example syntax.

    time: float
        Will be forwarded to value generator.
    userData: dict
        Will be forwarded to value generator.
    dofOffset: int[0]
        Offset for matrix index.
    """
    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    S_Dir = pg.matrix.ElementMatrix()
    S_Neu = pg.matrix.ElementMatrix()

    # if isinstance(rhs, np.ndarray):
    #     rhs = pg.Vector(rhs)

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        # print('val:', val)
        # du/dn = a(u0-u) || \beta du/dn + \alpha u = \gamma
        # combines to Matrix + au = RHS + au0

        u0 = None
        a = generateBoundaryValue(boundary, val, time, userData,
                                  expectList=True, nCoeff=nCoeff)

        try:
            if a.ndim == 2 and len(a) == boundary.nodeCount():
                a = a[0]
        except BaseException:
            # expecting [[a| a, u0 | a b g]_i] for i in boundary.nodes()
            print(boundary)
            print(a)
            print(a.ndim)
            pg.error("Can't interprete robin value.")

        if hasattr(a, '__iter__'):
            if len(a) == 1:
                a = a[0]
            elif len(a) == 2:
                u0 = a[1]
                a = a[0]
            elif len(a) == 3:
                alpha, beta, gamma = a[0], a[1], a[2]
                # a = [alpha, beta, gamma]
                if alpha != 0:
                    u0 = gamma/alpha
                else:
                    pg.warn('Robin boundary condition parameter alpha is zero, '
                            'falling back to Neumann condition.')
                    u0 = 0.0
                if beta != 0:
                    a = alpha/beta
                else:
                    pg.warn('Robin boundary condition parameter beta is zero, '
                            'please consider using Dirichlet instead.')
                    a = 0.0

        if a is not None and a != 0.0:
            S_Dir.u2(boundary)
            mat.add(S_Dir, f=a)
            # Sp *= p
            # S += Sp
        if u0 is not None and u0 != 0.0:
            S_Neu.u(boundary)
            rhs.add(S_Neu, a*u0, neg=False)


def assembleBC(bc, mesh, mat, rhs, time=None, userData={}, dofOffset=0,
               nCoeff=1):
    r"""Assemble all types of boundary conditions.

    Forward to appropriate assemble functions.

    Parameters
    ----------
    mesh : pg.Mesh
        mesh
    mat : pg.SparseMapMatrix
        system matrix
    rhs : vector
        right-hand side
    time : float|iterable
        time vector
    userData : dict
        user data passed to the assembly functions
    dofOffset : int
        offset for degree of freedom numbering

    Returns
    -------
    map{id: uDirichlet}: Map of index to Dirichlet value.
    """
    # we can't iterate because we want the following fixed order
    dirichletMap = {}
    bct = dict(bc)
    nDim = 1
    if mat is not None and mat.rows() == mesh.nodeCount() * mesh.dim():
        nDim = mesh.dim()

    if 'Neumann' in bct:
        assembleNeumannBC(rhs, parseArgToBoundaries(bct.pop('Neumann'), mesh),
                          nDim=nDim, time=time, userData=userData,
                          dofOffset=dofOffset,
                          nCoeff=nCoeff, dofPerCoeff=mesh.nodeCount())
    if 'Robin' in bct:
        assembleRobinBC(mat, parseArgToBoundaries(bct.pop('Robin'), mesh),
                        rhs=rhs, time=time, userData=userData,
                        dofOffset=dofOffset,
                        nCoeff=nCoeff, dofPerCoeff=mesh.nodeCount())
    if 'Dirichlet' in bct:
        uD = assembleDirichletBC(
            mat, parseArgToBoundaries(bct.pop('Dirichlet'), mesh),
            rhs=rhs, time=time, userData=userData,
            dofOffset=dofOffset,
            nCoeff=nCoeff, dofPerCoeff=mesh.nodeCount())
        dirichletMap.update(uD)

    if 'Nodes' in bct:
        # 'Nodes' : [list(Nodes), callable(Node)] ## for selected Nodes
        # 'Nodes' : callable(Node) ## for all nodes
        bc = bct.pop('Nodes')
        if isinstance(bc, list):
            nodes = bc[0]
            val = bc[1]
        else:
            nodes = mesh.nodes()
            val = bc

        nP = []
        if callable(val):
            for n in nodes:
                nP.append([n.id(), val(n)])
        else:
            pg.critical("Nodes boundary need a callable(Node)")

        uD = assembleDirichletBC(
            mat, [], nodePairs=nP,
            rhs=rhs, time=time, userData=userData, dofOffset=dofOffset,
            nCoeff=nCoeff, dofPerCoeff=mesh.nodeCount())
        dirichletMap.update(uD)

    if 'Node' in bct:
        uD = assembleDirichletBC(
            mat, [], nodePairs=bct.pop('Node'),
            rhs=rhs, time=time, userData=userData, dofOffset=dofOffset,
            nCoeff=nCoeff, dofPerCoeff=mesh.nodeCount())
        dirichletMap.update(uD)

    if len(bct.keys()) > 0:
        pg.warn("Unknown boundary condition[s]" +
                str(bct.keys()) + " will be ignored")

    return dirichletMap


def assembleLoadVector(mesh, f, userData={}):
    r"""Assemble the load vector. See createLoadVector.
    Maybe we will remove this
    """
    pg.deprecate('createLoadVector') # 20200115
    return createLoadVector(mesh, f, userData)


def createForceVector(mesh, f, userData={}):
    """Create a right hand side vector for vector valued solutions.

    Parameters
    ----------
    f: convertible
        List of rhs side options. Must be convertable to createLoadVector.
        See :py:mod:`createLoadVector`
    rhs: np.array()
        Squeezed vector of length mesh.nodeCount() * mesh.dimensions()
    """
    if not isinstance(f, list):
        pg.error("Create Force Vector need list of attribute f with an entry "
                 "for each dimension.")

    rhs = np.zeros(mesh.nodeCount() * mesh.dim())

    for i in range(mesh.dim()):
        rhs[i*mesh.nodeCount():(i+1)*mesh.nodeCount()] = \
            createLoadVector(mesh, f[i], userData)

    # rhs.reshape(mesh.nodeCount() * mesh.dim()) #contiguity not guarantied
    return rhs


def createLoadVector(mesh, f=1.0, userData={}):
    """Create right-hand-side vector from load/force vectors.

    Based on given mesh and load values (scalar solution) or
    force vectors (vector value solution).

    TODO
    ----
        * Callable for vector valued problems
        * Callable called dynamic on demand

    Parameters
    ----------
    f: float[1.0], array, callable(cell, [userData]), [f_x, f_y, f_z]

        * float will be assumed as constant for all cells
            like rhs = rhs(np.ones(mesh.cellCount() * f),
        * array of length mesh.cellCount() will be processed as load value for
            each cell: rhs = rhs(f),
        * array of length mesh.nodeCount() is assumed to be already processed
            correct: rhs = f
        * callable is evaluated on once for each cell and need to return a load
            value for each cell and can have  optional a userData dictionary:
            `f_cell = f(cell, [userData={}])`
            rhs = rhs(f(c, userData) for c in mesh.cells())
        * list with length of mesh.dimension() of float or array entries will
            create a squeezed rhs for vector valued problems
            rhs = squeeze([rhs(f[0]), rhs(f[1]), rhs(f[2])])

    Returns
    -------
    rhs: pg.Vector(mesh.nodeCount())
        Right-hand side load vector for scalar values or squeezed vector values
    """
    # f is dict('Node':callable, 'Cell': callable)
    if isinstance(f, dict):
        if 'Node' in f:
            fn = []
            if callable(f['Node']):
                for n in mesh.nodes():
                    fn.append(f['Node'](n, **userData))
            if hasattr(fn[0], '__iter__'):
                # result is vector valued
                return createLoadVector(mesh, list(np.array(fn).T),
                                        userData=userData)

            return createLoadVector(mesh, fn, userData=userData)

        elif 'Cell' in f:
            pg.error('Implement me!, createLoadVector()')

    # fix for the lazy
    if isinstance(f, int):
        f = float(f)

    # f is list [fx, fy, [fz]] for vector problems
    if isinstance(f, list) and len(f) == mesh.dim():
        return createForceVector(mesh, f, userData=userData)

    # f is list of array [f_0, f_1, ..., f_n] for scalar problems
    if isinstance(f, list) or hasattr(f, 'ndim'):
        if isinstance(f, list):
            rhs = np.zeros((len(f), mesh.nodeCount()))
            for i, fi in enumerate(f):
                userData['i'] = i
                rhs[i] = createLoadVector(mesh, fi, userData)
            return rhs

        elif f.ndim == 2:
            # assume rhs [n, nNodes] array is already a valid
            if len(f[0]) == mesh.nodeCount():
                return f

    rhs = pg.Vector(mesh.nodeCount(), 0)

    fArray = None

    if hasattr(f, '__len__'):
        if len(f) == mesh.cellCount():
            # scalar values for each cell
            fArray = f
        elif len(f) == mesh.nodeCount():
            # scalar values for each node
            fArray = f
        elif len(f) == mesh.nodeCount() * mesh.dim():
            # vector values for each node
            # maybe just for special cases with allready processed rhs
            return f

    elif callable(f) and not isinstance(f, pg.Vector):
        fArray = pg.Vector(mesh.cellCount())
        for c in mesh.cells():
            if userData is not None and userData.keys():
                fArray[c.id()] = f(c, userData)
            else:
                fArray[c.id()] = f(c)

    if fArray is None:
        fArray = cellValues(mesh, f, userData=userData)

    if len(fArray) == mesh.cellCount():
        b_l = pg.matrix.ElementMatrix()

        for c in mesh.cells():
            if fArray[c.id()] != 0.0:
                b_l.u(c)
                rhs.add(b_l, fArray[c.id()])

#            print("test reference solution:")
#            rhsRef = pg.Vector(mesh.nodeCount(), 0)
#            for c in mesh.cells():
#                b_l.u(c)
#                for i, idx in enumerate(b_l.idx()):
#                    rhsRef[idx] += b_l.row(0)[i] * fArray[c.id()]
#            np.testing.assert_allclose(rhs, rhsRef)
#            print("Remove revtest in assembleLoadVector after check")

    elif len(fArray) == mesh.nodeCount():
        # nodal values
        fA = pg.Vector(fArray)
        b_l = pg.matrix.ElementMatrix()
        for c in mesh.cells():
            b_l.u(c)
            # rhs.addVal(b_l.row(0) * fArray[b_l.idx()], b_l.idx())
            rhs.add(b_l, fA)
        # print("test reference solution:")
        # rhsRef = pg.Vector(mesh.nodeCount(), 0)
        # for c in mesh.cells():
        #     b_l.u(c)
        #     for i, idx in enumerate(b_l.idx()):
        #         rhsRef[idx] += b_l.row(0)[i] * fA[idx]

        # np.testing.assert_allclose(rhs, rhsRef)
        # print("Remove revtest in assembleLoadVector after check",
        #       sum(rhs), sum(rhsRef))

            # rhs = pg.Vector(fArray)

    else:
        raise Exception("Load vector have the wrong size: " +
                        str(len(fArray)))

    return rhs


def createStiffnessMatrix(mesh, a=None, isVector=False):
    r"""Create the Stiffness matrix.

    Calculates the Stiffness matrix :math:`{\bf S}` for the given mesh scaled
    with the per cell values a.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Arbitrary mesh to calculate the stiffness for.
        Type of base and shape functions depends on the cell types.

    a : iterable of type float, int, complex, :gimliapi:`GIMLI::Matrix`
        Per cell values., e.g., physical parameter.
        Length of a need to be `mesh.cellCount()`.
        If None given default is 1.

    isVector : bool [False]
        We want to solve for vector valued problems. Resulting SparseMatrix is
        a SparseMapMatrix and have the dimension
        `(nNodes x nDims, nNodes x nDims)` with `nNodes = mesh.nodeCount()` and
        `nDims = mesh.dimension()`.

    Returns
    -------
    A : :gimliapi:`GIMLI::SparseMatrix` | SparseMapMatrix
        Stiffness matrix, with real or complex values.
    """
    if mesh.cellCount() == 0:
        print(mesh)
        raise Exception("Mesh invalid")

    if isinstance(a, complex):
        a = np.full(mesh.cellCount(), a)

    if a is None:
        a = pg.Vector(mesh.cellCount(), 1.0)

    if isinstance(a, (float, int)):
        a = pg.Vector(mesh.cellCount(), a)

    A = None


    if isVector is False:
        if isinstance(a[0], (float, int, np.float64)):
            A = pg.matrix.SparseMatrix()
            A.fillStiffnessMatrix(mesh, a)
            return A

        dof = 0
        nDof = mesh.nodeCount()
    else:
        dof = mesh.nodeCount()
        nDof = mesh.nodeCount() * mesh.dimension()

    # if vector or scalar(Complex)
    if pg.isComplex(a[0]):
        isComplex = True
        A = pg.matrix.CSparseMapMatrix(nDof, nDof)
    else:
        isComplex = False
        A = pg.matrix.SparseMapMatrix(nDof, nDof)

    al = pg.core.ElementMatrix(dof=dof)

    if len(a) != mesh.cellCount():
        pg.error('Number of cell values need to match cell count')

    for c in mesh.cells():
        if isComplex is True:
            # al.gradU2(c, 1.0)
            al.ux2uy2uz2(c)
            A.add(al, f=a[c.id()], scale=1.0)
        else:
            if pg.isScalar(a[c.id()]):
                al.gradU2(c, a[c.id()])
                A.add(al)
            else:
                if hasattr(a[c.id()], 'voigtNotation'):
                    vN = a[c.id()].voigtNotation
                else:
                    vN = False

                # al.gradU2(c, a[c.id()], voigtNotation=vN)
                al.gradU2(c, np.array(a[c.id()]), voigtNotation=vN)
                A.add(al)

    # if isComplex is True:
    #     return pg.matrix.CSparseMatrix(A)

    return pg.matrix.asSparseMatrix(A)


def createMassMatrix(mesh, b=None):
    r"""Create the mass matrix.

    Calculates the Mass matrix (Finite element identity) for the given mesh.

    ..math::
            ...

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`

        Arbitrary mesh to calculate the mass element matrix.
        Type of base and shape functions depends on the cell types.

    b : array
        Per cell values. If None given default is 1.

    Returns
    -------
    A : :gimliapi:`GIMLI::SparseMatrix`
        Mass element matrix
    """
    # need callable here
    if b is None:
        b = pg.Vector(mesh.cellCount(), 1.0)
    elif not hasattr(b, '__iter__'):
        b = pg.Vector(mesh.cellCount(), b)

    B = pg.matrix.SparseMatrix()
    B.fillMassMatrix(mesh, b)
    return B

    # create matrix structure regarding the mesh
    # B.buildSparsityPattern(mesh)
    # define a local element matrix
    # B_l = pg.matrix.ElementMatrix()
    # for c in mesh.cells():
    #    B_l.u2(c)
    # # check if b[i] == B*b
    #    B_l *= b[c.id()]
    #    B += B_l
    # return B


def intDomain(u, mesh=None):
    r"""Return integral over nodal solution :math:`u`.

    .. math::
        \int_{\Omega} u

    TODO
    ----
        * refactor
        * better name?
        * Documentation
    """
    if mesh is not None:
        r = createLoadVector(mesh)
        return sum(r*u)
    pg.critical('Need a mesh to calculate the integral over domain')


def _feNorm(u, mat):
    """Create a norm within a Finite Element space.

    Create the Finite Element Norm with a preassembled system matrix.
    """
    if u.ndim == 2:
        u1 = np.array(u).T.reshape(u.shape[0]*u.shape[1])
        # same like:
        # u1 = pg.cat(u[:,0], u[:,1])
        return np.sqrt(pg.math.dot(u1, mat.mult(u1)))
    return np.sqrt(pg.math.dot(u, mat.mult(u)))


def normL2(u, mat=None, mesh=None, warn=True):
    r"""Create Lebesgue (L2) norm for finite element space.

    Note
    ----
        * wrong implementation .. better to come

    Find the L2 norm for a solution in the finite element space. :math:`u` exact solution,
    :math:`{\bf M}` Mass matrix, i.e., Finite element identity matrix.

    .. math::

        L2(f(x)) = || f(x) ||_{L^2} & = (\int |f(x)|^2 \mathrm{d}\:x)^{1/2} \\
                                    & \approx h (\sum |f(x)|^2 )^{1/2} \\
        L2(u) = || u ||_{L^2} & = (\int |u|^2 \mathrm{d}\:x)^{1/2} \\
                              & \approx (\sum M (u))^{1/2} \\
        e_{L2_rel} = \frac{L2(u)}{L2(u)} & =
                               \frac{(\sum M(u))^{1/2}}{(\sum M u)^{1/2}}

    The error for any approximated solution :math:`u_h` correlates to the L2
    norm of 'normL2(u - u_h, M)'. If you like relative values, you can also
    normalize this error with 'normL2(u - u_h, M) / normL2(u, M)*100'.

    Parameters
    ----------
    u: iterable
        Node based value to compute the L2 norm for.

    mat: Matrix
        Mass element matrix.

    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh with the FE space to generate M if necessary.

    warn: bool [True]
        Warn if you forget a matrix or mass matrix which leads to algebraic L2

    Returns
    -------
    ret: float
        :math:`L2(u)` norm.

    """
    if isinstance(mat, pg.Mesh):
        mesh = mat
        mat = None

    if mat is None and mesh is not None:
        mat = createMassMatrix(mesh)

    if mat is None:
        if warn:
            pg.warning("No mass matrix or a mesh here, to calculate L2-Norm. "
                   "Returning algebraic l2.")

        # M is Identity matrix
        return np.sqrt(np.dot(u, u))

    return _feNorm(u, mat)


def normH1(u, mat=None, mesh=None):
    r"""Create (H1) norm for the finite element space.

    Note
    ----
        * wrong implementation .. better to come

    Parameters
    ----------
    u : iterable
        Node based value to compute the H1 norm for.

    mat : Matrix
        Stiffness matrix.

    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh with the FE space to generate S if necessary.

    Returns
    -------
    ret : float
        :math:`H1(u)` norm.

    """
    if isinstance(mat, pg.Mesh):
        mesh = mat
        mat = None

    if mat is None and mesh is not None:
        mat = pg.solver.createStiffnessMatrix(mesh)

    if mat is None:
        raise Exception("Need Stiffness matrix or mesh to calculate H1 norm")

    return _feNorm(u, mat)


def solve(mesh, **kwargs):
    r"""Solve partial differential equation.

    This is a syntactic sugar convenience function for solving partial
    differential equation on a given mesh.
    Using the best guess method for the given parameter.
    Currently only Finite Element calculation using
    :py:mod:`pygimli.solver.solveFiniteElements`

    TODO
    :py:mod:`pygimli.solver.solveFiniteVolume`

    """
    return solveFiniteElements(mesh, **kwargs)


def solveFiniteElements(mesh, a=1.0, b=None, f=0.0, bc=None,
                        times=None, c=1.0, userData={},
                        verbose=False, **kwargs):
    r"""Solve partial differential equation with Finite Elements.

    This is a syntactic sugar convenience function for using the Finite Element
    functionality of the library core to solve partial differential equation
    (PDE) that match the following form:

    .. math::

        c \frac{\partial u}{\partial t} & = \nabla\cdot(a \nabla u)
        + b u + f(\mathbf{r},t)~~|~~\Omega_{\text{Mesh}}\\
        u & = h~~|~~\Gamma_{\text{Dirichlet}}\\
        a\frac{\partial u}{\partial \mathbf{n}} & =
        g~~|~~\Gamma_{\text{Neumann}}\\
        \alpha u + \beta\frac{\partial u}{\partial \mathbf{n}} & =
        \gamma~~|~~\Gamma_{\text{Robin}}\\
        \frac{\partial u}{\partial \mathbf{n}} & =
        \alpha(u_0-u)~~|~~\Gamma_{\text{Robin}}

    for the scalar :math:`u(\mathbf{r}, t)` or vector
    :math:`\mathbf(u)(\mathbf{r}, t)` solution at each node of a given mesh.
    The Domain :math:`\Omega` and the Boundary :math:`\Gamma` are defined
    through the mesh with appropriate boundary marker.

    To ensure vector solution, either set vector forces or at least one
    vector component boundary condition.

    TODO
    ----
    * unsteady ub and dub
    * 'Infinity' Boundary condition (u vanishes at infinity)

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    a: value | array | callable(cell, userData)
        Cell values of type float or complex can be scalar, anisotropy matrix
        or elastic tensor.
    b: value | array | callable(cell, userData) [None]
        Cell values. None means the term is unused.
    c: value | array | callable(cell, userData) [None]
        Scale the unsteady term, only for times is not None.
    f: value | array(cells) | array(nodes) | callable(args, kwargs)
        force values, for vector fields use (n x dim) values.
    bc: dict()
        Dictionary of boundary conditions.
        Current supported boundary conditions by dictionary keys:
        ``Dirichlet``, ``Neumann``, ``Robin``, ``Node``.

        The dictionary can contain multiple `key: Arg`. `Arg` will be parsed
        by :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray`

        If the dictionary key is ``Node`` then fixed values for single node
        indices can by be given. e.g., `bc={'Node': [nodeID, value]}`.
        Note this is only a shortcut
        for `bc={'Dirichlet': [mesh.node(nodeID), value]}`.

        The parameter `a` for Neumann boundary condition is choosen
        automatically from the diffusivity parameter `a` of the associated cell.

    times: array [None]
        Solve as time dependent problem for the given times.

    Keyword Arguments
    -----------------
    **kwargs
        u0: value | array | callable(pos, userData)
            Node values
        theta: float [1]
            * :math:`\theta = 0` means explicit Euler, maybe stable for
                :math:`\Delta t \quad\text{near}\quad h`
            * :math:`\theta = 0.5`, Crank-Nicolson scheme, maybe instable
            * :math:`\theta = 2/3`, Galerkin scheme
            * :math:`\theta = 1`, implicit Euler

            If unsure choose :math:`\theta = 0.5 + \epsilon` (probably stable).

        dynamic: bool [False]
            Boundary conditions for time depending problems will be considered
            dynamic for each time step.
        stats: bool
            Give some statistics.
        progress: bool
            Give some calculation progress.
        assembleOnly: bool
            Stops after matrix asssemblation.
            Returns the system matrix A and the rhs vector.
        fixPureNeumann: bool [auto]
            If set or detected automatic, we add the additional condition:
            :math:`\int_\Omega u dv = 0` making elliptic problems well-posed.
        rhs: iterable
            Pre assembled rhs. Will preferred on any f settings.
        ws: dict
            The WorkSpace is a dictionary that will get
            some temporary data during the calculation.
            Any keyvalue 'u' in the dictionary is used for the resulting array.
        vectorValued: bool (False)
            Solution forced to vector valued, in case the auto detection fails

    Returns
    -------
    u: array
        Returns the solution u either 1,n array for stationary problems or
        for m,n array for m time steps

    See also
    --------
    :ref:`tut:modelling` and :py:mod:`pygimli.solver.solve`

    Examples
    --------
    >>> # no need to import matplotlib, pygimli show does.
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> world = mt.createWorld(start=[-10, 0], end=[10, -10],
    ...                        marker=1, worldMarker=False)
    >>> c1 = mt.createCircle(pos=[0.0, -5.0], radius=3.0, area=.1, marker=2)
    >>> mesh = mt.createMesh([world, c1], quality=34.3)
    >>> u = pg.solver.solveFiniteElements(mesh, a={1: 100.0, 2: 1.0},
    ...                                   bc={'Dirichlet':{4: 1.0, 3: 0.0}})
    >>> ax = pg.show(mesh, u, showMesh=True)[0]
    >>> _ = pg.show(c1, ax=ax, fillRegion=False)
    """
    if bc is None:
        bc = {}

    workSpace = kwargs.pop('ws', dict())
    debug = kwargs.pop('debug', False)
    stats = kwargs.pop('stats', False)

    mesh.createNeighborInfos()
    if verbose:
        print("Mesh: ", str(mesh))

    # scalar solution default
    vectorValues = False
    dof = mesh.nodeCount()

    # check if force vector is a vector
    rhs = kwargs.pop('rhs', createLoadVector(mesh, f, userData=userData))

    # pg._g('###############')
    if len(rhs) > dof or kwargs.pop('vectorValued',
                                    _bcIsForVectorValues(bc, mesh)):
        if verbose:
            print("Solve vector valued.")
        vectorValues = True
        dof = mesh.nodeCount() * mesh.dimension()

    # pg._g('###############', dof)
    rhs.resize(dof)

    swatch = pg.core.Stopwatch(True)

    # check for material parameter
#   a = parseArgToArray(a, nDof=mesh.cellCount(), mesh=mesh, userData=userData)
    a = cellValues(mesh, a, userData=userData)
    isComplex = False
    if pg.utils.isComplex(a):
        isComplex = True
        rhs = np.array(rhs, dtype=complex)

    S = createStiffnessMatrix(mesh, a, isVector=vectorValues)
    M = None

    if b is not None and b != 0:
        b = cellValues(mesh, b, userData=userData)
        M = createMassMatrix(mesh, b)
        # pg.warn("check me")
        A = S - M
    else:
        A = S

    if times is None:
        if len(list(bc.items())) == 0 or \
           (len(list(bc.items())) == 1 and list(bc.keys())[0] == 'Neumann'):
            pn = True
        else:
            pn = False

        fixPureNeumann = kwargs.pop('fixPureNeumann', pn)

        assembleBC(bc, mesh, A, rhs, time=None, userData=userData)

        u = None
        if 'u' in workSpace:
            u = workSpace['u']

        singleForce = True

        if hasattr(rhs, 'ndim'):
            if rhs.ndim == 2:
                singleForce = False
                if u is None:
                    u = np.zeros(rhs.shape)
        else:
            if isinstance(a[0], complex):
                if u is None:
                    u = pg.CVector(rhs.size(), 0.0)
                rhs = pg.math.toComplex(rhs)
            else:
                if u is None:
                    u = pg.Vector(rhs.size(), 0.0)

        if fixPureNeumann is True:
            pg.info('Fixing pure Neumann boundary condition by forcing: '
                    'intDomain(u, mesh) = 0')
            r = createLoadVector(mesh)

            A = pg.BlockMatrix()
            A.add(S, 0, 0)
            A.add(r, 0, mesh.nodeCount())
            A.add(r, mesh.nodeCount(), 0, transpose=True)

            rhs = pg.cat(rhs, pg.Vector(1, 0))

        assembleTime = swatch.duration(True)

        if stats:
            stats.assembleTime = assembleTime

        if verbose:
            print("Assembling time: ", assembleTime)

        workSpace['Stiffness matrix'] = S
        workSpace['Mass matrix'] = M
        workSpace['System matrix'] = A
        workSpace['rhs'] = rhs

        if 'assembleOnly' in kwargs:
            return A, rhs

        if fixPureNeumann:
            if singleForce:
                uc = pg.solver.linSolve(A, rhs, 'scipy')
                u = uc[0:mesh.nodeCount()]
            else:
                pg.critical(
                    'Non-single force for pure Neumann not yet implemented')
        else:
            solver = pg.core.LinSolver(False)
            solver.setMatrix(A, 0)

            if singleForce:
                if isComplex is True:
                    # clean this up
                    rhs = pg.core.toComplex(rhs.real, rhs.imag)
                    u = solver.solve(rhs).array()
                else:
                    u = solver.solve(rhs)
            else:
                for i, r in enumerate(rhs):
                    u[i] = solver.solve(r)

        solverTime = swatch.duration(True)
        if verbose:
            if stats:
                stats.solverTime = solverTime
            print("Solving time: ", solverTime)

        if len(kwargs.keys()) > 0:
            pg.warn("Unused arguments", *kwargs)
        return u

    else:  # times given
        pg.solver.checkCFL(times, mesh, max(np.array(a).flatten()))

        if debug:
            print("start TL", swatch.duration())

        if c != 1.0:
            c = cellValues(mesh, c, userData=userData)

        M = createMassMatrix(mesh, c)
        F = createLoadVector(mesh, f)

        u0 = np.zeros(dof)
        if 'u0' in kwargs:
            u0 = parseArgToArray(kwargs['u0'], dof, mesh, userData)

        progress = None
        if 'progress' in kwargs:
            from pygimli.utils import ProgressBar
            progress = ProgressBar(its=len(times), width=40, sign='+')

        theta = kwargs.pop('theta', 1.0)
        dynamic = kwargs.pop('dynamic', False)

        if not dynamic:
            S = createStiffnessMatrix(mesh, a)
            assembleBC(bc, mesh, S, F, time=0.0, userData=userData)
            return crankNicolson(times, S, M, f=F,
                                 u0=u0, theta=theta,
                                 progress=progress)

        rhs = np.zeros((len(times), dof))
        # no time dependency for rhs so far ... TODO
        rhs[:] = F  # this is slow: optimize

        if debug:
            print("rhs", swatch.duration())
        U = np.zeros((len(times), dof))
        U[0, :] = u0

        # init state
        u = pg.Vector(dof, 0.0)

        if debug:
            print("u0", swatch.duration())

        measure = 0.
        for n in range(1, len(times)):
            swatch.reset()

            dt = times[n] - times[n-1]

            # previous timestep
            # print "i: ", i, dt, U[i - 1]

            swatch.reset()
            # (A + a*B)u is fastest,
            # followed by A*u + (B*u)*a and finally A*u + a*B*u and
            br = (M + S*(dt * (theta - 1.))) * U[n - 1] + \
                dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

            # print ('a',swatch.duration(True))
            # br = M * U[n - 1] - (A * U[n - 1]) * (dt*(1.0 - theta)) + \
            # dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

            # print ('br',swatch.duration(True))

            # br = M * U[n - 1] - (dt*(1.0 - theta)) * A * U[n - 1] + \
            # dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
            # print ('c',swatch.duration(True))

            measure += swatch.duration()

            A = M + S * dt * theta

            assembleBC(bc, mesh, A, br, time=times[n], userData=userData)

            if 'assembleOnly' in kwargs:
                return A, br

            # u = S/b
            t_prep = swatch.duration(True)
            solver = pg.core.LinSolver(A, verbose)
            solver.solve(br, u)

            if 'plotTimeStep' in kwargs:
                kwargs['plotTimeStep'](u, times[n])

            U[n, :] = np.asarray(u)

            if progress:
                progress.update(n, 't_prep: {0}ms t_step {1}s'.format(
                    pg.pf(t_prep*1000),
                    pg.pf(swatch.duration())))
        if debug:
            print("Measure(" + str(len(times)) + "): ",
                  measure, measure / len(times))
        return U


def checkCFL(times, mesh, vMax, verbose=False):
    r"""Check Courant-Friedrichs-Lewy condition.

    For advection and flow problems. CFL Number should be lower then 1 to
    ensure stability.

    .. math ::

        cMax = \Delta t \sum_n \left(\frac{vMax_n}{\Delta x_n}\right)

    with $n$ be number of spatial dimensions and $vMax_n$ and $x_n$ the dimensional component. For unstructured meshes we apply:

    .. math ::

        cMax = \Delta t \left(\frac{vMax}{min(mesh_h)}\right) * dim

    Parameters
    ----------
    """
    if pg.isScalar(times):
        dt = times
    else:
        dt = times[1] - times[0]

    dx = min(mesh.h())
    # min(entity.shape().h()
    # if mesh.dimension() == 1:
    #     dx = min(mesh.cellSizes())
    # else:
    #     dx = min(mesh.boundarySizes())
    c = (vMax * dt / dx) * mesh.dim()

    if c > 1:
        pg.warn("Courant-Friedrichs-Lewy Number:", c,
                "but should be lower 1 to ensure movement inside a cell "
                "per timestep. ("
                "vMax =", vMax,
                "dt =", dt,
                "dx =", dx,
                "dt <", dx/vMax,
                " | N > ", int(dt/(dx/vMax))+1, ")")
    if verbose:
        pg.info("Courant-Friedrichs-Lewy Number:", c)
    return c


def crankNicolson(times, S, I, f=None,
                  u0=None, theta=1.0, dirichlet=None,
                  solver=None, progress=None):
    """Solve time dependend problems using Crank-Nicolson scheme.

    Limitations so far:
        S = Needs to be constant over time (i.e. no change in coefficients)
        f = constant over time (would need assembling in every step)

    Args
    ----
    times: iterable(float)
        Timeteps to solve for. Give at least 2. Calculates len(times)-1 time
        steps. If times is not equidistant the matrix factorization is done
        for every time step.
    S: Matrix
        System matrix holds your discrete equations and boundary conditions
    I: Matrix
        Identity matrix (FD, FV) or Masselementmatrix (FE) to handle solution
        vector
    u0: iterable [None]
        Starting condition. zero or None if not given
    f: iterable (float) [None]
        External forces. Note f might also contain compensation values due to
        algebraic Dirichlet correction of S
    theta: float [1.0]
        * 0: Backward difference scheme (implicit)
        * 1: Forward difference scheme (explicit)
            strong time steps dependency .. will be unstable for to small values
        * 0.5: probably best tradeoff but can also be unstable

    dirichlet: dirichlet generator
        Genertor object to applay dirichlet boundary conditions
    solver: LinSolver [None]
        Provide a pre configured solver if you want some special.
    progress: Progress [None]
        Provide progress object if you want to see some.

    Returns
    -------
    np.ndarray:
        Solution for each time steps. First result step is input u0.
    """
    if len(times) < 2:
        raise BaseException("We need at least 2 times for "
                            "Crank-Nicolsen time discretization." + str(len(times)))

    swatches = pg.tictoc
    with swatches('Crank-Nicolson'):

        with swatches('prep'):

            if progress:
                timeMeasure = True

            dof = S.rows()

            rhs = np.zeros((len(times), dof))
            if f is not None:
                rhs[:] = f

            u = np.zeros((len(times), dof))

            if u0 is not None and not pg.isScalar(u0, 0):
                u[0, :] = u0

            if theta == 0:
                A = I.copy()

            if isinstance(solver, str):
                solver = pg.solver.LinSolver(solver=solver)
            if solver is None:
                solver = pg.solver.LinSolver(solver='scipy')


        dt = 0.0
        for n in range(1, len(times)):

            newDt = times[n] - times[n-1]

            if newDt < 1e-8 :
                pg.critical('Cannot find delta t for times', times)

            if abs(newDt - dt) > 1e-8:
                with swatches('factorize'):

                    ## new dt, so we need to factorize the matrix again
                    dt = newDt
                    # pg.info('dt', dt)

                    A = I + S * (dt * theta)

                    if dirichlet is not None:
                        dirichlet.apply(A, time=times[n])

                    solver.factorize(A)

                    St = None

            with swatches('build'):

                if theta == 0:
                    if St is None:
                        St = I - S * dt  # cache what's possible
                    b = St * u[n-1] + dt * rhs[n-1]
                elif theta == 1:
                    b = I * u[n-1] + dt * rhs[n]
                else:
                    if St is None:
                        St = I - S * (dt*(1.-theta))  # cache what's possible
                    b = St * u[n-1] + dt * ((1.0 - theta) * rhs[n-1] + theta * rhs[n])

            with swatches('dirichlet'):
                if dirichlet is not None:
                    dirichlet.apply(b, time=times[n])

            with swatches('solve'):
                u[n, :] = solver(b)

            if progress:
                progress.update(n)
                    # n, 't_prep: ' + pg.pf(timeAssemble[-1]*1000) + 'ms ' +
                    # 't_step: ' + pg.pf(timeSolve[-1]*1000) + 'ms')

            # if verbose and (n % verbose == 0):
            #     # print(min(u[n]), max(u[n]))
            #     print("timesteps:", n, "/", len(times),
            #           'runtime:', sw.duration(), "s",
            #           'assemble:', np.mean(timeAssemble),
            #           'solve:', np.mean(timeSolve))
    return u


class RungeKutta(object):
    """TODO DOCUMENT ME"""
    rk4a = [0.0,
            -567301805773.0/1357537059087.0,
            -2404267990393.0/2016746695238.0,
            -3550918686646.0/2091501179385.0,
            -1275806237668.0/842570457699.0]
    rk4b = [1432997174477.0/9575080441755.0,
            5161836677717.0/13612068292357.0,
            1720146321549.0/2090206949498.0,
            3134564353537.0/4481467310338.0,
            2277821191437.0/14882151754819.0]
    rk4c = [0.0,
            1432997174477.0/9575080441755.0,
            2526269341429.0/6820363962896.0,
            2006345519317.0/3224310063776.0,
            2802321613138.0/2924317926251.0]

    def __init__(self, solver, verbose=False):
        """TODO DOCUMENT_ME"""
        self.solver = solver
        self.verbose = verbose
        self.order = 5
        self.dt = None
        self.time = 0
        self.tMax = None
        self.u = None
        self.resU = None
        self.nSteps = 0

    def run(self, u0, dt, tMax=1):
        """TODO DOCUMENT_ME"""
        self.start(u0, dt, tMax)

        for _ in range(self.nSteps):
            self.step()

        return self.u

    def start(self, u0, dt, tMax=1):
        """TODO DOCUMENT_ME"""
        self.nSteps = int(np.ceil(tMax/dt))
        self.dt = dt
        self.time = 0
        self.tMax = tMax
        self.u = deepcopy(u0)
        self.resU = deepcopy(u0)

        if isinstance(self.resU, list):
            for r in self.resU:
                r *= 0.0
        else:
            self.resU *= 0.0

    def step(self):
        """TODO DOCUMENT ME"""
        if self.time + self.dt > self.tMax:
            self.dt = self.tMax - self.time

        if self.order == 1:
            # explicit Euler
            k1 = self.solver.explicitRHS(self.u,
                                         self.time)
            self.u += self.dt * k1

        elif self.order == 3:
            k1 = self.solver.explicitRHS(self.u, self.time)
            k1 = self.u + self.dt * k1

            k2 = self.solver.explicitRHS(k1, self.time)
            k2 = (3*self.u + k1 + self.dt*k2)/4

            k3 = self.solver.explicitRHS(k2, self.time)

            self.u = (self.u + 2*k2 + 2*self.dt*k3)/3

        elif self.order == 4:
            # classical 4 step Runge-Kutta rk4
            k1 = self.solver.explicitRHS(self.u,
                                         self.time)
            k2 = self.solver.explicitRHS(self.u + self.dt/2 * k1,
                                         self.time + self.dt/2)
            k3 = self.solver.explicitRHS(self.u + self.dt/2 * k2,
                                         self.time + self.dt/2)
            k4 = self.solver.explicitRHS(self.u + self.dt * k3,
                                         self.time + self.dt)
            self.u += 1./6. * self.dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        elif self.order == 5:
            # low storage Version of rk4
            for jRK in range(5):
                tLocal = self.time + self.rk4c[jRK] * self.dt

                rhs = self.solver.explicitRHS(self.u, tLocal)

                if isinstance(self.resU, list):
                    for i in range(len(self.resU)):
                        self.resU[i] = self.rk4a[jRK] * self.resU[i] + \
                            self.dt * rhs[i]

                        self.u[i] += self.rk4b[jRK] * self.resU[i]
                else:
                    self.resU = self.rk4a[jRK] * self.resU + self.dt * rhs
                    self.u += self.rk4b[jRK] * self.resU

        self.time += self.dt
        return self.u


if __name__ == "__main__":
    pass
