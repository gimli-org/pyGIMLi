import pygimli as pg


def streamline(mesh, field, startCoord, dLengthSteps, dataMesh=None,
               maxSteps=1000, verbose=False, coords=(0, 1)):
    """Create a streamline from start coordinate and following a vector field in up and down direction.
    """
    xd, yd, vd = streamlineDir(mesh, field, startCoord, dLengthSteps,
                               dataMesh=dataMesh, maxSteps=maxSteps,
                               down=True,  verbose=verbose, coords=coords)

    c = mesh.findCell(startCoord)

    if c is not None:
        c.setValid(True)

    xu, yu, vu = streamlineDir(mesh, field, startCoord, dLengthSteps,
                               dataMesh=dataMesh, maxSteps=maxSteps,
                               down=False, verbose=verbose, coords=coords)

    return xd + xu[1:], yd + yu[1:], vd + vu[1:]


def streamlineDir(mesh, field, startCoord, dLengthSteps, dataMesh=None,
                  maxSteps=150, down=True, verbose=False, coords=(0, 1)):
    """
        down = -1, up = 1, both = 0
    """
    xd = []
    yd = []
    vd = []

    pot = None
    vx = None
    vy = None
    isVectorData = False

    if isinstance(field, pg.PosVector):
        field = field.array()

    if hasattr(field[0], '__len__'):
        if abs(max(field[:, 0])) == 0 and abs(max(field[:, 1]) == 0):
            raise Exception("No data range streamline: min/max == ",
                                min(field[:, 0]))

        vx = pg.Vector(field[:, 0])
        vy = pg.Vector(field[:, 1])

        isVectorData = True
    else:
        if min(field) == max(field):
            raise Exception("No scalar data range for any gradients "
                                " to draw a streamline: min/max == ",
                                min(field))

        if dataMesh is not None:
            if len(field) == dataMesh.nodeCount():
                pot = pg.Vector(field)
            elif len(field) == dataMesh.cellCount():
                pot = pg.core.cellDataToPointData(dataMesh, field)
            else:
                print(len(field), dataMesh)
                raise Exception("Data length (%i) for streamline is "
                    "neighter nodeCount (%i) nor cellCount (%i)" %
                    (len(field), dataMesh.nodeCount(), dataMesh.nodeCount()))
        else:
            if len(field) == mesh.nodeCount():
                pot = pg.Vector(field)
            elif len(field) == mesh.cellCount():
                pot = pg.core.cellDataToPointData(mesh, field)
            else:
                print(len(field), dataMesh)
                raise Exception("Data length (%i) for streamline is "
                    "neighter nodeCount (%i) nor cellCount (%i)" %
                    (len(field), mesh.nodeCount(), mesh.nodeCount()))

    direction = 1
    if down:
        direction = -1

    # search downward
    pos = pg.RVector3(startCoord)
    c = mesh.findCell(startCoord)
    dLength = c.center().dist(c.node(0).pos()) / dLengthSteps

    # stream line starting point
    if c is not None:
        xd.append(pos[coords[0]])
        yd.append(pos[coords[1]])
        vd.append(-1)

    lastC = c
    lastU = -direction * 1e99

    d = None
    while c is not None and len(xd) < maxSteps:

        # valid .. temporary check if there is already a stream within the cell
        if not c.valid():
            break

        if isVectorData:
            u = 0.
            if len(vx) == mesh.cellCount():
                d = pg.RVector3(vx[c.id()], vy[c.id()])
            elif len(vx) == mesh.nodeCount():
                d = pg.RVector3(c.pot(pos, vx), c.pot(pos, vy))
            elif dataMesh:
                cd = dataMesh.findCell(pos)
                if cd is None:
                    raise Exception("Cannot find " + str(pos) + " dataMesh")
                if len(vx) == dataMesh.cellCount():
                    d = pg.RVector3(vx[cd.id()], vy[cd.id()])
                elif len(vx) == dataMesh.nodeCount() and \
                        len(vy) == dataMesh.nodeCount():
                    d = pg.RVector3(cd.pot(pos, vx), cd.pot(pos, vy))
                else:
                    print(dataMesh, len(vx), len(vy))
                    raise Exception("data size wrong")
            else:
                print("mesh:", mesh, len(vx), len(vy))
                raise Exception("Data length neighter node size or cell size.")
        else:
            if dataMesh:
                cd = dataMesh.findCell(pos)
                if not cd:
                    break

                d = cd.grad(pos, pot)
                u = cd.pot(pos, pot)
            else:
                d = c.grad(pos, pot)
                u = c.pot(pos, pot)

        # always go u down
        dAbs = d.length()
        #print("cell:", c.id(), u, d, dAbs)

        if dAbs == 0.0:
            #print(d, "check this in streamlineDir(",)
            break

        if down:
            if u > lastU:
                break
        else:
            if u < lastU:
                break

        # * min(1.0, ((startCoord - pos).length()))
        pos += direction * d / dAbs * dLength
        c = mesh.findCell(pos, False)

        # Change cell here .. set old cell to be processed
        if c is not None:
            xd.append(pos[coords[0]])
            yd.append(pos[coords[1]])
            # set the stating value here
            if vd[0] == -1:
                vd[0] = dAbs
            vd.append(dAbs)

            ## check for degenerating stream
            if len(xd) > 2:
                pos0 = pg.Pos(xd[-3], yd[-3])
                pos1 = pg.Pos(xd[-2], yd[-2])
                pos2 = pg.Pos(xd[-1], yd[-1])
                if (pos0.dist(pos2) < pos0.dist(pos1)):
                    pg.debug('degenerating stream aborted')
                    break

            # If the new cell is different from the current we move into the
            # new cell and make the last to be invalid ..
            # the last active contains a stream element
            if c.id() != lastC.id():
                lastC.setValid(False)
                lastC = c
                dLength = c.center().dist(c.node(0).pos()) / dLengthSteps
        else:
            # There is no new cell .. the last active contains a stream element
            lastC.setValid(False)

        lastU = u
        if verbose:
            print(pos, u)

    # Stream line has stopped and the current cell (if there is one) ..
    # .. contains a stream element
    if c is not None:

        c.setValid(False)

    if down:
        xd.reverse(), yd.reverse(), vd.reverse()

    return xd, yd, vd
