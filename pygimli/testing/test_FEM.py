#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

import pygimli as pg
import numpy as np

class TestFiniteElementBasics(unittest.TestCase):
    def test_Poisson(self):
        """
        \d² u / d x² + f = 0

        u = sin(x) ## maybe a good test case for harmonic (spectral) base functions TODO
        f = sin(x)

        u = -cos(x)
        f = -cos(x)

        ### the following are already tested in Helmholtz
        u = x
        f = 0

        u = x*x
        f = -2
        """
        pass

    def test_Helmholtz(self):
        """
            d² u / d x² + k u + f = 0
            k = 2
            a) P1(exact)
                u = x
                f = -k x
            b) P2(exact)
                u = x*x
                f = -(2 + 2x*x)
        """
        h = np.pi/2 / 21
        x = np.arange(0.0, np.pi/2, h)
        mesh = pg.createGrid(x)

        ### test a)
        k = 2.0
        u = lambda _x: _x
        f = lambda _x, _k: -_k * u(_x)

        x = pg.x(mesh)
        dirichletBC = {1: u(min(x)), 2: u(max(x))}
        uFEM = pg.solve(mesh, a=1, b=k, f=f(x, k),
                        bc={'Dirichlet': dirichletBC})

        # pg.plt.plot(x, uFEM, '.')
        # pg.plt.plot(pg.sort(x), u(pg.sort(x)))
        # pg.wait()

        np.testing.assert_allclose(uFEM, u(x))

        ### test b)
        u = lambda _x: _x * _x
        f = lambda _x, _k: -(2. + k * u(_x))

        mesh = mesh.createP2()
        x = pg.x(mesh)
        dirichletBC = {1: u(min(x)), 2: u(max(x))}
        uFEM = pg.solve(mesh, a=1, b=k, f=f(x, k),
                        bc={'Dirichlet': dirichletBC})
        np.testing.assert_allclose(uFEM, u(x), atol=1e-6)

        # pg.plt.plot(x, uFEM, '.')
        # pg.plt.plot(pg.sort(x), u(pg.sort(x)))
        # pg.wait()

    def test_Neumann_BC(self):
        def _test_(mesh, p2=False, show=False):
            """
                \Laplace u = 0
                x = [0, -1], y, z
                int_0^1 u = 0

                du/dn =  1 (xMin) (inflow on -x)
                du/dn = -1 (xMax) (outflow on +x)
                du/dn = 0 (rest)
                u = 0.5 -x          linear solution, i.e., exact with p1

                du/dn = -1 (xMin) (outflow -x)
                du/dn = -1 (xMax) (outflow +x)
                du/dn = 0 (rest)
                u = -1/6 + x -x²    quadratic solution, i.e., exact with p2
            """
            uExact = lambda x, a, b, c: a + b *x + c * x**2
            bc={'Neumann': {1: 1.0, 2: -1.0}}
            uE = uExact(pg.x(mesh), 0.5, -1.0, 0.0)

            if p2 is True:
                mesh = mesh.createP2()
                bc['Neumann'][1] = -1.0
                uE = uExact(pg.x(mesh), -1/6, 1.0, -1.0)

            u = pg.solve(mesh, a=1, f=0, bc=bc, verbose=0)
            v = pg.solver.grad(mesh, u)

            if show:
                if mesh.dim() == 1:
                    idx = np.argsort(pg.x(mesh))
                    fig, ax = pg.plt.subplots()
                    ax.plot(pg.x(mesh)[idx], uE[idx], label='exact')
                    ax.plot(pg.x(mesh)[idx], u[idx], 'o', label='FEM')

                    model, response = pg.frameworks.fit(uExact,
                                                    u[idx], x=pg.x(mesh)[idx])
                    print(model)
                    ax.plot(pg.x(mesh)[idx], response, 'x', label='FIT')
                    ax.grid(True)
                    ax.legend()
                elif mesh.dim() == 2:
                    pg.show(mesh, u, label='u')
                    ax, _ = pg.show(mesh, pg.abs(v), label='abs(grad(u))')
                    pg.show(mesh, v, showMesh=1, ax=ax)

                pg.info("int Domain:", pg.solver.intDomain(u, mesh))
                pg.info("int Domain:", pg.solver.intDomain([1.0]*mesh.nodeCount(), mesh), sum(mesh.cellSizes()))

                pg.wait()
            ## test du/dn of solution and compare with Neumann BC
            for m, val in bc['Neumann'].items():

                for b in mesh.boundaries(mesh.boundaryMarkers() == m):
                    ## for non Tailor Hood Elements, the gradient is only
                    # known at the cell center so the accuracy for the
                    # gradient depends on the distance boundary to cell
                    # center. Accuracy = du/dx(dx) = 1-2x = 2 * dx
                    c = b.leftCell()
                    dx = c.center().dist(b.center())
                    dudn = b.norm(c).dot(v[c.id()])
                    if p2:
                        np.testing.assert_allclose(dudn, val, atol=2 * dx)
                        #print(dudn, val)
                    else:
                        np.testing.assert_allclose(dudn, val)

            np.testing.assert_allclose(pg.solver.intDomain([1.0]*\
                                                    mesh.nodeCount(), mesh),
                                       sum(mesh.cellSizes()))
            np.testing.assert_allclose(pg.solver.intDomain(u, mesh), 0.0,
                                       atol=1e-8)
            np.testing.assert_allclose(np.linalg.norm(u-uE), 0.0,
                                       atol=1e-8)

            return v

        # 1D
        x = np.linspace(0, 1, 101)
        _test_(pg.createGrid(x=x), p2=False, show=False)
        _test_(pg.createGrid(x=x), p2=True, show=False)

        # # 2D grid
        x = np.linspace(0, 1, 11)
        _test_(pg.createGrid(x=x, y=x), p2=False, show=False)
        _test_(pg.createGrid(x=x, y=x), p2=True, show=False)

        # 2D tri
        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(
                                start=[0, 0], end=[1, 1], worldMarker=False),
                                        area=0.05)

        _test_(mesh, p2=False, show=False)
        _test_(mesh, p2=True, show=False)

        # 3D prism
        mesh = pg.meshtools.extrudeMesh(mesh, np.linspace(0, 1, 11))
        _test_(mesh, p2=False)
        _test_(mesh, p2=True)

        # 3D quads
        x = np.linspace(0, 1, 5)
        _test_(pg.createGrid(x=x, y=x, z=x), p2=False)
        #_test_(pg.createGrid(x=x, y=x, z=x), p2=True) ## fails .. need check!

        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(
                                       start=[0, 0, 0], end=[1, 1, 1], worldMarker=False),
                                       area=0.1)

        # 3D tet
        _test_(mesh, p2=False)
        _test_(mesh, p2=True)

        # pg.show(mesh)
        # pg.wait()

    def test_Dirichlet_BC(self):
        """
        """
        def _testP2_(mesh, show=False):
            """ Laplace u - 2 = 0
                with u(x) = x² + x/(xMax-xMin) + xMin/(xMax-xMin)
                and u(xMin)=0 and u(xMax)=1
                Test for u_h === u(x) for P2 base functions
            """
            meshP2 = mesh.createP2()
            u = pg.solve(meshP2, f=-2, bc={'Dirichlet': {'1,2': 2}},)

            xMin = mesh.xMin()
            xSpan = (mesh.xMax() - xMin)

            if show:
                if mesh.dim()==1:
                    pg.plt.figure()
                    x = pg.x(mesh)
                    ix = np.argsort(x)
                    pg.plt.plot(x[ix], x[ix]**2 - (xSpan/2)**2 + 2)
                elif mesh.dim() > 1:
                    pg.show(meshP2, u, label='u = x**2')

           
            uE = pg.x(meshP2)**2 - (xSpan/2)**2 +2
            np.testing.assert_allclose(u, uE)
            np.testing.assert_allclose(0.0, np.linalg.norm(u-uE), atol=1e-8)

        def _testP1_(mesh, show=False, followP2=True):
            """ Laplace u = 0
                with u = x and u(x=min)=0 and u(x=max)=1
                Test for u == exact x for P1 base functions
            """
            #print('b', mesh, mesh.cell(0))
            u = pg.solve(mesh, a=1, b=0, f=0,
                         bc={'Dirichlet': {1: 0, 2: 1}})

            if show:
                if mesh.dim()==1:
                    pg.plt.figure()
                    x = pg.x(mesh)
                    ix = np.argsort(x)
                    pg.plt.plot(x[ix], u[ix])
                elif mesh.dim() > 1:
                    pg.show(mesh, u, label='u')

            xMin = mesh.xMin()
            xSpan = (mesh.xMax() - xMin)
            np.testing.assert_allclose(u, (pg.x(mesh)-xMin) / xSpan)

            if followP2:
                _testP2_(mesh, show)
            return u


        # 1D
        _testP1_(pg.createGrid(x=np.linspace(-1, 1, 11)), show=False)
        _testP1_(pg.createGrid(x=np.linspace(-3, 3, 11)), show=False)

        # 2D reg quad
        _testP1_(pg.createGrid(x=np.linspace(-2, 2, 11),
                               y=np.linspace(0, 1, 11)), show=False)
        # 2D reg quad scale
        _testP1_(pg.createGrid(x=np.linspace(-0.04, 0.04, 11),
                               y=np.linspace(-0.4, 0, 11)))
        # 3D reg quad
        # _testP1_(pg.createGrid(x=np.linspace(-2, 2, 11),
        #                        y=np.linspace( 0, 1, 11),
        #                        z=np.linspace( 0, 1, 11)),
        #                        followP2=False
        #                        )
        #return
        # check why P2 fails here,

        # 2D tri
        mesh = pg.meshtools.createMesh(pg.meshtools.createWorld(start=[-1, -4],
                                                                end=[1, -6], worldMarker=0), area=0.1)
        # mesh.setBoundaryMarkers(np.array([0,1,3,2,4])[mesh.boundaryMarkers()])
        _testP1_(mesh, show=False)

        # 3D prism
        mesh = pg.meshtools.extrudeMesh(mesh, np.linspace(0, 1, 11))
        _testP1_(mesh, show=False)

        # 3D tet
        mesh = pg.meshtools.createMesh(pg.meshtools.createCube(size=[4, 4, 4],
                                                               boundaryMarker=9,
                                                               area=1.1))

        for b in mesh.boundaries(mesh.boundaryMarkers() == 9):

            if b.norm()[0] == -1:
                b.setMarker(1)
            elif b.norm()[0] == +1:
                b.setMarker(2)

        _testP1_(mesh, show=False)
        #pg.wait()

        # grid = pg.createGrid(x=np.linspace(-2, 2, 11), y=np.linspace(0, 1, 11))
        # grid.rotate([0, 0, np.pi/4])
        #can't find proper test for this rotated case
        #v = _testP1_(grid, show=False) #2D reg - rotated

    def test_Robin_BC(self):
        """
            Linear function should result in exact solution with linear base
            u(x) = Ax + B on x = [0, .. ,1]

            Dirichlet BC: u(x) = g
            ----------------------
            u(0) = B
            u(1) = A + B

            Neumann BC: du/dn(x) = f (non unique, needs calibration point)
            ------------------------------------------------------------
            du/dn(0) = -A (n = [-1, 0, 0])
            du/dn(1) =  A (n = [ 1, 0, 0])

            Robin BC v1: du/dn(x) = a(u0-u(x))
            ----------------------------------
            du/dn(0) = -A = a (u0 - B)     (a=1, u0=-A/a +B)
            du/dn(1) =  A = a (u0 - (A+B)) (a=1, u0= A/a +A+B)

            Robin BC v2: b du/dn(x) + a u(x) = g
            ------------------------------------
            du/dn(0) = -A + a B = g (a=1, b=1.0, g=-A+a*B)
            du/dn(1) =  A + a (A+B) = g  (a=1, b=1.0, g= A + a*(A+B))
        """
        return 
        show=False
        A = 1
        B = 3
        x = np.linspace(0, 1, 11)
        uExact = lambda x: A*x + B
        if show:
            pg.plt.plot(x, uExact(x), label='uExact')

        mesh = pg.createGrid(x=x)
        u = pg.solve(mesh, bc={'Dirichlet':{1:B, 2:A+B}})
        np.testing.assert_allclose(u, uExact(x))

        if show:
            pg.plt.plot(x, u, 'o', label='u Dirichlet')

        n = mesh.findNearestNode([0.2, 0.0])
        u = pg.solve(mesh, bc={'Node': [n, uExact(mesh.node(n).x())],
                               'Neumann':{1:-A, 2:A}})
        np.testing.assert_allclose(u, uExact(x))

        if show:
            pg.plt.plot(x, u, '.', label='u Neumann')

        a = 0.5
        u1 = -A/a +B
        u2 =  A/a +A+B

        u = pg.solve(mesh, bc={'Robin':{1: [a, u1],
                                        2: [a, u2]},})
        np.testing.assert_allclose(u, uExact(x))
        if show:
            pg.plt.plot(x, u, '-^', label='u Robin-v1')

        al = 0.5
        ga1 = -A + al * B
        ga2 =  A + al * (A+B)
        u = pg.solve(mesh, bc={'Robin':{1: [al, 1.0, ga1],
                                        2: [al, 1.0, ga2],
                                        }
                                })
        np.testing.assert_allclose(u, uExact(x))
        if show:
            pg.plt.plot(x, u, 'v', label='u Robin-v2')
            pg.plt.legend()


    def testElementMatrix(self):
        a = pg.core.ElementMatrix()


if __name__ == '__main__':

    # test = TestFiniteElementBasics()
    # test.test_Dirichlet_BC()
    # test.test_Neumann_BC()
    # test.test_Robin_BC()

    unittest.main()
