r"""
Basics of Finite Element Analysis
---------------------------------

This tutorial covers the first steps into Finite Element computation
referring the *M* (Modelling) in *pyGIMLi*.

We will not dig into deep details about the theory of the Finite Elements Analysis (FEA)
here, as this can be found in several books, e.g., :cite:`Zienkiewicz1977`.

Anyhow, there is a little need for theory to understand what it means
to use FEA for the solution of a boundary value problem.
So we start with some basics.

Assuming Poisson's equation as a simple partial differential problem
to be solved for the sought scalar field :math:`u(\mathbf{r})` within
a modelling domain :math:`\mathbf{r}\in\Omega`
with a non zero right hand side function :math:`f`.

.. math::

   - \Delta u & = f \quad{\mathrm{in}}\quad~\Omega\\
            u & = g \quad{\mathrm{on}}\quad\partial\Omega\;.

The Laplace operator :math:`\Delta = \nabla\cdot\nabla` given by the divergence
of the gradient, is the sum of the second partial derivatives of the field
:math:`u(\mathbf{r})` with respect to the Cartesian coordinates
in 1D space :math:`\mathbf{r} = (x)`, in 2D :math:`\mathbf{r} = (x, y)`, or 3D
space :math:`\mathbf{r} = (x, y, z)`.
On the boundary :math:`\partial\Omega` of the domain, we want
known values of :math:`u=g` as so called Dirichlet boundary conditions.

A common approach to solve this problem is the method of weighted residuals.
The base assumption states that an approximated solution :math:`u_h\approx u` will
only satisfy the differential equation with a rest :math:`R`: :math:`\Delta u_h + f = R`. If we choose some weighting functions :math:`w`,
we can try to minimize the resulting residuum over our modelling domain as:

.. math::

    \int_{\Omega} R w = 0\;,

which leads to:

.. math::

    \int_{\Omega} - \Delta u_h w = \int_{\Omega} f w\;.

It is preferable to eliminate the second derivative in the Laplace operator,
either due to integration by parts or by applying the product rule and
Gauss's law.
This leads to the so called weak formulation:

.. math::

    \int_{\Omega} \nabla u_h \nabla w - \int_{\partial \Omega}\mathbf{n}\nabla u_h w & = \int_{\Omega} f w \\
    \int_{\Omega} \nabla u_h \nabla w & = \int_{\Omega} f w + \int_{\partial \Omega}\frac{\partial u_h}{\partial\mathbf{n}} w\;.

We can solve these integrals after choosing an appropriate basis
for an approximate solution :math:`u_h` as:

.. math::

    u_h = \sum_i \mathrm{u}_i N_i\quad\text{with}\quad i = 0\ldots\mathcal{N}\;.

The latter fundamental FEA relation discretizes the continuous solution :math:`u_h`
into a discrete values :math:`\mathrm{u} = \{\mathrm{u}_i\}` for a number of :math:`i = 0\ldots\mathcal{N}` discrete points, usually called nodes.

The basis functions :math:`N_i` can be understood as interpolation rules for
the discrete solution between adjacent nodes and will be chosen later.

Now we can set the unknown weighting functions to be the same as the basis
functions :math:`w=N_j` with :math:`j=0\ldots\mathcal{N}` (Galerkin method)

.. math::

    \int_{\Omega} \sum_i \mathrm{u}_i \nabla N_i \nabla N_j\\ & = \int_{\Omega} f_j N_j + \int_{\partial \Omega} h N_j
    \quad \text{with}\quad h\\ & = \frac{\partial u}{\partial \mathbf{n}}

this can be rewritten with :math:`h=0` as:

.. math::

    \mathrm{A} \mathrm{u} &= \mathrm{b} \\
    & \text{with} \\
    \mathrm{A} & = \{\mathrm{a}_{i,j}\} = \int_{\Omega}\nabla N_i \nabla N_j \quad\text{known as 'Stiffness matrix'}\\
    \mathrm{b} & = \{\mathrm{b}_j\} = \int_{\Omega} f_j N_j \quad\text{known as 'Load vector'}

The solution of this linear system of equations leads to the
discrete solution :math:`\mathrm{u} = \{\mathrm{u}_i\}` for all
:math:`i=1\ldots\mathcal{N}` nodes inside the modelling domain.

For the practical part, the choice of the nodes is crucial. If we choose too
little, the accuracy of the sought solution might be too small. If we choose too
many, the dimension of the system matrix will be too large, which leads to
higher memory consumption and calculation times.

To define the nodes, we discretize our modelling domain into cells, or the
eponymous elements. Cells are basic geometric shapes like triangles or
hexahedrons and are constructed from the nodes and collected in a mesh. See the
tutorials about the mesh basics (:ref:`tut:basics`). In summary, the discrete
solutions of the differential equation using FEA on a specific mesh are defined
on the node positions of the mesh.

The chosen mesh cells also define the base functions and the integration rules
that are necessary to assemble the stiffness matrix and the load vector and will
be discussed in a different tutorial (TOWRITE link here).

To finally solve our little example we still need to handle the application of
the boundary condition :math:`u=g` which is called Dirichlet condition. Setting
explicit values for our solution is not covered by the general Galerkin weighted
residuum method but we can solve it algebraically. We reduce the linear system
of equations by the known solutions :math:`g={g_k}` for all :math:`k` nodes on
the affected boundary elements: (maybe move this to the BC tutorial)

.. math::

    \mathrm{A_D}\cdot\mathrm{u} &= \mathrm{b_D} \\
    & \text{with} \\
    \mathrm{A_D} & = \{\mathrm{a}_{i,j}\}\quad\forall i, j ~\notin~ k ~\text{and}~1~\text{for}~i,j \in k\\
    \mathrm{b_D} & = \{\mathrm{b}_j\} - \mathrm{A}\cdot\mathrm{g}\quad\forall j \notin k~\text{and}~g_k~\text{for}~j \in k

Now we have all parts together to assemble :math:`\mathrm{A_D}` and
:math:`\mathrm{b_D}` and finally solve the given boundary value problem.

It is usually a good idea to test a numerical approach with known solutions.
To keep things simple we create a modelling problem from the reverse direction.
We choose a solution, calculate the right hand side function
and select the domain geometry suitable for nice Dirichlet values.

.. math::

    u(x,y) & = \operatorname{sin}(x)\operatorname{sin}(y)\\
    - \Delta u & = f(x,y) = 2 \operatorname{sin}(x)\operatorname{sin}(y)\\
    \Omega \in I\!R^2 & \quad \text{on}\quad 0 \leq x \leq 2\pi,~~  0 \leq y \leq 2\pi \\
    u & = g = 0 \quad \text{on}\quad \partial \Omega

We now can solve the Poison equation applying the FEA capabilities of pygimli
and compare the resulting approximate solution :math:`\mathrm{u}`
with our known exact solution :math:`u(x,y)`.
"""

import numpy as np
import pygimli as pg

###############################################################################
# We start to define the modelling domain and functions for the exact solution
# and the values for the load vector.
# The desired mesh of our domain will be a grid with equidistant spacing in
# x and y directions.
#
domain = pg.createGrid(x=np.linspace(0.0, 2*np.pi, 25),
                       y=np.linspace(0.0, 2*np.pi, 25))

uExact = lambda pos: np.sin(pos[0]) * np.sin(pos[1])
f = lambda cell: 2.0 * np.sin(cell.center()[0]) * np.sin(cell.center()[1])

###############################################################################
# We use the existing shortcut functions for the assembling of the basic FEA
# system matrices and vectors.
# The implemented parts of the solving process are supposed
# to be dimension independent. You only need to find a valid mesh with the
# supported element types.
#
A = pg.solver.createStiffnessMatrix(domain)
b = pg.solver.createLoadVector(domain, f)

###############################################################################
# To apply the boundary condition we first need to identify all boundary
# elements. The default grid applies the following boundary marker on the
# outermost boundaries: 1 (left), 2(right), 3(top), and 4(bottom).
boundaries = pg.solver.parseArgToBoundaries({'1,2,3,4': 0.0}, domain)

###############################################################################
# `parseArgToBoundaries` is a helper function to collect a list of
# tupels (Boundary element, value), which can be used to apply the Dirichlet
# conditions.
#
pg.solver.assembleDirichletBC(A, boundaries, b)

###############################################################################
# The approximate solution :math:`\mathrm{u}` can then be found as the
# solution of the linear system of equations.
#
u = pg.solver.linSolve(A, b)

###############################################################################
# The resulting scalar field can displayed with the `pg.show` shortcut.
#
ax, cb = pg.show(domain, u, label='Approximated solution $u$', nLevs=7)

###############################################################################
# For analyzing the accuracy for the approximation we apply the
# L2 norm for the finite element space :py:mod:`pygimli.solver.normL2` for a
# set of different solutions with decreasing cell size. Instead of using the
# the single assembling steps again, we apply our Finite Element shortcut function
# :py:mod:`pygimli.solver.solve`.
#
domain = pg.createGrid(x=np.linspace(0.0, 2*np.pi, 3),
                       y=np.linspace(0.0, 2*np.pi, 3))

h = []
l2 = []
for i in range(5):
    domain = domain.createH2()
    u_h = pg.solve(domain, f=f, bc={'Dirichlet':{'1:5': 0}})
    u = np.array([uExact(_) for _ in domain.positions()])
    l2.append(pg.solver.normL2(u - u_h, domain))
    h.append(min(domain.boundarySizes()))
    print("NodeCount: {0}, h:{1}m, L2:{2}%".format(domain.nodeCount(),
                                                   h[-1], l2[-1]))

ax,_ = pg.show()
ax.loglog(h, l2, 'o-')
ax.set_ylabel('Approximation error: $L_2$ norm')
ax.set_xlabel('Cell size $h$ (m)')
ax.grid()

###############################################################################
# We calculated the examples before for a homogeneous material parameter a=1,
# but we can apply any heterogeneous values to0. One way is to create a list of
# parameter values, one for each cell of the domain. Currently the values for
# each cell can be of type float, complex, or real valued anisotropy or
# constitutive matrix. For illustration we show a calculation with an
# anisotropic material. We simply use the same setting as above and assume a
# -45 degree dipping angle in the left and 45 degree dipping in the right part
# of the domain. Maybe we will find someday a more meaningful example. If you
# have an idea please don't hesitate to share.
#
a = [None]*domain.cellCount()
for c in domain.cells():
    if c.center()[0] < np.pi:
        a[c.id()] = pg.solver.createAnisotropyMatrix(lon=1.0, trans=10.0,
                                               theta=-45/180 * np.pi)
    else:
        a[c.id()] = pg.solver.createAnisotropyMatrix(lon=1.0, trans=10.0,
                                               theta=45/180 * np.pi)

u = pg.solve(domain, a=a, f=f, bc={'Dirichlet':{'*': 0}})
ax, cb = pg.show(domain, u, label='Solution $u$ for anisotropic $a$', nLevs=7)
