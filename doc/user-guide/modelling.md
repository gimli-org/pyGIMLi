---
file_format: mystnb
kernelspec:
  name: python3
---

# Modelling

As the Finite Element analysis (FEA) is the most commonly used numerical method used,
this tutorial covers this method, referring the *M* (Modelling) in *pyGIMLi*.

## Finite Elements for Poisson equation

Details on the theory of the Finite Elements Analysis (FEA) can be found in several books,
e.g., {cite}`Zienkiewicz1977`.
Nevertheless, we provide a brief overview of the main concepts and ideas behind FEA to solve
boundary value problems.

Assume Poisson's equation as the simplest partial differential equation (PDE)
to be solved for the scalar field $ u(\mathbf{r}) $ within a modelling domain ${r}\in\Omega$
with a source $f$.

$$ - \nabla \cdot a \nabla u = f \quad{\mathrm{in}}~\Omega $$

<!-- $$\\ u =  \quad{\mathrm{on}}~\partial\Omega\ $$ -->

The generalized Poisson operator $\Delta = \nabla\cdot\nabla$, i.e., the divergence of the flow
(gradient times some conductivity), is a second-order partial derivative of the field
$u(\mathbf{r})$ in Cartesian coordinates, i.e., in 1D $ \mathbf{r} = (x) $,
in 2D $ {r} = (x, y) $, or in 3D space $ \mathbf{r} = (x, y, z) $.
On the boundary $\partial\Omega$ of the domain, we assume known values of $u=u_B$ (Dirichlet
boundary conditions) or gradients $\partial u/\partial n=g_B$.

A common approach to solve this problem is the method of weighted residuals.
An approximated solution $ u_h\approx u$ satisfies the PDE with a remainder $R=\Delta u_h + f$.
We choose some weighting functions $w$ and minimize $R$ over our modelling domain:

$$ \int_{\Omega} w R = 0\; $$

which leads to

$$ \int_{\Omega} - w \nabla \cdot a \nabla u_h = \int_{\Omega} w f $$

It is preferable to eliminate the second derivative in the Laplace operator,
either through to integration by parts or by applying the product rule and
Gauss's law. This leads to the so-called weak formulation:

$$  \int_{\Omega} a \nabla u_h \cdot \nabla w - \int_{\partial \Omega} a w \mathbf{n} \cdot \nabla u_h = \int_{\Omega} w f $$

$$ \int_{\Omega} a \nabla u_h \cdot \nabla w  = \int_{\Omega} w f + \int_{\partial \Omega} a w \frac{\partial u_h}{\partial\mathbf{n}} $$

We choose a function basis to approximate $u_h$:

$$ u_h = \sum_i^\mathcal{N} \mathrm{u}_i N_i $$

This fundamental relation discretizes the continuous solution $u_h$ by a vector of discrete coefficients
$\mathrm{u} = \{\mathrm{u}_i\} $ for a number $\mathcal{N}$ degrees of freedom (dof).
The basis functions $N_i$ can be understood as interpolation rules describing the solution on the whole domain, spanned by nodes, edges or faces in a mesh.
Mostly, $N_i$ represent nodal basis functions (e.g. hat functions) so that the $\mathrm{u}_i$ correspond to the solution at nodes.

Now we can set the unknown weighting functions to be identical to the basis
functions $w=N_j$ (Galerkin method) so that $ \forall j=0\ldots\mathcal{N} $ holds

$$ \int_{\Omega} a \sum_i \nabla N_i \cdot \nabla N_j \mathrm{u}_i =
   \int_{\Omega} f_j N_j + \int_{\partial \Omega} g N_j
   \quad \text{with}\quad g  = \frac{\partial u}{\partial \mathbf{n}} $$

For $g=0$ (natural Neumann boundary conditions) this can be rewritten as

$$ \mathrm{A} \mathrm{u} = \mathrm{b} $$

with $\mathrm{A}$ referred to as *Stiffness matrix*

$$ \mathrm{A} = \{\mathrm{a}_{i,j}\} = \int_{\Omega}\nabla N_i \cdot \nabla N_j $$

and the right hand $\mathrm{b}$ side known as *load vector*

$$ \mathrm{b} = \{\mathrm{b}_j\} = \int_{\Omega} f_j N_j $$

The solution of this linear system of equations leads to the
discrete solution $ \mathrm{u} = \{\mathrm{u}_i\} $ for all
$ i=1\ldots\mathcal{N}$ dofs spanning the modelling domain.

The choice of the dofs is crucial. If we choose too few, the accuracy of the sought solution might be pool. If we choose too many, the dimension of the system matrix will become large, leading to higher memory consumption and computation times.

To define the nodes, we discretize our modelling domain into cells, or the
eponymous elements. Cells are basic geometric shapes like triangles or
hexahedrons and are constructed from the nodes and collected in a mesh.
For more details, refer to the [Meshes](meshes.md) section.

To complete the solution for the small example, we still need to apply the boundary condition $u=u_B$ which is known as the Dirichlet condition. Setting
explicit values for our solution is not covered by the general Galerkin weighted
residuum method but we can solve it algebraically. We reduce the linear system
of equations by the known solutions $u_B={u_k}$ for all $k$ nodes on
the affected boundary elements:

$$ \mathrm{A_D}\cdot\mathrm{u} = \mathrm{b_D} \\
     \text{with } \mathrm{A_D} = \{\mathrm{a}_{i,j}\}\quad\forall i, j ~\notin~ k ~\text{and}~1~\forall~i,j \in k\\
    \mathrm{b_D}  = \{\mathrm{b}_j\} - \mathrm{A}\cdot\mathrm{g}\quad\forall j \notin k~\text{and}~u_k~\forall~j \in k  $$

Now we have all parts for assembling $\mathrm{A_D}$ and
$\mathrm{b_D}$ and finally solve the given boundary value problem for the interior points.

It is usually a good idea to test a numerical approach with known solutions.
To keep things simple, we create a modelling problem from the reverse direction (method of manufactured solution).
We choose a solution, calculate the right hand side function
and select the domain geometry suitable for nice Dirichlet values.

$$  u(x,y) = \operatorname{sin}(x)\operatorname{sin}(y) $$
$$ - \Delta u = f(x,y) = 2 \operatorname{sin}(x)\operatorname{sin}(y) $$
$$  \Omega \in I\!R^2  \quad \text{on}\quad 0 \leq x \leq 2\pi,~~  0 \leq y \leq 2\pi $$
$$ u  = g = 0 \quad \text{on}\quad \partial \Omega  $$

We now can solve the Poison equation applying the FEA capabilities of pygimli
and compare the resulting approximate solution $\mathrm{u}$
with our known exact solution $u(x,y)$.

## Parameterizing a mesh with physical properties

After importing the necessary modules

```{code-cell}
:tags: [hide-cell]

import numpy as np
import pygimli as pg

from pygimli.solver import solve
from pygimli.viewer import show
from pygimli.viewer.mpl import drawStreams
import pygimli.meshtools as mt
```

we utilize the `meshtools` module to generate a mesh with different regions that can be attributed by physical properties.

`createWorld` creates the definition for the modelling domain. `worldMarker=True` indicates the boundary conditions for the Earths surface and the subsurface. We assume layer boundaties at `y=-10` and `y=-30` so what we have a world with three different markers for the three layers.

```{code-cell}
world = mt.createWorld(start=[-50, 0], end=[50, -50], layers=[-10, -30],worldMarker=True)
```

We create two circular anomalies and assign the markers 4 and 5:

```{code-cell}
block_1 = mt.createCircle(pos=[-5, -3.], radius=[4, 1], marker=4,
                          boundaryMarker=10, area=0.1)
block_2  = mt.createCircle(pos=[10, -3.], radius=[4, 1], marker=5,
                          boundaryMarker=10, area=0.1)
```

The geometry definitions are merged into a Piecewise-Linear Complex (PLC) and plotted using `pg.show`, the keywords `markers=True` and `boundaryMarkers=True` show how the regions are numbered.

```{code-cell}
geom = world + block_1 + block_2
ax, cb = pg.show(geom, markers=True, boundaryMarkers=False)
```

Create a mesh for the finite element modelling with appropriate mesh quality.

```{code-cell}
mesh = mt.createMesh(geom, quality=34)
```

You can also print the amount of markers that are in the mesh by using the following commands:

```{code-cell}
number_of_cellmarkers = list(set(mesh.cellMarkers()))
print(number_of_cellmarkers)
```

Create a map to set resistivity values for the five regions
 [[regionNumber, resistivity], [regionNumber, resistivity], [...]

```{code-cell}
rhomap = [[1, 100.],
          [2, 75.],
          [3, 50.],
          [4, 150.],
          [5, 25]]
```

Here we assigned a different resistivity value to each part of the mesh using its cell markers and we can view the resistivity distribution with the following command:

```{code-cell}
ax, cb = pg.show(mesh, data=rhomap, label=pg.unit('res'), markers=True)
```
We will now show the different options to set boundary conditions and get these ready to then simulate and model the data.

## Boundary conditions (BC)

:::{admonition} Definition of Boundary Conditions
:class: tip

- Boundary marker (-1) : surface boundary conditions - `pg.core.MARKER_BOUND_HOMOGEN_NEUMANN`
- Boundary marker (-2) : mixed-boundary conditions - `pg.core.MARKER_BOUND_MIXED`
- Boundary marker ( >= 1 ) : no-flow boundaries
:::

As shown in [meshes section](meshes.md), pyGIMLi automatically assigns boundaries when using `mt.createWorld()`. However, you can assign BCs to different elements of your PLC or mesh.

There are different ways of specifying BCs. They can be maps from markers to values, explicit functions or implicit (lambda) functions. We use the example of the Poisson equation on the unit square and specify different boundary conditions on the four sides.

- The boundary 1 (left) and 2 (right) are directly mapped to the values 1 and 2.
- On side 3 (top) a lambda function 3+x is used (p is the boundary position and p[0] its x coordinate).
- On side 4 (bottom) a function uDirichlet is used that simply returns 4 in this example but can compute anything as a function of the individual boundaries b.

```{code-cell}
def uDirichlet(boundary):
    """Return a solution value for a given boundary.
        Scalar values are applied to all nodes of the boundary."""
    return 4.0

dirichletBC = {1: 1,                                           # left
               2: 2.0,                                         # right
               3: lambda boundary: 3.0 + boundary.center()[0], # bottom
               4: uDirichlet}                                  # top
```

The boundary conditions are passed using the BC keyword dictionary `dirichletBC`.

```{code-cell}
grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 21),
                     y=np.linspace(-1.0, 1.0, 21))
u = solve(grid, f=1., bc={'Dirichlet': dirichletBC})

# Note that showMesh returns the created figure ax and the created colorBar.
ax, cbar = show(grid, data=u, label='Solution $u$')

show(grid, ax=ax)

ax.text(1.02, 0, '$u=2$', va='center', ha='left',  rotation='vertical')
ax.text(-1.01, 0, '$u=1$', va='center', ha='right', rotation='vertical')
ax.text(0, 1.01, '$u=4$', ha='center')
ax.text(0, -1.01, '$u=3+x$', ha='center', va='top')

ax.set_title('$\\nabla\cdot(1\\nabla u)=1$')

ax.set_xlim([-1.1, 1.1])  # some boundary for the text
ax.set_ylim([-1.1, 1.1]);
```

Alternatively we can define the gradients of the solution on the boundary, i.e., Neumann type BC. This is done with another dictionary {marker: value} and passed by the bc dictionary.

```{code-cell}
neumannBC = {1: -0.5,  # left
             4: 2.5}  # bottom

dirichletBC = {3: 1.0}  # top

u = solve(grid, f=0., bc={'Dirichlet': dirichletBC, 'Neumann': neumannBC})
```

Note that on boundary 2 (right) has no BC explicitly applied leading to default (natural) BC that are of homogeneous Neumann type $\frac{\partial u}{\partial n}=0$

```{code-cell}
ax = show(grid, data=u, filled=True, orientation='vertical',
          label='Solution $u$',
          levels=np.linspace(min(u), max(u), 14), hold=True)[0]

# Instead of the grid we now want to add streamlines to show the gradients of
# the solution (i.e., the flow direction).
drawStreams(ax, grid, u)
ax.text(0.0, 1.01, '$u=1$',
        horizontalalignment='center')  # top -- 3
ax.text(-1.0, 0.0, '$\partial u/\partial n=-0.5$',
        va='center', ha='right', rotation='vertical')  # left -- 1
ax.text(0.0, -1.01, '$\partial u/\partial n=2.5$',
        ha='center', va='top')  # bot -- 4
ax.text(1.01, 0.0, '$\partial u/\partial n=0$',
        va='center', ha='left', rotation='vertical')  # right -- 2

ax.set_title('$\\nabla\cdot(1\\nabla u)=0$')

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1]);
```

Its also possible to force single nodes to fixed values too.

```{code-cell}
u = solve(grid, f=1., bc={'Node': [grid.findNearestNode([0.0, 0.0]), 1.0]})
np.testing.assert_approx_equal(u[grid.findNearestNode([0.0, 0.0])], 1.0, significant=10)

ax, _ = pg.show(grid, u, logScale=False, label='Solution $u$',)
_ = pg.show(grid, ax=ax)
```
