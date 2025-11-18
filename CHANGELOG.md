# Changelog

##
## v1.1 (06/05/2020)
Release presented at EGU2020: https://doi.org/10.5194/egusphere-egu2020-18751

``` bash
conda install -c gimli -c conda-forge pygimli=1.1 # Win, Linux and Mac!
```

#### Main changes:
- Improved and cleaned frameworks: https://www.pygimli.org/pygimliapi/_generated/pygimli.frameworks.html
- `pygimli.physics.ert.ERTManager` with full BERT functionality
- Repository restructuring (Python part now at root level)
- 3D viewer based on PyVista: https://www.pygimli.org/pygimliapi/_generated/pygimli.viewer.pv.html
- Improved API documentation and website

##### API renaming & namespace clean-up
- `pg.mplviewer` - > `pg.viewer.mpl`
- `pg.RVector` -> `pg.Vector`
- `pg.RMatrix` -> `pg.Matrix`
- New module `pg.math` (cos, cot, sin, det, exp10, median, randn, rms, etc.)
- New module `pg.matrix` (Vector, BlockMatrix, etc.)

##### Default changes in viewer
- `logScale=False`
- `colorBar=True`

---

## v1.0.12 (28/10/2019)
``` bash
conda install -c gimli -c conda-forge pygimli=1.0.12 # Win & Linux!
```

#### bug

- [**bug**] Problem using fatray [#188](https://github.com/gimli-org/gimli/issues/188)
- [**bug**] pg.solver.solveFiniteElements(mesh, a=1.0, b=0...) does not work if b is set as a field [#181](https://github.com/gimli-org/gimli/issues/181)

---

## v1.0.11 (03/04/2019)
``` bash
conda install -c gimli -c conda-forge pygimli=1.0.11 # currently Linux only
````

#### bug

- [**bug**] Examples Raypaths, Layered Example not working [#167](https://github.com/gimli-org/gimli/issues/167)

#### closed

- [**closed**] Trying the Field data inversion (“Koenigsee”) example [#164](https://github.com/gimli-org/gimli/issues/164)
- [**closed**] Question -ModuleNotFoundError: No module named 'pygimli' [#160](https://github.com/gimli-org/gimli/issues/160)

---

## v1.0.10 (04/01/2019)
### Changes

- improved handling of region markers during mesh creation (PR #157 by @m-weigand)
- fixed marker bug for coverage mapping (https://github.com/gimli-org/gimli/commit/145dbad348e61908bbde37beec6175b2bf5f15f5)
- Matplotlib 3.0 support

### Closed issues

#### bug

- [**bug**] drawField does not work for very small values [#136](https://github.com/gimli-org/gimli/issues/136)
- [**bug**] potential wrong reference to const & RMatrix() const in python bindings [#61](https://github.com/gimli-org/gimli/issues/61)

#### closed

- [**closed**] Parameter check in createRectangle [#158](https://github.com/gimli-org/gimli/issues/158)
- [**closed**] Mesh creation: marker parameter not always honored [#152](https://github.com/gimli-org/gimli/issues/152)

#### feature request

- [**feature request**] Support vtk STRUCTURED_GRID Dataset [#87](https://github.com/gimli-org/gimli/issues/87)

### Misc

#### fixes
- polyCreate* that does not use leftDirection. 
- wrong Neumann BC scaling, 

#### adds
- dataContainerERT.addFourPointData accept now also values
- pygimli/core/datacontainer.py for easy extentions
- createCylinder(OOC: polytools), createMesh() accecpt 3D PLC and now tries a systemcall for tetgen
- dataIndex(), dataContainer.sortSensorIndex() returns now sorting indieces.
- isGeometry flag for core.Mesh. If set meshtool automatic check for duplicated nodes.
- BlockMatrix to sparseMap conversion
- argument for solver.linSolver to change solving backend

#### doc
- slight changes to interpolation api docs































---

## v1.0.9 (30/10/2018)
``` bash
conda install -c gimli pygimli=1.0.9 # currently Linux only
````

#### New functionality
- **Support for secondary nodes in traveltime calculations** (following Giroux & Larouche, 2013; http://dx.doi.org/10.1016/j.cageo.2012.12.005)
    - New mesh method `mesh.createSecondaryNodes(numberOfSecNodes)`
    - New argument when setting mesh in RefractionManager `rst.setMesh(mesh, numberOfSecNodes)`

- New method to visualize raypaths: `rst.showRayPaths()` (https://github.com/gimli-org/gimli/commit/a89eb2c620b3e211e04d63961b4067faedf8b323, [API Doc](https://www.pygimli.org/pygimliapi/_generated/pygimli.physics.traveltime.html#pygimli.physics.traveltime.Refraction.showRayPaths))

#### Other changes
- Fix RValue conversion bug (https://github.com/gimli-org/gimli/commit/172c31f063b8818aa7f195e0432412c13adff593)
- Automatic detection of inline backend (improves `pg.show()` usage in Jupyter Notebooks) (https://github.com/gimli-org/gimli/commit/207680e86b8a521c735684f673767974202abfc6)
- Allow for `cmap` and other `kwargs` in `rst.showVA()` (https://github.com/gimli-org/gimli/commit/096d02db953f4c1e0224fdefb478c97160de26b2)
- Fix compatibility issue with pytest 3.8 (https://github.com/gimli-org/gimli/commit/ba5930794ff8f6223e9a429377062e874732c830)


---

## v1.0.8 (27/09/2018)
#### New functionality

- automatic type conversion for Numpy scalars (https://github.com/gimli-org/gimli/commit/34041b03241f51d6ee16af144c70acb2e5596dce)
- automatic type conversion numpy.ndarray to BVector
- expose abs(RVector), abs(R3Vector) and abs(CVector). No need for pg.abs anymore.
- expose log(RVector). No need for pg.log anymore.
- add STL reading support for multiple soldis
- add basic operatores +/- for IVector

#### Closed

- [**closed**] pygimli.frameworks.harmfit does not work with numpy integers for nc [#143](https://github.com/gimli-org/gimli/issues/143)

#### Changes
- mt.mergePLC recognize node on edge touching (not vice versa)
- show now checks for used 2d coordinates (x,y) or (x,z)
- mesh.createNodeWithCheck can now check for edge touching and split the edge on hit
- renamed: mesh.extract to mesh.createSubMesh



---

## v1.0.7 (04/08/2018)

#### bug

- [**bug**] Zero or one-based sensor indexing for traveltime data container? [#141](https://github.com/gimli-org/gimli/issues/141)

#### closed

- [**closed**] Retrieve pyGIMLi version is 'untagged' [#137](https://github.com/gimli-org/gimli/issues/137)
- [**closed**] marker view with pg.show() [#135](https://github.com/gimli-org/gimli/issues/135)

#### fixes:
- pure singe region now build correct interregion constraints without region file
- passive body cem problem in bert fop
- *Vector < operator
- mesh bounding box problem for post createNode

### add:
- __repr__ for pg.Node.
- pg.load for *.collect files
- circular flag for patchValMap
- DataContainer.remove(BVector)
- showNodes flag to drawPLC  
- boundingbox syntax sugar

### remove:
- DataContainer.filter


---

## v1.0.6 (03/04/2018)
#### new functionality
- `pg.show(geom, markers=True)` to inspect region and boundary markers of PLCs (also works for meshes)
- `mt.interpolateAlongCurve`, add simple extrapolation for to little curve bounds.
- `pg.warn(*args), pg.info(*args), pg.debug(*args)` etc. easier logging interface #14
- `Pos=Vector3` alias
- `mt.createGrid`  proxy function
- static `ndim` property to `pg.Vector`, `pg.Matrix` etc.
- convenience calls: `pg.x`, `pg.y`, and `pg.z` to return appropriate arrays for (data, mesh, R3Vector, ndarray, ..)
- core: improve mesh entities selection via (BVector), e.g., `mesh.cells(mesh.cellMarkers() == 1)` 
is the same as `mesh.cells(pg.find(mesh.cellMarkers() == 1))`
- mplviewer new flag: showBoundary
- `pg.getConfigPath()` returns path to user config
- add more flexible way to define callables for boundary conditions

#### other changes
- `pg.test()` compatible with numpy 1.14 (and 1.13)
- documentation compatible with Sphinx > 1.6 (new `-i` flag in sphinx-autogen)
- increased visibility of hyperlinks in documentation
- removed io/ path #124

#### backward compatibility

- [**backward compatibility**] Backward compatibility of pg.interpolate [#131](https://github.com/gimli-org/gimli/issues/131)

#### bug

- [**bug**] Plotting issue for multiple subplots in Jupyter Notebook [#125](https://github.com/gimli-org/gimli/issues/125)
- fix implementation of Neumann boundary conditions.

#### enhancement

- [**enhancement**] pg.meshtools.cellDataToNodeData does not take scalar fields [#120](https://github.com/gimli-org/gimli/issues/120)
- CellBrowser: Highlight style, click now toggle, data and mesh can be changed

#### closed issues

- [**closed**] Drawing a colorbar with pg.mplviewer.drawModel [#130](https://github.com/gimli-org/gimli/issues/130)
- [**closed**] 'MemoryError' for pb.invert [#129](https://github.com/gimli-org/gimli/issues/129)
- [**closed**] FDEM dataset [#122](https://github.com/gimli-org/gimli/issues/122)
- [**closed**] petro.resistivityArchie interpolation returns 'corrupt mesh' error [#107](https://github.com/gimli-org/gimli/issues/107)
- [**closed**] API doc not building with sphinx 1.6.3 [#74](https://github.com/gimli-org/gimli/issues/74)
---

## v1.0.5 (22/12/2017)
Install with Ana- or Miniconda:

``` bash
conda install pygimli=v1.0.5
```
#### Fixes:
- automatic RValue conversion for scalar numpy data types to pg.core long 
- Crank-Nicolson solver is now working again with default theta setting

#### building and distribution

- [**building and distribution**] Make setup.py and conda's meta.yaml aware of pg.__version__ [#110](https://github.com/gimli-org/gimli/issues/110)

---

## v1.0.4 - Start of consequent versioning (22/12/2017)
Install with Ana- or Miniconda:

``` bash
conda install pygimli=v1.0.4
```
