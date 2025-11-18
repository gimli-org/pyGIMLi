# encoding: utf-8
"""
Mesh generation and modification.

.. note::

    Although we discriminate here between grids (structured meshes) and meshes
    (unstructured), both objects are treated the same internally.
"""

from pygimli.core import (createMesh1D, createMesh1DBlock,
                          createMesh2D, createMesh3D)

from .grid import (createGrid, createGridPieShaped,
                   appendBoundary,
                   appendBoundaryGrid,
                   appendTriangleBoundary,  # why is this in grid.py?
                   appendTetrahedronBoundary)

from .mapping import (cellDataToBoundaryData, cellDataToNodeData,
                      fillEmptyToCellArray, interpolate, interpolateAlongCurve,
                      nodeDataToBoundaryData, nodeDataToCellData,
                      tapeMeasureToCoordinates, extract2dSlice)

from .mesh import (checkMeshConsistency, convert, convertMeshioMesh,
                   convertHDF5Mesh, createMesh,
                   createParaMesh, createParaMesh2DGrid, createMeshFromHull,
                   exportFenicsHDF5Mesh, exportHDF5Mesh,
                   exportSTL, extrudeMesh, merge2Meshes, mergeMeshes,
                   readFenicsHDF5Mesh, readGmsh, readHDF5Mesh,
                   readHydrus2dMesh, readHydrus3dMesh, readSTL, readTetgen,
                   readTriangle, readMeshIO, refineHex2Tet, refineQuad2Tri,
                   toSubsurface, fromSubsurface, extractUpperSurface2dMesh)

from .polytools import createParaDomain2D  # keep for backward compatibility
from .polytools import (createCircle, createCube, createSphere, createCylinder, createFacet,
                        createSurface,
                        createLine, createParaMeshPLC,
                        createParaMeshSurface,
                        createParaMeshPLC3D,
                        createPolygon, merge,
                        createRectangle, createWorld, exportPLC, mergePLC,
                        mergePLC3D, readPLC, syscallTetgen, extrude)

from .quality import quality

#  This is neither functional nor good practice  #  why?
#  __all__ = [name for name in dir() if '_' not in name]

__all__ = [
    'appendTriangleBoundary',
    'appendTetrahedronBoundary',
    'createMesh',
    'createMeshFromHull',
    'readGmsh',
    'readTriangle',
    'readTetgen',
    'readHydrus2dMesh',
    'readHydrus3dMesh',
    'readHDF5Mesh',
    'readFenicsHDF5Mesh',
    'readSTL',
    'readMeshIO',
    'refineQuad2Tri',
    'mergeMeshes',
    'merge2Meshes',
    'createParaMesh',
    'createParaMesh2DGrid',
    'createPolygon',
    'createRectangle',
    'createWorld',
    'createCircle',
    'createLine',
    'createParaMeshPLC',
    'createParaMeshPLC3D',
    'convertHDF5Mesh',
    'exportHDF5Mesh',
    'exportFenicsHDF5Mesh',
    'extrudeMesh',
    'mergePLC',
    'readPLC',
    'writePLC',
    'exportPLC',
    'createParaDomain2D',  # keep for backward compatibility
    'quality'
]
