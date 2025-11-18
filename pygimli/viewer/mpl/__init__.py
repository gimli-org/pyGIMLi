# -*- coding: utf-8 -*-
"""Drawing functions using Matplotlib."""


# are the following suitable for a drawing package?
from .utils import (hold,
                    wait,
                    noShow,
                    updateFig,
                    updateAxes,
                    renameDepthTicks,
                    insertUnitAtNextLastTick,
                    saveFigure,
                    saveAxes,
                    adjustWorldAxes,
                    createAnimation,
                    saveAnimation,
                    setOutputStyle,
                    setPlotStuff,
                    plotLines,
                    twin, createTwinX, createTwinY,
                    isInteractive,
                    registerShowPendingFigsAtExit)

from .boreholes import BoreHole, BoreHoles, create_legend

from .colorbar import (createColorBar,
                       createColorBarOnly,
                       updateColorBar,
                       addCoverageAlpha,
                       autolevel,
                       cmapFromName,
                       findAndMaskBestClim,
                       setCbarLevels,
                       setMappableData)
from .meshview import (CellBrowser,
                       createMeshPatches,
                       createTriangles,
                       drawBoundaryMarkers,
                       drawField,
                       drawMesh,
                       drawMeshBoundaries,
                       drawModel,
                       drawParameterConstraints,
                       drawPLC,
                       drawSelectedMeshBoundaries,
                       drawSelectedMeshBoundariesShadow,
                       drawSensors,
                       drawStreamLines,
                       drawStreams)

# do not default import
# from .overlayimage import (cacheFileName,
#                            deg2MapTile,
#                            getMapTile,
#                            mapTile2deg,
#                            underlayMap,
#                            underlayBKGMap)

from .matrixview import (drawBlockMatrix, drawSparseMatrix)

# TODO example scripts for the following and refactor is needed
# maybe ploter should named show or draw
from .dataview import (drawSensorAsMarker,  # dups to meshview??
                       showVecMatrix,
                       drawVecMatrix,
                       showValMapPatches,
                       drawValMapPatches,
                       showDataMatrix,
                       drawDataMatrix,
                       generateMatrix,
                       showDataContainerAsMatrix,
                       patchMatrix, # deprectated (Naming)
                       patchValMap, # deprectated (Naming)
                       plotDataContainerAsMatrix, # need renaming
                       plotMatrix, # deprectated (Naming)
                       plotVecMatrix,# deprectated (Naming)
                       )

# which of these do we actually need?
from .modelview import (drawModel1D,
                        #showmymatrix,  # needed ?
                        draw1DColumn, # needed or redundant ?
                        draw1dmodel,   # needed or redundant ?
                        show1dmodel,  # needed or redundant ?
                        #draw1dmodelErr,  # needed or redundant ?
                        #draw1dmodelLU,  # needed or redundant ?
                        showStitchedModels,
                        showfdemsounding)

__all__ = [
    "BoreHole", "BoreHoles", "create_legend", "addCoverageAlpha", "autolevel",
    "cmapFromName",
    "createColorBar", "createColorBarOnly", "updateColorBar",
    "findAndMaskBestClim", "setCbarLevels", "saveFigure", "saveAxes",
    "setMappableData", "drawSensorAsMarker", "generateMatrix", "patchMatrix",
    "patchValMap", "plotDataContainerAsMatrix", "plotMatrix", "plotVecMatrix",
    "CellBrowser", "createMeshPatches",
    "createTriangles", "draw1DColumn", "drawField", "drawMesh",
    "drawMeshBoundaries", "drawModel", "hold", "wait",
    "setOutputStyle", "setPlotStuff", "createAnimation", "saveAnimation",
    "drawParameterConstraints", "drawPLC", "drawSelectedMeshBoundaries",
    "drawSelectedMeshBoundariesShadow", "drawSensors", "drawStreamLines",
    "drawStreams", "insertUnitAtNextLastTick", "plotLines", "updateAxes"]
# "cacheFileName", "deg2MapTile", "getMapTile", "mapTile2deg", "underlayMap"]

