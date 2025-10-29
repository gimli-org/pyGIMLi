GIMLab - GIMLi MatLab toolboxes

This matlab directory provides functions around GIMLi to read, write and plot different stuff like matrices, meshes etc. However, they are very old (formerly called DCmatlab) out of a time when we switched to Python and therefore they are not developed or tested anymore. Please do no expect anything like the Python module pyGIMLi that brings not only the Python bindings to the C++ library GIMLi, but a whole universe around it. 
Therefore we strongly recommend using Python instead.

The functions are organized in several toolboxes:

base    - base functions like matrix input/output, string representation, numerical base
eqsolve - equation solvers for different kinds, mostly for least squares problems
mesh    - functions for mesh input/output and manipulation
plot    - plotting functions from 1d to 3d
poly    - import, export and manipulation of PLC' (poly file) for mesh input
trans   - transformation functions