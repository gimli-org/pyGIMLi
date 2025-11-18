#
# Find the readroc includes and library
#
# This module defines
# READPROC_INCLUDE_DIR, where to find tiff.h, etc.
# READPROC_LIBRARIES, the libraries to link against to use readproc.
# READPROC_FOUND, If false, do not try to use readproc.

include (FindPackageHandleStandardArgs)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_path (READPROC_INCLUDE_DIR proc/readproc.h)
    find_library (READPROC_LIBRARIES NAMES proc procps)
    find_package_handle_standard_args (readproc DEFAULT_MSG READPROC_LIBRARIES READPROC_INCLUDE_DIR)
endif()


# FIND_PATH(Readproc_INCLUDE_DIR proc/readproc.h
#   /usr/local/include
#   /usr/include
# )

# # On unix system, debug and release have the same name
# FIND_LIBRARY(Readproc_LIBRARIES procps
#             /usr/local/lib
#             /usr/lib
# )

# # Standard package handling
# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(readproc
#     "readproc could not be found."
#     Readproc_LIBRARIES
#     Readproc_INCLUDE_DIR)

# mark_as_advanced(
#   Readproc_INCLUDE_DIR
#   Readproc_LIBRARIES
# )