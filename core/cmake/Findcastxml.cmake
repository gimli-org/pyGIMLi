#
# Find the castxml executable
#
# This module defines:
#   CASTXML_EXECUTABLE   - path to castxml
#   CASTXML_FOUND        - TRUE if castxml was found
#   CASTXML_VERSION      - version string (e.g. "0.6.5")

find_program(CASTXML_EXECUTABLE
    NAMES
        castxml
    PATHS
        ${EXTERNAL_DIR}/bin
    NO_DEFAULT_PATH
)

if (CASTXML_EXECUTABLE)
    set(CASTXML_FOUND TRUE)

    # Retrieve version string
    execute_process(
        COMMAND "${CASTXML_EXECUTABLE}" --version
        OUTPUT_VARIABLE _castxml_version_output
        ERROR_VARIABLE  _castxml_version_output
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (_castxml_version_output MATCHES "castxml version ([0-9]+\\.[0-9]+\\.?[0-9]*)")
        set(CASTXML_VERSION "${CMAKE_MATCH_1}")
    else()
        set(CASTXML_VERSION "unknown")
    endif()

    message(STATUS "Found castxml: ${CASTXML_EXECUTABLE} (version ${CASTXML_VERSION})")

    # Optional minimum version check
    if (Findcastxml_FIND_VERSION)
        if (CASTXML_VERSION VERSION_LESS Findcastxml_FIND_VERSION)
            message(WARNING "castxml version ${CASTXML_VERSION} is less than required ${Findcastxml_FIND_VERSION}")
            if (Findcastxml_FIND_VERSION_EXACT OR Findcastxml_FIND_REQUIRED)
                set(CASTXML_FOUND FALSE)
            endif()
        endif()
    endif()
else()
    message(STATUS "NOT Found castxml executable: cannot build pygimli.")
endif(CASTXML_EXECUTABLE)

mark_as_advanced(CASTXML_EXECUTABLE CASTXML_VERSION)
