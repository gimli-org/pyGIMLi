set(TARGET_NAME pygimli)

CONFIGURE_FILE(getwheelname.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/getwheelname.py @ONLY
)
execute_process(COMMAND "${Python_EXECUTABLE}" getwheelname.py
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE PY_WHEELFILE
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    ERROR_QUIET
)

set(WHEELFILE ${PY_WHEELFILE})

if (STATUS EQUAL 1)
    cprint(yellow "Can't determine wheel filenname. Skipping wheel build. "
                  "If you want to build wheels? python -m pip install wheel")
else()

    cprint(green "********** Prepare wheel build for: ${TARGET_NAME} ***********")

    function(get_versions versionstring major minor patch)
        string(REGEX REPLACE "^([vV])([0-9]*)([.][0-9]*[.][0-9]*-?.*)$" "\\2" numbers ${versionstring} )
        set(major ${numbers} PARENT_SCOPE)
        string(REGEX REPLACE "^([vV][0-9]*[.])([0-9]*)([.][0-9]*-?.*)$" "\\2" numbers ${versionstring} )
        set(minor ${numbers} PARENT_SCOPE)
        string(REGEX REPLACE "^([vV][0-9]*[.][0-9]*[.])([0-9]*)(-?.*)$" "\\2" numbers ${versionstring} )
        set(patch ${numbers} PARENT_SCOPE)
    endfunction()

    get_versions(${GIMLI_GIT_TAG} major minor patch)
    set(PGCORE_VERSION_WHL ${major}.${minor})

    message(STATUS "Detected Python wheel filename: ${WHEELFILE}")
    CONFIGURE_FILE(setup.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/setup.py @ONLY
    )

    message(STATUS "Copying python files to ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/${TARGET_NAME}")
    file(GLOB_RECURSE PYGIMLI_FILES RELATIVE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
    foreach (filename ${PYGIMLI_FILES})
        #message(STATUS "Copying ${filename} to ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/${TARGET_NAME}/${target}")
        configure_file("${filename}" "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/${filename}" COPYONLY)
    endforeach (filename)

    add_custom_target(whl${TARGET_NAME}PackageBuild DEPENDS whlpgcore
        COMMAND "${Python_EXECUTABLE}" -m pip wheel --quiet . --wheel-dir=${WHEELHOUSE} --find-links ${WHEELHOUSE}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Building python wheel package (${WHEELFILE}) for ${TARGET_NAME}"
    )

    add_custom_target(whl${TARGET_NAME} DEPENDS whl${TARGET_NAME}PackageBuild)

    ## test the wheel: make whlTest
    # update the environment with VIRTUAL_ENV variable (mimic the activate script) and search again
    add_custom_target(whl${TARGET_NAME}TestInstall DEPENDS whl${TARGET_NAME} DEPENDS whlpgcoreTestRun
        COMMAND "${Python3_EXECUTABLE}" -m pip uninstall -y ${TARGET_NAME}
        COMMAND "${Python3_EXECUTABLE}" -m pip install ${WHEELHOUSE}/${PY_WHEELFILE} --quiet --find-links ${WHEELHOUSE}
        COMMENT "Installing ${TARGET_NAME} ((${WHEELFILE}) in virtual test environment"
    )
    add_custom_target(whl${TARGET_NAME}TestRun DEPENDS whl${TARGET_NAME}TestInstall
        COMMAND ${Python3_EXECUTABLE} -c "import pygimli as pg; pg.version(cache=False)"
        VERBATIM
        COMMENT "Testing ${TARGET_NAME} installation (${WHEELFILE}) in virtual test environment. ${Python3_EXECUTABLE} -c \'import pygimli as pg; pg.version()\'"
    )
    add_custom_target(whl${TARGET_NAME}Test DEPENDS whl${TARGET_NAME}TestRun)

    add_custom_target(whlTest DEPENDS whl${TARGET_NAME}TestRun)
    add_custom_target(whl DEPENDS whl${TARGET_NAME})

    ## deploy the wheel: make whlpygimliDeploy
    add_custom_target(whl${TARGET_NAME}Deploy DEPENDS whl${TARGET_NAME}Test
        COMMAND "${Python_EXECUTABLE}" -m twine check ${WHEELHOUSE}/${WHEELFILE}
        COMMAND "${Python_EXECUTABLE}" -m twine upload --repository pypi ${WHEELHOUSE}/${WHEELFILE}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Upload python wheel package ${WHEELFILE} to 'pypi' repository"
    )
endif()