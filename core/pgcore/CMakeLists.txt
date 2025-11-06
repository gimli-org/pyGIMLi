set(TARGET_NAME pgcore)

CONFIGURE_FILE(getwheelname.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/getwheelname.py @ONLY
)
execute_process(COMMAND "${Python_EXECUTABLE}" getwheelname.py
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE PY_WHEELFILE
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    ERROR_QUIET
)

if (STATUS EQUAL 1)
    cprint(yellow "Can't determine wheel filenname. Skipping wheel build. "
                  "If you want to build wheels? python -m pip install wheel")
else()

    cprint(green "********** Prepare wheel build for: ${TARGET_NAME} ***********")

    message(STATUS "Detected Python wheel filename: ${PY_WHEELFILE}")
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

    string(REPLACE "linux" "manylinux_2_28" PY_WHEELFILE_OUTPLATFORM ${PY_WHEELFILE})
    message(STATUS "Desired Python wheel filename after audtitwheel repair: ${PY_WHEELFILE_OUTPLATFORM}")

    ## POST Build

    ### PREPARE the wheel: make whl${TARGET_NAME}
    add_custom_target(whl${TARGET_NAME}CopyLibs)
    set(WHEEL_LIBRARIES_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/libs/ )

    add_custom_command(TARGET whl${TARGET_NAME}CopyLibs
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${WHEEL_LIBRARIES_DIR}
        COMMENT "Create target dir to copy libraries: ${WHEEL_LIBRARIES_DIR}"
    )
    if (OS_VAR STREQUAL "almalinux" OR OS_VAR STREQUAL "Windows" OR OS_VAR STREQUAL "Darwin")
        ## all other libs ar copied by auditwheel, delvewheel or delocate
        set(BINDEPS_LIBS $<TARGET_FILE:_pygimli_>)
    else()
        ## copy important libs only for in source dev
        set(BINDEPS_LIBS
                #$<TARGET_FILE:${Boost_PYTHON_LIBRARY}>
                $<TARGET_FILE:gimli>
                $<TARGET_FILE:_pygimli_>)
    endif()
    foreach(file ${BINDEPS_LIBS})
        add_custom_command(TARGET whl${TARGET_NAME}CopyLibs
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND}
                            -E copy_if_different ${file} ${WHEEL_LIBRARIES_DIR}
                    #DEPENDS ${file}
                    COMMENT "Copy libraries: ${file} ${WHEEL_LIBRARIES_DIR}"
        )
    endforeach(file)

    ### PACKAGE the wheel
    add_custom_target(whl${TARGET_NAME}PackageBuild DEPENDS whl${TARGET_NAME}CopyLibs
        COMMAND "${Python_EXECUTABLE}" -m pip wheel --quiet --wheel-dir=${WHEELHOUSE} .
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Building python wheel package for ${TARGET_NAME}"
    )
    ### REPAIR the wheel
    if (OS_VAR STREQUAL "almalinux")
        add_custom_target(whl${TARGET_NAME}PackageRepair DEPENDS whl${TARGET_NAME}PackageBuild
            COMMAND "auditwheel" repair ${WHEELHOUSE}/${PY_WHEELFILE}
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMENT "Repairing python wheel package at: ${CMAKE_BINARY_DIR}"
        )
        add_custom_target(whl${TARGET_NAME} DEPENDS whl${TARGET_NAME}PackageRepair)
    elseif (OS_VAR STREQUAL "Windows")
        add_custom_target(whl${TARGET_NAME}PackageRepair DEPENDS whl${TARGET_NAME}PackageBuild
            COMMAND "delvewheel" repair ${WHEELHOUSE}/${PY_WHEELFILE} --add-path ${CMAKE_BINARY_DIR}/bin/
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMENT "Repairing python wheel package at: ${CMAKE_BINARY_DIR}"
        )
        add_custom_target(whl${TARGET_NAME} DEPENDS whl${TARGET_NAME}PackageRepair)
    elseif (OS_VAR STREQUAL "Darwin")
        add_custom_target(whl${TARGET_NAME}PackageRepair DEPENDS whl${TARGET_NAME}PackageBuild
            COMMAND "delocate-wheel" -v ${WHEELHOUSE}/${PY_WHEELFILE}
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMENT "Repairing python wheel package at: ${CMAKE_BINARY_DIR}"
        )
        add_custom_target(whl${TARGET_NAME} DEPENDS whl${TARGET_NAME}PackageRepair)
    else()
        cprint(yellow "Warning! Skipping auditwheel because of wrong linux toolchain (${OS_VAR})."
                    "For manylinux target, build on Almalinux 8 (docker image).")
        add_custom_target(whl${TARGET_NAME} DEPENDS whl${TARGET_NAME}PackageBuild)
    endif()
    ### FINAL TARGET make whlpgcore

    # copy wheel to dist folder
     add_custom_target(whl${TARGET_NAME}Dist DEPENDS whl${TARGET_NAME}
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/dist"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${WHEELHOUSE}/${PY_WHEELFILE}" "${CMAKE_BINARY_DIR}/dist/${PY_WHEELFILE}"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Copy python wheel package ${PY_WHEELFILE} to 'dist' path"
    )

    ## test the wheel: make whlTest
    add_custom_target(whl${TARGET_NAME}TestInstall DEPENDS whl${TARGET_NAME}Dist
        COMMAND "${Python3_EXECUTABLE}" -m pip uninstall -y ${TARGET_NAME}
        COMMAND "${Python3_EXECUTABLE}" -m pip install dist/${PY_WHEELFILE} --find-links ${WHEELHOUSE}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Installing ${TARGET_NAME} in virtual test environment"
    )
    add_custom_target(whl${TARGET_NAME}TestRun DEPENDS whl${TARGET_NAME}TestInstall
        COMMAND ${Python3_EXECUTABLE} -c "import pgcore; print(pgcore.versionStr())"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        VERBATIM
        COMMENT "Testing ${TARGET_NAME} installation in virtual test environment. ${Python3_EXECUTABLE} -c \'import pgcore; print(pgcore.versionStr())\'"
    )
    add_custom_target(whl${TARGET_NAME}Test DEPENDS whl${TARGET_NAME}TestRun)

    add_custom_target(whl${TARGET_NAME}Deploy DEPENDS whl${TARGET_NAME}Test
        COMMAND "${Python_EXECUTABLE}" -m twine check ${WHEELHOUSE}/${PY_WHEELFILE}
        COMMAND "${Python_EXECUTABLE}" -m twine upload --repository pypi ${WHEELHOUSE}/${PY_WHEELFILE}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Upload python wheel package ${PY_WHEELFILE} to 'pypi' repository"
    )

endif()
