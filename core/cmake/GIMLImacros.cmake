################################################################################
# Macro definitions used by GIMLI's cmake build
################################################################################
function(cprint)
    list(GET ARGV 0 Color)
    list(REMOVE_AT ARGV 0)
    execute_process(COMMAND 
    cmake -E cmake_echo_color --${Color} ${ARGV})


#   if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
#     list(REMOVE_AT ARGV 0)
#     _message(${MessageType} "${BoldRed}${ARGV}${ColourReset}")
#   elseif(MessageType STREQUAL WARNING)
#     list(REMOVE_AT ARGV 0)
#     _message(${MessageType} "${BoldYellow}${ARGV}${ColourReset}")
#   elseif(MessageType STREQUAL AUTHOR_WARNING)
#     list(REMOVE_AT ARGV 0)
#     _message(${MessageType} "${BoldCyan}${ARGV}${ColourReset}")
#   elseif(MessageType STREQUAL STATUS)
#     list(REMOVE_AT ARGV 0)
#     _message(${MessageType} "${ARGV}")
#     _message(STATUS "${ARGV}")
#     #_message(${MessageType} "${Green}${ARGV}${ColourReset}")
#   else()
#     _message("${ARGV}")
#   endif()
endfunction()


macro(add_python_module PYTHON_MODULE_NAME SOURCE_DIR EXTRA_LIBS OUTDIR)

    set(PYTHON_TARGET_NAME "_${PYTHON_MODULE_NAME}_")

    file(GLOB ${PYTHON_MODULE_NAME}_SOURCE_FILES ${SOURCE_DIR}/*.cpp)
    list(SORT ${PYTHON_MODULE_NAME}_SOURCE_FILES)

    set_source_files_properties(${PYTHON_MODULE_NAME}_SOURCE_FILES
                                 PROPERTIES GENERATED TRUE)

    include_directories(BEFORE ${SOURCE_DIR})
    include_directories(${Python_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated/)
    
    add_definitions(-DPYGIMLI)
    add_definitions(-DBOOST_PYTHON_NO_PY_SIGNATURES)
    add_definitions(-DBOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY)
    
    add_library(${PYTHON_TARGET_NAME} MODULE ${${PYTHON_MODULE_NAME}_SOURCE_FILES})
    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES PREFIX "")
        
    #TODO check!! (python3-config --extension-suffix)

    if (APPLE)
        # set(GIMLI_LIBRARY "${CMAKE_BINARY_DIR}/${LIBRARY_INSTALL_DIR}/libgimli.dylib")
        target_link_libraries(${PYTHON_TARGET_NAME} "-bundle -undefined dynamic_lookup")
    elseif (WIN32)
        # set(GIMLI_LIBRARY "${CMAKE_BINARY_DIR}/bin/libgimli.dll")
        set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES SUFFIX ".pyd")
    else()
        # set(GIMLI_LIBRARY "${CMAKE_BINARY_DIR}/${LIBRARY_INSTALL_DIR}/libgimli.so")
    endif()
    
    # target_link_libraries(${PYTHON_TARGET_NAME} ${GIMLI_LIBRARY}) 
    # target_link_libraries(${PYTHON_TARGET_NAME} $<TARGET_FILE:gimli>) 
    target_link_libraries(${PYTHON_TARGET_NAME} gimli) 
    target_link_libraries(${PYTHON_TARGET_NAME} ${Boost_PYTHON_LIBRARY})
    #target_link_libraries(${PYTHON_TARGET_NAME} ${Python_LIBRARIES})

    if (Python_Development.Module_FOUND)
        #target_link_libraries(${PYTHON_TARGET_NAME} PRIVATE Python::Module)
        target_link_libraries(${PYTHON_TARGET_NAME} Python::Module)
    else()
        target_link_libraries(${PYTHON_TARGET_NAME} ${Python_LIBRARIES})
    endif()


    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY_DEBUG ${OUTDIR})
    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${OUTDIR})
    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${OUTDIR})
    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
                        LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${OUTDIR})
    
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_CLANGXX)
	    set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES
                            COMPILE_FLAGS "-fvisibility=hidden -Wno-unused-value -Wno-infinite-recursion"
                                )
        
        if (WIN32 AND ADDRESSMODEL EQUAL "64")
            set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES 
                                DEFINE_SYMBOL "MS_WIN64")
        endif()
    endif()
    
    
    #--copy pattern files to build folder--
    ## needed?
    # set(PYTHON_IN_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
    # set(PYTHON_OUT_PATH "${CMAKE_BINARY_DIR}/package")

    # file(GLOB_RECURSE PYTHON_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    #                 "${PYTHON_MODULE_NAME}/*.py"
    #                 "${PYTHON_MODULE_NAME}/*.png"
    #                 "${PYTHON_MODULE_NAME}/*.xrc"
    #                 "${PYTHON_MODULE_NAME}/*.fbp")

    

    # foreach(file ${PYTHON_FILES})

    #     #message ("${PYTHON_IN_PATH}/${file} ${PYTHON_OUT_PATH}/${file}")
    #     add_custom_command(
    #         COMMAND
    #             cmake -E copy_if_different
    #             ${PYTHON_IN_PATH}/${file}
    #             ${PYTHON_OUT_PATH}/${file}
    #         DEPENDS "${PYTHON_IN_PATH}/${file}"
    #         TARGET
    #             copy_python
    #         VERBATIM
    #         COMMENT
    #             "Updating python file: ${file}"
    #     )
    # endforeach(file)


#----install-----------------------
#     foreach(file ${PYTHON_FILES})
#         get_filename_component( path_name "${file}" PATH )
#         #file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pygimli/${file} DESTINATION ${path_name})
#         install( FILES "${file}" DESTINATION ${path_name} )
#     endforeach(file)
#
#     install(TARGETS ${PYTHON_TARGET_NAME} LIBRARY DESTINATION "${PYTHON_MODULE_NAME}/")
endmacro()

function(find_python_module module)
    string(TOUPPER ${module} module_upper)
    if(NOT PY_${module_upper})
        if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
            set(${module}_FIND_REQUIRED TRUE)
        endif()
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND "${Python_EXECUTABLE}" "-c"
            "import re, ${module}; print(re.compile('\\__init__.py.*').sub('',${module}.__file__))"
            RESULT_VARIABLE _${module}_status
            OUTPUT_VARIABLE _${module}_location
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        if(NOT _${module}_status)
            set(PY_${module_upper} ${_${module}_location} CACHE STRING
                "Location of Python module ${module}")
            message(STATUS "Find python module ${module}")
            message(STATUS "Result: ${_${module}_status}")
            message(STATUS "Output: ${_${module}_location}")
        else(NOT _${module}_status)
            message(STATUS "Find python module ${module} fails:")
            message(STATUS "Result: ${_${module}_status}")
            message(STATUS "Output: ${_${module}_location}")
            message(STATUS "Maybe you can provide the location by stetting PY_${module_upper}")
        endif(NOT _${module}_status)

        find_package_handle_standard_args(${module}
                                      FOUND_VAR ${module}_FOUND
                                      REQUIRED_VARS PY_${module_upper}
                                      )
    else(NOT PY_${module_upper})
        set(${module}_FOUND TRUE)
    endif(NOT PY_${module_upper})

    if (${module}_FOUND)
        set( ${module}_FOUND ${${module}_FOUND} CACHE INTERNAL ${module}_FOUND)
        set( ${module}_LOC ${_${module}_location} CACHE INTERNAL ${module}_LOC)
    endif()
endfunction(find_python_module)


macro(findBuildTools)
    message(STATUS "checking for some build tools ...")
    #unzip try cmake -E tar
    #find_package(Tar REQUIRED)  ${CMAKE_COMMAND} -E tar "cfvz"
    find_program(PATCH_TOOL NAMES patch  REQUIRED)
    find_program(SED_TOOL NAMES sed REQUIRED)
    find_package(Wget REQUIRED)
    find_package(Git REQUIRED)
endmacro(findBuildTools)


macro(find_or_build_package package get_package)

    set (extraMacroArgs ${ARGN})

    # Did we get any optional args?
    list(LENGTH extraMacroArgs numExtraArgs)
    if (${numExtraArgs} GREATER 0)
        list(GET extraMacroArgs 0 optionalArg)
        set(foceLocal True)
    else()
        set(foceLocal False)
    endif()

    string(TOUPPER ${package} upper_package)
    find_or_build_package_check(${package} ${get_package} ${upper_package}_FOUND ${foceLocal})
endmacro()

macro(build_package package get_package)
    findBuildTools()

    message(STATUS "building ${package} from foreign sources into ${THIRDPARTY_DIR}" )

    file(MAKE_DIRECTORY ${THIRDPARTY_DIR})

    if (J)
        set(ENV{PARALLEL_BUILD} ${J})
    endif()

    if (NOT get_package)
        set(get_package ${lower_package})
    endif()

    execute_process(
        COMMAND
            bash ${PROJECT_SOURCE_DIR}/core/scripts/buildThirdParty.sh ${get_package}
        WORKING_DIRECTORY
            ${THIRDPARTY_DIR}
    )
endmacro()

macro(find_or_build_package_check 
        package 
        get_package 
        checkVar 
        forceLocal
        )

    message(STATUS "** Find or build ${package} at: ${checkVar} force: ${forceLocal}")
    find_package(${package})
    message(STATUS "Found: ${${package}_FOUND}")

    string(TOUPPER ${package} upper_package)
    string(TOLOWER ${package} lower_package)

    set (FORCE_LOCAL_REBUILD 0)

    message(STATUS "Local build ${package} forced: ${forceLocal}")

    if ($ENV{CLEAN})
        if(${forceLocal} OR ${package}_LOCAL)
            set(FORCE_LOCAL_REBUILD 1)
            set(ENV{CLEAN} 1)
            message(STATUS "Rebuild forced for: ${package}")
        endif()
    endif()

    if (NOT ${checkVar} OR ${FORCE_LOCAL_REBUILD})

        build_package(${package} ${get_package})

        message(STATUS "checking again for ${package} ...")
		find_package(${package})
        message(STATUS "Found: ${${package}_FOUND}")

        set(${package}_LOCAL 1 CACHE INTERNAL "this package was build local")
    else()
        message(STATUS "** Find or build ${package} done.")
    endif()
endmacro()


