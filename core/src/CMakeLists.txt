set(libgimli_TARGET_NAME gimli)

link_directories(${libgimli_BINARY_DIR})

file (GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp" "bert/*.cpp")
file (GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.h" "bert/*.h")

if (NOT MSVC)
	list(REMOVE_ITEM SOURCE_FILES wingetopt.cpp)
endif(NOT MSVC)

###########################################################################
#                       SHARED LIBRARY SET UP
###########################################################################
add_library(${libgimli_TARGET_NAME} SHARED ${SOURCE_FILES} ${HEADER_FILES})

include (GenerateExportHeader)

generate_export_header( ${libgimli_TARGET_NAME}
             BASE_NAME ${libgimli_TARGET_NAME}
             EXPORT_MACRO_NAME DLLEXPORT
             EXPORT_FILE_NAME ${libgimli_TARGET_NAME}_export.h
             STATIC_DEFINE LIBGIMLI_BUILT_AS_STATIC
)

if (TRIANGLE_FOUND)
    include_directories(${Triangle_INCLUDE_DIR})
    target_link_libraries(${libgimli_TARGET_NAME} ${Triangle_LIBRARIES})
endif (TRIANGLE_FOUND)

if (CHOLMOD_FOUND)
    include_directories(${CHOLMOD_INCLUDE_DIRS})
    target_link_libraries(${libgimli_TARGET_NAME} ${CHOLMOD_LIBRARIES})
endif (CHOLMOD_FOUND)

if (UMFPACK_FOUND)
    include_directories(${UMFPACK_INCLUDE_DIRS})
    target_link_libraries(${libgimli_TARGET_NAME} ${UMFPACK_LIBRARIES})
endif (UMFPACK_FOUND)

if (PYTHON_FOUND)
    include_directories(${Python_INCLUDE_DIRS})

    if (NOT APPLE)
        # needed since the later doesnt seems to work
        target_link_libraries(${libgimli_TARGET_NAME} ${Python_LIBRARIES})

        if (Python_Development.Module_FOUND)
            target_link_libraries(${libgimli_TARGET_NAME} Python::Module)
        endif()
    endif()

endif (PYTHON_FOUND)

if (Boost_INCLUDE_DIR)
    include_directories(${Boost_INCLUDE_DIR})
endif()

if (CONDA_BUILD)
    include_directories(${CONDA_PREFIX})
endif()

if (OpenBLAS_INCLUDE_DIR)
    include_directories(${OpenBLAS_INCLUDE_DIR})
    target_link_libraries(${libgimli_TARGET_NAME} ${BLAS_LIBRARIES})

    # if (NOT APPLE)
    #     target_link_libraries(${libgimli_TARGET_NAME} openblas
    # endif()
endif()

if (USE_BOOST_THREAD)
    target_link_libraries(${libgimli_TARGET_NAME} ${Boost_THREAD_LIBRARIES})
    target_link_libraries(${libgimli_TARGET_NAME} ${Boost_SYSTEM_LIBRARIES})
    message(STATUS "+++++++++++++++++++++++++++++++++ ${Boost_THREAD_LIBRARIES}")

    if (Boost_THREAD_FOUND)
        target_link_libraries(${libgimli_TARGET_NAME} ${Boost_THREAD_LIBRARIES})
    endif (Boost_THREAD_FOUND)

    if (Boost_SYSTEM_FOUND)
        target_link_libraries(${libgimli_TARGET_NAME} ${Boost_SYSTEM_LIBRARIES})
    endif (Boost_SYSTEM_FOUND)
endif(USE_BOOST_THREAD)

## copy libgimli.so INSOURCE so we can use and develop there directly
set(INSOURCE_LIBRARIES_DIR ${CMAKE_SOURCE_DIR}/pygimli/core/libs/ )
add_custom_command(TARGET ${libgimli_TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${INSOURCE_LIBRARIES_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${libgimli_TARGET_NAME}> ${INSOURCE_LIBRARIES_DIR}
    COMMENT "copy $<TARGET_FILE:gimli> to ${INSOURCE_LIBRARIES_DIR}"
)

# add_custom_target(copyLibGIMLi_INSOURCE DEPENDS gimli
# COMMAND
#     cmake -E copy_if_different $<TARGET_FILE:gimli> ${OUTDIR}
# COMMENT "copy $<TARGET_FILE:gimli> to ${OUTDIR}"
# )

if (OpenMP_CXX_FOUND)
    target_link_libraries(${libgimli_TARGET_NAME} OpenMP::OpenMP_CXX)
endif()

if (APPLE)
    target_link_libraries(${libgimli_TARGET_NAME} "-undefined dynamic_lookup")
endif()

if (WIN32)
    target_link_libraries(${libgimli_TARGET_NAME} psapi)

    #if (USE_BOOST_THREAD)
        #file(COPY ${Boost_THREAD_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
        #file(COPY ${Boost_SYSTEM_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)

        #install(PROGRAMS ${Boost_THREAD_LIBRARIES} DESTINATION bin)
        #install(PROGRAMS ${Boost_SYSTEM_LIBRARIES} DESTINATION bin)
    #endif(USE_BOOST_THREAD)

	if(${BLAS_LIBRARIES})
		file(COPY ${BLAS_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
		if (LAPACK_LIBRARIES)
			file(COPY ${LAPACK_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
		endif()
		install(PROGRAMS ${BLAS_LIBRARIES} DESTINATION bin)
		install(PROGRAMS ${LAPACK_LIBRARIES} DESTINATION bin)
	endif(${BLAS_LIBRARIES})

elseif (UNIX AND NOT APPLE)
	if (READPROC_FOUND)
        #message(STATUS " add:" ${READPROC_LIBRARIES})
        target_link_libraries(${libgimli_TARGET_NAME} ${CMAKE_DL_LIBS})
		include_directories(${READPROC_INCLUDE_DIR})
		target_link_libraries(${libgimli_TARGET_NAME} ${READPROC_LIBRARIES})
	endif (READPROC_FOUND)

	install(TARGETS ${libgimli_TARGET_NAME} LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})

endif(WIN32)

# install(FILES ${HEADER_FILES} ${PROJECT_BINARY_DIR}/${libgimli_TARGET_NAME}_export.h
#         DESTINATION ${INCLUDE_INSTALL_DIR})

