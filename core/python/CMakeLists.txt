ADD_CUSTOM_TARGET(pggen ALL)
ADD_CUSTOM_COMMAND(
    COMMAND
        ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/_pygimli_.cache
    COMMAND
        #${CMAKE_COMMAND} -E env "PYTHONPATH=$ENV{PYTHONPATH}${PYGCCXML_PATH}" ${Python_EXECUTABLE} will work on cmake3.1
        ${Python_EXECUTABLE}
        #python
    ARGS
        ${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py
--extra-path=${PYGCCXML_PATH} --extra-includes=${Boost_INCLUDE_DIR}
--caster=${CASTER_EXECUTABLE} --clang=${CLANG_EXECUTABLE}
    TARGET
        pggen
    POST_BUILD
    #DEPENDS
        #set(ENV{PYTHONPATH} $ENV{PYTHONPATH}${PYGCCXML_PATH})
        #${CMAKE_CURRENT_SOURCE_DIR}/generate_pygimli_code.py
        #${CMAKE_CURRENT_SOURCE_DIR}/pygimli.h
        #ENV PYTHONPATH $ENV{PYTHONPATH}${PYGCCXML_PATH} APPEND
   )

#HACK!! but we need this to force cmake looking for our generated binding files
ADD_CUSTOM_TARGET(fixpygimlisrc ALL)
ADD_CUSTOM_COMMAND(
    COMMAND
        ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    TARGET
        fixpygimlisrc
    POST_BUILD
   )

ADD_CUSTOM_TARGET(pgbuild ALL)
ADD_CUSTOM_COMMAND(
    #test -z $J && ${CMAKE_MAKE_PROGRAM} -j 1 -C ${CMAKE_BINARY_DIR} pg  || echo "no 1"
    COMMAND
        if [ -z $J ]; then ${CMAKE_MAKE_PROGRAM} -j 1 -C ${CMAKE_BINARY_DIR} pg \; else ${CMAKE_MAKE_PROGRAM} -j $J -C ${CMAKE_BINARY_DIR} pg \; fi
    TARGET
        pgbuild
    POST_BUILD
   )

find_file (HAVE_DUMMY dummy.cpp PATHS ${CMAKE_CURRENT_BINARY_DIR}/generated/)
if (NOT HAVE_DUMMY)
    file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/generated/dummy.cpp "")
endif()

if(WIN32)
    ADD_CUSTOM_TARGET(pgcopydll ALL)
    ADD_CUSTOM_COMMAND(
        COMMAND
            cp ${CMAKE_BINARY_DIR}/bin/*.dll ${CMAKE_SOURCE_DIR}/pygimli/core
        TARGET
            pgcopydll
        COMMENT
            "Updating needed dll. build/bin/*dll -> gimli/pygimli/core"
    )
endif(WIN32)

add_python_module(pygimli "${CMAKE_CURRENT_BINARY_DIR}/generated/" "gimli" "${CMAKE_SOURCE_DIR}/pygimli/core/libs")
include_directories(${Python_NumPy_INCLUDE_DIRS})

## copy libgimli.so INSOURCE so we can use and develop there directly
##add_dependencies(pggen copyLibGIMLi_INSOURCE)
#add_dependencies(copy_python pggen)

add_dependencies(fixpygimlisrc pggen)
add_dependencies(pgbuild fixpygimlisrc)

if(WIN32)
    add_dependencies(pgcopydll pgbuild)
    add_dependencies(pygimli pgcopydll)
else(WIN32)
    add_dependencies(pygimli pgbuild)
endif(WIN32)

# for calling the build only

add_custom_target(pg ALL)

add_dependencies(pg _pygimli_)

ADD_CUSTOM_TARGET(pgtest
    COMMAND ${Python_EXECUTABLE} -c 'import pygimli\; pygimli.test(show=False, abort=True) '
)

#add_dependencies(check pgtest)
#add_custom_target(checkpg ALL)

# add_test(
#     NAME
#         pgtest:pygimli_tests
#     COMMAND
#         make pgtest
#     )

set(SCRIPTS pytripatch.py; meshmerge.py; meshconvert.py; pygi; datamerge.py)

foreach(program ${SCRIPTS})
    install(PROGRAMS apps/${program} DESTINATION bin)
endforeach(program)

if (WIN32)
    get_filename_component(Boost_PYTHON_ARCHIVE ${Boost_PYTHON_LIBRARY} LAST_EXT)
    # message("############################## " ${Boost_PYTHON_ARCHIVE})

    if (Boost_PYTHON_ARCHIVE MATCHES ".a")
        # message("############################## " ${Boost_PYTHON_LIBRARY})
        string(REGEX REPLACE "[.]a" "" Boost_PYTHON_DLL ${Boost_PYTHON_LIBRARY})
        # message("############################## " ${Boost_PYTHON_DLL})
        file(COPY ${Boost_PYTHON_DLL} DESTINATION ${CMAKE_BINARY_DIR}/bin)
        install(PROGRAMS ${Boost_PYTHON_DLL} DESTINATION bin)
    else ()
        #file(COPY ${Boost_PYTHON_LIBRARY} DESTINATION ${CMAKE_BINARY_DIR}/bin)
        #file(COPY ${Boost_THREAD_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
        #file(COPY ${Boost_SYSTEM_LIBRARIES} DESTINATION ${CMAKE_BINARY_DIR}/bin)
        install(PROGRAMS ${Boost_PYTHON_LIBRARY} DESTINATION bin)
    endif()
endif()
