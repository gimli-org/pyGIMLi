file (GLOB unitTest_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.h)

#PURIFYCOMMAND:FILEPATH=c:/Progra~1/Rational/common/purify.exe
#MEMORYCHECK_COMMAND:FILEPATH=/home/kitware/local/bin/valgrind
#set(WITH_MEMCHECK TRUE)
#set(WITH_COVERAGE TRUE)

set(TARGET_NAME gimliUnitTest)

add_executable(${TARGET_NAME} unitTest.cpp ${unitTest_HEADERS})

if (CPPUNIT_FOUND)

add_definitions(-DHAVE_UNITTEST)
    
    include_directories(${CPPUNIT_INCLUDE_DIR})
    
    target_link_libraries(${TARGET_NAME} gimli)
    target_link_libraries(${TARGET_NAME} ${CPPUNIT_LIBRARIES})
    
    # target_link_libraries(${TARGET_NAME} ${Python_LIBRARIES})
        
    # if (Python_Development.Module_FOUND)
    #     target_link_libraries(${TARGET_NAME} Python::Module)
    # else()
    #     target_link_libraries(${TARGET_NAME} ${Python_LIBRARIES})
    # endif()

    add_dependencies(${TARGET_NAME} gimli)
    
    ADD_CUSTOM_TARGET(gtest DEPENDS ${TARGET_NAME}
        COMMAND ${CMAKE_BINARY_DIR}/bin/${TARGET_NAME}
    )
    
    #add_dependencies(check gtest)

    # add_test(
    #     NAME 
    #         gtest:gimli_unit_tests
    #     COMMAND
    #         #cmake "--build . --target gtest"
    #         make gtest
    #     )

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/unittest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
   
    INSTALL_FILES(${CMAKE_CURRENT_SOURCE_DIR} unittest.sh)

endif()
