cmake_minimum_required( VERSION 3.1 )

project( EucclhydRemapProject VERSION 1.0 LANGUAGES CXX )

# Disable in-source builds to prevent source tree corruption.
if( " ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}" )
  message( FATAL_ERROR "FATAL: In-source builds are not allowed.
                       You should create a separate directory for build files." )
endif()

enable_testing()

set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

add_subdirectory( "./cppnabla" )

add_executable( eucclhydremap EucclhydRemap.cc )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    target_link_libraries( eucclhydremap 
        PRIVATE
            cppnabla
            stdc++fs
    )
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    target_link_libraries( eucclhydremap 
        PRIVATE
            cppnabla
    )
else()
    message( FATAL_ERROR "Not configured to be built on Windows" )
endif()

# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/tests" )

set( TEST_NAMES "SOD_X_bi_mat_Superbee_PenteBorne_mix;SOD_X_bi_mat_MinMod_simple;SOD_X_bi_mat_Superbee_PenteBorne_simple;SOD_X_bi_mat_Superbee_simple")

foreach( TEST_NAME ${TEST_NAMES} )
    add_test( NAME ${TEST_NAME}
            COMMAND "${TEST_DIR}/launch_test.sh" $<TARGET_FILE:eucclhydremap> "${TEST_DIR}/${TEST_NAME}" )
    set_tests_properties( ${TEST_NAME} PROPERTIES
                        ENVIRONMENT OMP_NUM_THREADS=1 )
endforeach()
