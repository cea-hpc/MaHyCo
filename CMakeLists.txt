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
add_subdirectory( "./eucclhyd_remap" )
add_subdirectory( "./tests" )

add_executable( mahyco )
target_sources( mahyco 
    PRIVATE
        main.cc
)
target_link_libraries( mahyco PRIVATE eucclhydremap )