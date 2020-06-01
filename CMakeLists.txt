cmake_minimum_required( VERSION 3.1 )

project( EucclhydRemapProject VERSION 1.0 LANGUAGES CXX )

set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

add_subdirectory( "./cppnabla" )

add_executable( eucclhydremap EucclhydRemap.cc )
target_link_libraries( eucclhydremap 
    PRIVATE
        cppnabla
        stdc++fs
)
