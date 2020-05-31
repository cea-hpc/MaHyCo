cmake_minimum_required( VERSION 3.1 )

project( EucclhydRemapProject VERSION 1.0 LANGUAGES CXX )

set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

list( APPEND CMAKE_PREFIX_PATH "/home/anaxi/opt/lib/cmake/Kokkos" )
find_package( Kokkos REQUIRED )

set( ENV{KOKKOS_HOME} "/home/anaxi/opt" )
add_subdirectory( ${CMAKE_SOURCE_DIR}/../libcppnabla ./libcppnabla )

add_executable( eucclhydremap EucclhydRemap.cc )
target_link_libraries(eucclhydremap 
    PRIVATE
        Kokkos::kokkos
        cppnabla
        stdc++fs
)
