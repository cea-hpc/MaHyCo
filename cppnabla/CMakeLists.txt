set( LIBRARY_NAME "cppnabla" )

set( KOKKOS_INSTALL_PATH "Undef" CACHE STRING "Path to the Kokkos installation directory" )

message( STATUS "KOKKOS_INSTALL_PATH is ${KOKKOS_INSTALL_PATH}" )
if( ${KOKKOS_INSTALL_PATH} STREQUAL "Undef" )
    message( FATAL_ERROR "Please specify the KOKKOS installation directory (-DKOKKOS_INSTALL_PATH)" )
endif()

list( APPEND CMAKE_PREFIX_PATH "${KOKKOS_INSTALL_PATH}" )
find_package( Kokkos REQUIRED )

if( ${Kokkos_FOUND} )
    message( STATUS "Kokkos found" )
endif()

add_library( ${LIBRARY_NAME} )

target_sources( ${LIBRARY_NAME} PRIVATE
    "types/LinearAlgebraFunctions.cc"
    "types/Matrix.cc"
    "types/MultiArray.cc"
    "mesh/CartesianMesh2D.cc"
    "mesh/CartesianMesh2DGenerator.cc"
    "mesh/FileWriter.cc"
    "mesh/VtkFileWriter2D.cc"
    "mesh/PvdFileWriter2D.cc"
    "utils/Utils.cc"
    "utils/Timer.cc"
)

target_include_directories( ${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_link_libraries( ${LIBRARY_NAME}
PUBLIC
    Kokkos::kokkos
)