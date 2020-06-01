set( LIBRARY_NAME "cppnabla" )

list( APPEND CMAKE_PREFIX_PATH "/home/anaxi/opt/lib/cmake/Kokkos" )
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