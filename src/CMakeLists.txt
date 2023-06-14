set(EXAMPLE_NAME Mahyco)
project(Mahyco LANGUAGES C CXX)
include(../mahyco.common.cmake)

message(STATUS "SOURCES_SUPLEMENTAIRES=${CINETIQUE_SRC}")

target_sources(${EXAMPLE_NAME} PRIVATE MahycoAnnexe.cc)
target_sources(${EXAMPLE_NAME} PRIVATE PrepareRemap.cc)


if (IS_DIRECTORY ${CINETIQUE_SRC}) 
enable_language( Fortran )
file( GLOB SOURCES  ${CINETIQUE_SRC}/*.f90 )
add_library(cinetique SHARED ${SOURCES})
target_link_libraries(cinetique PUBLIC)
set(CINETIQUE_LIB cinetique)

add_library(Autre eos/autre/AutreEOSService.cc)
target_include_directories(Autre PUBLIC .)
target_link_libraries(Autre  PUBLIC   arcane_core ${CINETIQUE_LIB})
arcane_target_add_axl(Autre eos/autre/AutreEOS)
set(AUTRE_EOS Autre)
else()
set(AUTRE_EOS "")
endif()

add_library(PerfectGas eos/perfectgas/PerfectGasEOSService.cc)
target_include_directories(PerfectGas PUBLIC .)
target_link_libraries(PerfectGas PUBLIC arcane_core)

add_library(StiffenedGas eos/stiffenedgas/StiffenedGasEOSService.cc)
target_include_directories(StiffenedGas PUBLIC .)
target_link_libraries(StiffenedGas PUBLIC arcane_core)


add_library(RIDER casTest/RIDER/RIDERService.cc)
target_include_directories(RIDER PUBLIC .)
target_link_libraries(RIDER PUBLIC arcane_core)

add_library(SOD casTest/SOD/SODService.cc)
target_include_directories(SOD PUBLIC .)
target_link_libraries(SOD PUBLIC arcane_core)

add_library(SEDOV casTest/SEDOV/SEDOVService.cc)
target_include_directories(SEDOV PUBLIC .)
target_link_libraries(SEDOV PUBLIC arcane_core)

add_library(OTHER casTest/OTHER/OTHERService.cc)
target_include_directories(OTHER PUBLIC .)
target_link_libraries(OTHER PUBLIC arcane_core)

add_library(ONECELL casTest/ONECELL/ONECELLService.cc)
target_include_directories(ONECELL PUBLIC .)
target_link_libraries(ONECELL PUBLIC arcane_core)

add_library(POINTTRIPLE casTest/POINTTRIPLE/POINTTRIPLEService.cc)
target_include_directories(POINTTRIPLE PUBLIC .)
target_link_libraries(POINTTRIPLE PUBLIC arcane_core)

add_library(ADVECTION casTest/ADVECTION/ADVECTIONService.cc)
target_include_directories(ADVECTION PUBLIC .)
target_link_libraries(ADVECTION PUBLIC arcane_core)

add_library(CHOCBULLE casTest/CHOCBULLE/CHOCBULLEService.cc)
target_include_directories(CHOCBULLE PUBLIC .)
target_link_libraries(CHOCBULLE PUBLIC arcane_core)

add_library(INSTART casTest/INSTART/INSTARTService.cc)
target_include_directories(INSTART PUBLIC .)
target_link_libraries(INSTART PUBLIC arcane_core)

add_library(VORTEX casTest/VORTEX/VORTEXService.cc)
target_include_directories(VORTEX PUBLIC .)
target_link_libraries(VORTEX PUBLIC arcane_core)

add_library(Remap Remap/RemapADIService.cc Remap/RemapALEService.cc Remap/Remap-ordre-3.cc Remap/UtilesRemap.cc Remap/RemapADIFinal.cc Remap/UtilesRemapALE.cc Remap/DualPhaseRemap.cc)
target_include_directories(Remap PUBLIC .)
target_link_libraries(Remap PUBLIC arcane_core)

#include_directories(${MyLib_INCLUDE_DIRS})

target_link_libraries(${EXAMPLE_NAME} PRIVATE PerfectGas StiffenedGas ${AUTRE_EOS} RIDER POINTTRIPLE ADVECTION SOD SEDOV OTHER CHOCBULLE INSTART ONECELL VORTEX Remap)

arcane_target_add_axl(RIDER casTest/RIDER/RIDER)
arcane_target_add_axl(SOD casTest/SOD/SOD)
arcane_target_add_axl(SEDOV casTest/SEDOV/SEDOV)
arcane_target_add_axl(OTHER casTest/OTHER/OTHER)
arcane_target_add_axl(ONECELL casTest/ONECELL/ONECELL)
arcane_target_add_axl(POINTTRIPLE casTest/POINTTRIPLE/POINTTRIPLE)
arcane_target_add_axl(ADVECTION casTest/ADVECTION/ADVECTION)
arcane_target_add_axl(CHOCBULLE casTest/CHOCBULLE/CHOCBULLE)
arcane_target_add_axl(INSTART casTest/INSTART/INSTART)
arcane_target_add_axl(VORTEX casTest/VORTEX/VORTEX)

arcane_target_add_axl(PerfectGas eos/perfectgas/PerfectGasEOS)
arcane_target_add_axl(StiffenedGas eos/stiffenedgas/StiffenedGasEOS)

arcane_target_add_axl(Remap Remap/RemapADI)
arcane_target_add_axl(Remap Remap/RemapALE)

# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E ${CMAKE_CURRENT_BINARY_DIR}/Mahyco)


# Pour les test NONREG
# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )


set( TEST_NAMES_PARA_4
        "CAS_2D_Lag_BiSodCaseX"
        "CAS_2D_Lag_BiSodCaseY"
		"CAS_2D_Lag_SodCaseY"
		"CAS_Lag_BiSodCaseZ"
        "CAS_Lag_SodCaseX"
        "CAS_BiSodCaseY"
		"CAS_BiSodCaseZ"
		"CAS_SodCaseY"
        "CAS_SodCaseZ"
        "CAS_SodCaseX_ordre_3"
        "CAS_BiSodCaseX_ordre_3"
        "CAS_SodCaseY_ordre_3"
        "CAS_SodCaseZ_ordre_3"
		"CAS_BiSodCaseY_ordre_3"
        "CAS_BiSodCaseZ_ordre_3"
		"CAS_BiSodCaseX_pente_borne_minmod_debar_1"
		"CAS_BiSodCaseX_pente_borne_minmod_debar_2"
		"CAS_BiSodCaseX_pente_borne_superbee_debar_2"
		"CAS_BiSodCaseY_pente_borne_minmod_debar_2"
		"CAS_BiSodCaseY_pente_borne_superbee_debar_2"
		)
set( TEST_NAMES_PARA_8		
		"CAS_MonoRiderDeformationTimeReverse"
		"CAS_MonoRiderVortexTimeReverse"
		"CAS_MonoRiderVortex"
		"CAS_RiderDeformation"
		"CAS_RiderDeformationTimeReverse"
		"CAS_RiderRotation"
		"CAS_RiderTx"
		"CAS_RiderTy"
		"CAS_RiderVortex"
		"CAS_RiderVortexTimeReverse"
		"CAS_MonoRiderTy_pb_simple_superbee"
		"CAS_MonoRiderTy_pb_mixte_superbee_vanleer"
		"CAS_Sedov"
		)
set( TEST_NAMES
		"CAS_Lag_BiSodCaseX"
        "CAS_2D_Lag_SodCaseZ"
        "CAS_Lag_BiSodCaseY"
		"CAS_Lag_SodCaseY"
        "CAS_Lag_SodCaseZ"
		"CAS_BiSodCaseX"
		"CAS_MonoRiderDeformation"
		"CAS_MonoRiderRotation"
		"CAS_MonoRiderT45"
		"CAS_MonoRiderTy"
		"CAS_2D_Lag_SodCaseX_newton"
		"CAS_2D_Sedov"
		"CAS_2D_ALE_Sedov_euler"
		"CAS_2D_ALE_BiSedov_euler"
		)

set( TEST_NAMES_PR
        "CAS_2D_Lag_SodCaseX"
        "CAS_SodCaseX"
		"CAS_MonoRiderTx"
		"CAS_RiderT45"
		)
		
set( PARA4 "para_4")
set( PARA8 "para_8")
set( SEQ "seq")
set( SEQPR "seq_pr")

# Pour définir le launcher MPI
set(MPI_LAUNCHER mpiexec CACHE STRING "Définit le launcher MPI" )

# Pour définir la base du répertoire temporaire d'exécution des tests de non-regression
set(BASETMPDIR_NONREG "/tmp" CACHE PATH "Définit la base du répertoire temporaire d'exécution des tests de non-regression" )

# protection-reprise
foreach( TEST_NAME ${TEST_NAMES_PR} )
 add_test(NAME "${SEQPR}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${SEQPR} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach()

# sur 4 coeurs
foreach( TEST_NAME ${TEST_NAMES_PARA_4} )
 add_test(NAME "${PARA4}_${TEST_NAME}" COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${PARA4} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach()

# sur 8 coeurs
foreach( TEST_NAME ${TEST_NAMES_PARA_8} )
 add_test(NAME "${PARA8}_${TEST_NAME}" COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${PARA8} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach()

# sequentiel
foreach( TEST_NAME ${TEST_NAMES} )
 add_test(NAME "${SEQ}_${TEST_NAME}" COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${SEQ} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach()

