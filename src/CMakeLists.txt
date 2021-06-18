set(EXAMPLE_NAME Mahyco)
project(Mahyco LANGUAGES C CXX)
include(../mahyco.common.cmake)

target_sources(${EXAMPLE_NAME} PRIVATE MahycoAnnexe.cc)
target_sources(${EXAMPLE_NAME} PRIVATE PhaseRemapFinal.cc)
target_sources(${EXAMPLE_NAME} PRIVATE PrepareRemap.cc)

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


add_library(Remap Remap/RemapADIService.cc Remap/RemapALEService.cc Remap/Remap-ordre-3.cc Remap/UtilesRemap.cc Remap/UtilesRemapALE.cc Remap/DualPhaseRemap.cc)
target_include_directories(Remap PUBLIC .)
target_link_libraries(Remap PUBLIC arcane_core)



target_link_libraries(${EXAMPLE_NAME} PRIVATE PerfectGas StiffenedGas RIDER SOD SEDOV OTHER Remap)

arcane_target_add_axl(RIDER casTest/RIDER/RIDER)
arcane_target_add_axl(SOD casTest/SOD/SOD)
arcane_target_add_axl(SEDOV casTest/SEDOV/SEDOV)
arcane_target_add_axl(OTHER casTest/OTHER/OTHER)

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

# protection-reprise
foreach( TEST_NAME ${TEST_NAMES_PR} )
 add_test(NAME "${TEST_NAME}" COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${TEST_DIR}/../build/src/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${SEQPR})
endforeach()

# sur 4 coeurs
foreach( TEST_NAME ${TEST_NAMES_PARA_4} )
 add_test(NAME ${TEST_NAME} COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${TEST_DIR}/../build/src/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${PARA4} )
endforeach()

# sur 8 coeurs
foreach( TEST_NAME ${TEST_NAMES_PARA_8} )
 add_test(NAME ${TEST_NAME} COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${TEST_DIR}/../build/src/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${PARA8} )
endforeach()

# sequentiel
foreach( TEST_NAME ${TEST_NAMES} )
 add_test(NAME "${TEST_NAME}" COMMAND  ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${TEST_DIR}/../build/src/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${SEQ})
endforeach()

