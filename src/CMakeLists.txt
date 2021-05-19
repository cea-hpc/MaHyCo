set(EXAMPLE_NAME Mahyco)
project(Mahyco LANGUAGES C CXX)
include(../mahyco.common.cmake)

target_sources(${EXAMPLE_NAME} PRIVATE MahycoAnnexe.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Initialisations/InitSOD.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Initialisations/InitRider.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Initialisations/InitSEDOV.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Initialisations/InitUnitTest.cc)
target_sources(${EXAMPLE_NAME} PRIVATE PhaseRemapFinal.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Remap/PhaseRemap.cc)
target_sources(${EXAMPLE_NAME} PRIVATE PrepareRemap.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Remap/Remap-ordre-3.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Remap/UtilesRemap.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Remap/DualPhaseRemap.cc)

add_library(PerfectGas eos/perfectgas/PerfectGasEOSService.cc)
target_include_directories(PerfectGas PUBLIC .)
target_link_libraries(PerfectGas PUBLIC arcane_core)

add_library(StiffenedGas eos/stiffenedgas/StiffenedGasEOSService.cc)
target_include_directories(StiffenedGas PUBLIC .)
target_link_libraries(StiffenedGas PUBLIC arcane_core)

target_link_libraries(${EXAMPLE_NAME} PRIVATE PerfectGas StiffenedGas)

arcane_target_add_axl(PerfectGas eos/perfectgas/PerfectGasEOS)
arcane_target_add_axl(StiffenedGas eos/stiffenedgas/StiffenedGasEOS)

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

