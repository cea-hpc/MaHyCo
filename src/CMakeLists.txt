set(EXAMPLE_NAME Mahyco)
project(Mahyco LANGUAGES C CXX)
include(../mahyco.common.cmake)

target_sources(${EXAMPLE_NAME} PRIVATE MahycoAnnexe.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Initialisations/InitSOD.cc)
target_sources(${EXAMPLE_NAME} PRIVATE Initialisations/InitRider.cc)
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
set(HOME /home/perlatj/MAHYCO_ARCANE/Mahyco)
# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )

set( TEST_NAMES "CAS_2D_Lag_BiSodCaseX"
     		"CAS_2D_Lag_BiSodCaseY"
		"CAS_2D_Lag_BiSodCaseZ"
     		"CAS_2D_Lag_SodCaseX"
		"CAS_2D_Lag_SodCaseY"
     		"CAS_2D_Lag_SodCaseZ"
		"CAS_Lag_BiSodCaseX"
     		"CAS_Lag_BiSodCaseY"
		"CAS_Lag_BiSodCaseZ"
     		"CAS_Lag_SodCaseX"
		"CAS_Lag_SodCaseY"
     		"CAS_Lag_SodCaseZ"
		"CAS_MonoRiderDeformation"
		"CAS_MonoRiderDeformationTimeReverse"
		"CAS_MonoRiderRotation"
		"CAS_MonoRiderT45"
		"CAS_MonoRiderTx"
		"CAS_MonoRiderTy"
		"CAS_MonoRiderVortex"
		"CAS_MonoRiderVortexTimeReverse"
		"CAS_RiderDeformation"
		"CAS_RiderDeformationTimeReverse"
		"CAS_RiderRotation"
		"CAS_RiderT45"
		"CAS_RiderTx"
		"CAS_RiderTy"
		"CAS_RiderVortex"
		"CAS_RiderVortexTimeReverse")

foreach( TEST_NAME ${TEST_NAMES} )
 add_test(NAME ${TEST_NAME} COMMAND  ${HOME}/NONREGRESSSION/launch_test.sh ${HOME}/src/Mahyco ${HOME}/NONREGRESSSION/${TEST_NAME})


#    add_test( NAME ${TEST_NAME}
#              COMMAND "${TEST_DIR}/launch_test.sh" $<TARGET_FILE:mahyco> "${TEST_DIR}/${TEST_NAME}" )
#    set_tests_properties( ${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1 )
endforeach()




if(MPIEXEC_EXECUTABLE)
#  add_test(NAME mahyco_NR1_p4
#    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 Mahyco ../../NONREGRESSSION/CAS_2D_Lag_BiSodCaseX/Donnees.arc
#    )
  # Teste le lancement avec 'arcane_run'
#  add_test(NAME mahyco_NR1_p4_run
#    COMMAND ${RUN_COMMAND} -n 4 ../../NONREGRESSSION/CAS_2D_Lag_BiSodCaseX/Donnees.arc
#    )
endif()
