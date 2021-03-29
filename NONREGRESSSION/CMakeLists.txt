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
    add_test( NAME ${TEST_NAME}
              COMMAND "${TEST_DIR}/launch_test.sh" $<TARGET_FILE:mahyco> "${TEST_DIR}/${TEST_NAME}" )
    set_tests_properties( ${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1 )
endforeach()
