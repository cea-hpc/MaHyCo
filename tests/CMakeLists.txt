# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )

set( TEST_NAMES "SOD_X_bi_mat_Superbee_PenteBorne_mix"
                "SOD_X_bi_mat_MinMod_simple"
                "SOD_X_bi_mat_Superbee_PenteBorne_simple"
                "SOD_X_bi_mat_Superbee_simple"
                "SOD_Y_bi_mat_MinMod_simple"
                "SOD_Y_bi_mat_Superbee_PenteBorne_mixte"
                "SOD_Y_bi_mat_Superbee_PenteBorne_pure"
                "SOD_Y_bi_mat_VanLeer_PenteBorne_pure"
                "SOD_X_Mono_mat_MinModG_simple"
                "SOD_X_Mono_mat_MinMod_simple"
                "SOD_X_Mono_mat_SuperBee_simple"
                "SOD_X_Mono_mat_VanLeer_simple" )

foreach( TEST_NAME ${TEST_NAMES} )
    add_test( NAME ${TEST_NAME}
              COMMAND "${TEST_DIR}/launch_test.sh" $<TARGET_FILE:mahyco> "${TEST_DIR}/${TEST_NAME}" )
    set_tests_properties( ${TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1 )
endforeach()