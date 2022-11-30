# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E  $<TARGET_FILE:Mahyco>)

# Pour les test NONREG
# Test name is also test directory
set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}")

set(TEST_NAMES_PARA_4
		"CAS_2D_Lag_BiSodCaseX"
		"CAS_2D_Lag_BiSodCaseY"
		"CAS_2D_Lag_SodCaseY"
		"CAS_Lag_BiSodCaseZ"
		"CAS_Lag_SodCaseX"
		"CAS_BiSodCaseX"
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
set(TEST_NAMES_PARA_8
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
set(TEST_NAMES
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

set(TEST_NAMES_PR
		"CAS_2D_Lag_SodCaseX"
		"CAS_SodCaseX"
		"CAS_MonoRiderTx"
		"CAS_RiderT45"
		)

set(TEST_NAMES_CUDA_1
		"CAS_BiSodCaseX"
		"CAS_BiSodCaseY"
		"CAS_BiSodCaseZ"
		"CAS_Lag_BiSodCaseX"
		"CAS_RiderRotation"
		)

set(TEST_NAMES_CUDA_4
		"CAS_BiSodCaseX"
		)

set(PARA4 "para_4")
set(PARA8 "para_8")
set(SEQ "seq")
set(SEQPR "seq_pr")
set(CUDA1 "cuda_1")
set(CUDA4 "cuda_4")

# Pour définir le launcher MPI
set(MPI_LAUNCHER mpiexec CACHE STRING "Définit le launcher MPI")

# Pour définir la base du répertoire temporaire d'exécution des tests de non-regression
set(BASETMPDIR_NONREG "/tmp" CACHE PATH "Définit la base du répertoire temporaire d'exécution des tests de non-regression")

# protection-reprise
foreach (TEST_NAME ${TEST_NAMES_PR})
	add_test(NAME "${SEQPR}_${TEST_NAME}" COMMAND ${TEST_DIR}/launch_test.sh  $<TARGET_FILE:Mahyco> ${TEST_DIR}/${TEST_NAME} ${SEQPR} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

# sur 4 coeurs
foreach (TEST_NAME ${TEST_NAMES_PARA_4})
	add_test(NAME "${PARA4}_${TEST_NAME}" COMMAND ${TEST_DIR}/launch_test.sh  $<TARGET_FILE:Mahyco> ${TEST_DIR}/${TEST_NAME} ${PARA4} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

# sur 8 coeurs
foreach (TEST_NAME ${TEST_NAMES_PARA_8})
	add_test(NAME "${PARA8}_${TEST_NAME}" COMMAND ${TEST_DIR}/launch_test.sh  $<TARGET_FILE:Mahyco> ${TEST_DIR}/${TEST_NAME} ${PARA8} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

# sequentiel
foreach (TEST_NAME ${TEST_NAMES})
	add_test(NAME "${SEQ}_${TEST_NAME}" COMMAND ${TEST_DIR}/launch_test.sh  $<TARGET_FILE:Mahyco> ${TEST_DIR}/${TEST_NAME} ${SEQ} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

if (WANT_CUDA)
	# mono-processus 1 GPU
	foreach (TEST_NAME ${TEST_NAMES_CUDA_1})
		add_test(NAME "${CUDA1}_${TEST_NAME}" COMMAND ${TEST_DIR}/launch_test.sh  $<TARGET_FILE:Mahyco> ${TEST_DIR}/${TEST_NAME} ${CUDA1} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
	endforeach ()

	# sur 4 GPUs (1 GPU par processus)
	foreach (TEST_NAME ${TEST_NAMES_CUDA_4})
		add_test(NAME "${CUDA4}_${TEST_NAME}" COMMAND ${TEST_DIR}/launch_test.sh  $<TARGET_FILE:Mahyco> ${TEST_DIR}/${TEST_NAME} ${CUDA4} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
	endforeach ()
endif ()

