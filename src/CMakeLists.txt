#set(CMAKE_CXX_FLAGS "-Werror -Wall -Wextra -pedantic -Wnon-virtual-dtor -g " CACHE STRING "" FORCE)

add_library(mahyco_config INTERFACE)
if (WANT_PROF_ACC)
    target_compile_definitions(mahyco_config INTERFACE PROF_ACC)
    # Pour avoir des #warning explicites
    target_compile_definitions(mahyco_config INTERFACE P4GPU_HAS_WARNING_INFO)
endif ()

find_package(Arcane 3.7 REQUIRED)

# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

set(MAHYCO_EXPORT_TARGET ${PROJECT_NAME}Targets)

add_executable(Mahyco
        MahycoAnnexe.cc
        MahycoConfig.cmake.in
        MahycoModule.cc
        MahycoModule.h
        main.cc
        PrepareRemap.cc
        TypesMahyco.h)

arcane_generate_axl(Mahyco)
arcane_add_arcane_libraries_to_target(Mahyco)
target_include_directories(Mahyco PRIVATE . ${CMAKE_CURRENT_BINARY_DIR})
target_sources(Mahyco PRIVATE MahycoAnnexe.cc)
target_sources(Mahyco PRIVATE PrepareRemap.cc)

install(TARGETS Mahyco EXPORT ${MAHYCO_EXPORT_TARGET})
configure_file(Mahyco.config ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
configure_file(Mahyco.arc ${CMAKE_CURRENT_BINARY_DIR} @ONLY)

add_subdirectory(casTest)

# Equation of State
add_subdirectory(eos)

# Remap
add_subdirectory(Remap)

# Cartesian
add_subdirectory(cartesian)

# Services communs pour les accélérateurs
add_subdirectory(accenv)

# Message Passing
add_subdirectory(msgpass)

# Pour la gestion des accélérateurs et des sources à compiler avec le
# compilateur pour accélérateur
arcane_accelerator_enable()
arcane_accelerator_add_source_files(MahycoModule.cc)
arcane_accelerator_add_source_files(MahycoAnnexe.cc)
arcane_accelerator_add_source_files(PrepareRemap.cc)
arcane_accelerator_add_to_target(Mahyco)

target_link_libraries(Mahyco PRIVATE mahyco_config casTest Remap cartesian eos accenv msgpass)

# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E ${CMAKE_CURRENT_BINARY_DIR}/Mahyco)


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
    add_test(NAME "${SEQPR}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${SEQPR} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

# sur 4 coeurs
foreach (TEST_NAME ${TEST_NAMES_PARA_4})
    add_test(NAME "${PARA4}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${PARA4} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

# sur 8 coeurs
foreach (TEST_NAME ${TEST_NAMES_PARA_8})
    add_test(NAME "${PARA8}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${PARA8} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

# sequentiel
foreach (TEST_NAME ${TEST_NAMES})
    add_test(NAME "${SEQ}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${SEQ} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
endforeach ()

if (WANT_CUDA)
    # mono-processus 1 GPU
    foreach (TEST_NAME ${TEST_NAMES_CUDA_1})
        add_test(NAME "${CUDA1}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${CUDA1} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
    endforeach ()

    # sur 4 GPUs (1 GPU par processus)
    foreach (TEST_NAME ${TEST_NAMES_CUDA_4})
        add_test(NAME "${CUDA4}_${TEST_NAME}" COMMAND ${TEST_DIR}/../NONREGRESSION/launch_test.sh ${CMAKE_CURRENT_BINARY_DIR}/Mahyco ${TEST_DIR}/../NONREGRESSION/${TEST_NAME} ${CUDA4} ${MPI_LAUNCHER} ${BASETMPDIR_NONREG})
    endforeach ()
endif ()


# Configuration
set(config_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
set(namespace "Alien::")

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Use:
#   * PROJECT_VERSION
write_basic_package_version_file(
        "${version_config}" COMPATIBILITY SameMajorVersion
)

set(MAHYCO_MODULE_PATH "${config_install_dir}")

# Configure '<PROJECT-NAME>Config.cmake'
# Use variables:
#   * TARGETS_EXPORT_NAME
#   * PROJECT_NAME
configure_package_config_file(
        "MahycoConfig.cmake.in"
        "${project_config}"
        INSTALL_DESTINATION "${config_install_dir}"
        PATH_VARS MAHYCO_MODULE_PATH
)

set(Mahyco_project_config "${project_config}" CACHE FILEPATH "Mahyco project config" FORCE)
get_filename_component(Mahyco_config_build_dir "${project_config}" DIRECTORY CACHE)
set(Mahyco_config_install_dir "${config_install_dir}" CACHE STRING "Mahyco config install dir" FORCE)
set(Mahyco_VERSION "${Mahyco_VERSION}" CACHE STRING "Mahyco version" FORCE)

install(FILES "${project_config}" "${version_config}"
        DESTINATION "${config_install_dir}"
        )

install(EXPORT ${MAHYCO_EXPORT_TARGET}
        NAMESPACE ${namespace}
        DESTINATION ${config_install_dir})
