add_library(msgpass Algo1SyncDataD.cc
        Algo1SyncDataD.h
        Algo1SyncDataDH.cc
        Algo1SyncDataDH.h
        ComputeAndSync.h
        EnumerateEnvAndSync.h
        IAlgo1SyncData.h
        IsCommDeviceAware.cc
        MeshVariableSynchronizerList.cc
        MeshVariableSynchronizerList.h
        MsgPassInit.cc
        PackTransfer.cc
        PackTransfer.h
        SyncAndCompute.h
        SyncBuffers.cc
        SyncBuffers.h
        SyncEnvIndexes.cc
        SyncEnvIndexes.h
        SyncItems.cc
        SyncItems.h
        VarSyncAlgo1.cc
        VarSyncAlgo1.h
        VarSyncMng.cc
        VarSyncMng.h
        VarSyncMngOptions.h)

target_include_directories(msgpass PRIVATE ..)
target_link_libraries(msgpass PUBLIC Arcane::arcane_core)
# Pour MPI
#find_package(MPI)
#if(MPI_FOUND)
#  target_compile_definitions(msgpass PRIVATE OMPI_SKIP_MPICXX)
#  target_compile_definitions(msgpass PRIVATE MSG_PASS_HAS_MPI)
#  target_include_directories(msgpass PRIVATE "${MPI_INCLUDE_PATH}")
#  find_path(P4GPU_MPI_EXT NAMES "mpi-ext.h" PATHS "${MPI_C_HEADER_DIR}")
#  if(P4GPU_MPI_EXT)
#    target_compile_definitions(msgpass PRIVATE MSG_PASS_HAS_MPI_EXT)
#    target_link_libraries(msgpass PRIVATE mpi)
#  endif()
#endif()
# On force is_comm_device_aware() == true
target_compile_definitions(msgpass PRIVATE MSG_PASS_HAS_COMM_DEVICE_AWARE)
# fin MPI
install(TARGETS msgpass EXPORT ${MAHYCO_EXPORT_TARGET})

arcane_accelerator_enable()
arcane_accelerator_add_source_files(VarSyncMng.cc)
arcane_accelerator_add_source_files(IsCommDeviceAware.cc)
arcane_accelerator_add_source_files(MsgPassInit.cc)
arcane_accelerator_add_source_files(PackTransfer.cc)
arcane_accelerator_add_source_files(SyncItems.cc)
arcane_accelerator_add_source_files(SyncEnvIndexes.cc)
arcane_accelerator_add_source_files(SyncBuffers.cc)
arcane_accelerator_add_source_files(MeshVariableSynchronizerList.cc)
arcane_accelerator_add_source_files(Algo1SyncDataD.cc)
arcane_accelerator_add_source_files(Algo1SyncDataDH.cc)
arcane_accelerator_add_source_files(VarSyncAlgo1.cc)
arcane_accelerator_add_to_target(msgpass)
