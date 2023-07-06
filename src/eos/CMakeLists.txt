add_library(PerfectGas
perfectgas/PerfectGasEOSService.cc
perfectgas/PerfectGasEOSService.h
IEquationOfState.h)

target_include_directories(PerfectGas PRIVATE ..)
target_link_libraries(PerfectGas PUBLIC Arcane::arcane_core)
target_link_libraries(PerfectGas PRIVATE mahyco_config)
install(TARGETS PerfectGas EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(StiffenedGas
  stiffenedgas/StiffenedGasEOSService.cc
  stiffenedgas/StiffenedGasEOSService.h
  IEquationOfState.h)

target_include_directories(StiffenedGas PRIVATE ..)
target_link_libraries(StiffenedGas PUBLIC Arcane::arcane_core)
target_link_libraries(StiffenedGas PRIVATE mahyco_config)
install(TARGETS StiffenedGas EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(StdPerfectGas
stdperfectgas/StdPerfectGasEOSService.cc
stdperfectgas/StdPerfectGasEOSService.h
stdperfectgas/PhyVarType.cc
stdperfectgas/PhyVarType.h
stdperfectgas/EosTypes.h
IEquationOfState.h)

target_include_directories(StdPerfectGas PRIVATE ..)
target_link_libraries(StdPerfectGas PUBLIC Arcane::arcane_core)
target_link_libraries(StdPerfectGas PRIVATE mahyco_config)
install(TARGETS StdPerfectGas EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(StdPerfectGasAcc1
stdperfectgasacc1/StdPerfectGasAcc1EOSService.cc
stdperfectgasacc1/StdPerfectGasAcc1EOSService.h
stdperfectgasacc1/PhyVarType.cc
stdperfectgasacc1/PhyVarType.h
stdperfectgasacc1/EosTypes.h
stdperfectgasacc1/Mallocator.h
IEquationOfState.h)

target_include_directories(StdPerfectGasAcc1 PRIVATE ..)
target_link_libraries(StdPerfectGasAcc1 PUBLIC Arcane::arcane_core)
target_link_libraries(StdPerfectGasAcc1 PRIVATE mahyco_config)
find_package(Boost)
if(Boost_FOUND)
  target_compile_definitions(StdPerfectGasAcc1 PRIVATE MAHYCO_HAS_BOOST)
  target_link_libraries(StdPerfectGasAcc1 PRIVATE Boost::headers)
endif()
install(TARGETS StdPerfectGasAcc1 EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(StdPerfectGasAcc2
stdperfectgasacc2/StdPerfectGasAcc2EOSService.cc
stdperfectgasacc2/StdPerfectGasAcc2EOSService.h
stdperfectgasacc2/PhyVarType.cc
stdperfectgasacc2/PhyVarType.h
stdperfectgasacc2/PhyVars.cc
stdperfectgasacc2/PhyVars.h
stdperfectgasacc2/EosTypes.h
IEquationOfState.h)

target_include_directories(StdPerfectGasAcc2 PRIVATE ..)
target_link_libraries(StdPerfectGasAcc2 PUBLIC Arcane::arcane_core)
target_link_libraries(StdPerfectGasAcc2 PRIVATE mahyco_config)
install(TARGETS StdPerfectGasAcc2 EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(QnewtStdPerfectGas
qnewt_stdperfectgas/QnewtStdPerfectGasEOSService.cc
qnewt_stdperfectgas/QnewtStdPerfectGasEOSService.h
qnewt_stdperfectgas/PhyVarType.cc
qnewt_stdperfectgas/PhyVarType.h
qnewt_stdperfectgas/EosTypes.h
qnewt_stdperfectgas/Utils.h
IEquationOfState.h)

target_include_directories(QnewtStdPerfectGas PRIVATE ..)
target_link_libraries(QnewtStdPerfectGas PUBLIC Arcane::arcane_core)
target_link_libraries(QnewtStdPerfectGas PRIVATE mahyco_config)
install(TARGETS QnewtStdPerfectGas EXPORT ${MAHYCO_EXPORT_TARGET})

arcane_accelerator_enable()

arcane_accelerator_add_source_files(perfectgas/PerfectGasEOSService.cc)
arcane_accelerator_add_source_files(stiffenedgas/StiffenedGasEOSService.cc)

#arcane_accelerator_add_source_files(stdperfectgas/StdPerfectGasEOSService.cc)
#arcane_accelerator_add_source_files(stdperfectgas/PhyVarType.cc)

arcane_accelerator_add_source_files(stdperfectgasacc1/StdPerfectGasAcc1EOSService.cc)
arcane_accelerator_add_source_files(stdperfectgasacc1/PhyVarType.cc)

arcane_accelerator_add_source_files(stdperfectgasacc2/StdPerfectGasAcc2EOSService.cc)
arcane_accelerator_add_source_files(stdperfectgasacc2/PhyVarType.cc)
arcane_accelerator_add_source_files(stdperfectgasacc2/PhyVars.cc)

#arcane_accelerator_add_source_files(qnewt_stdperfectgas/QnewtStdPerfectGasEOSService.cc)
#arcane_accelerator_add_source_files(qnewt_stdperfectgas/PhyVarType.cc)

arcane_accelerator_add_to_target(PerfectGas)
arcane_accelerator_add_to_target(StiffenedGas)
#arcane_accelerator_add_to_target(StdPerfectGas)
arcane_accelerator_add_to_target(StdPerfectGasAcc1)
arcane_accelerator_add_to_target(StdPerfectGasAcc2)
#arcane_accelerator_add_to_target(QnewtStdPerfectGas)

arcane_target_add_axl(PerfectGas perfectgas/PerfectGasEOS)
arcane_target_add_axl(StiffenedGas stiffenedgas/StiffenedGasEOS)
arcane_target_add_axl(StdPerfectGas stdperfectgas/StdPerfectGasEOS)
arcane_target_add_axl(StdPerfectGasAcc1 stdperfectgasacc1/StdPerfectGasAcc1EOS)
arcane_target_add_axl(StdPerfectGasAcc2 stdperfectgasacc2/StdPerfectGasAcc2EOS)
arcane_target_add_axl(QnewtStdPerfectGas qnewt_stdperfectgas/QnewtStdPerfectGasEOS)

add_library(eos INTERFACE)
target_link_libraries(eos INTERFACE 
    PerfectGas 
    StiffenedGas 
    StdPerfectGas
    StdPerfectGasAcc1
    StdPerfectGasAcc2
    QnewtStdPerfectGas
    )
