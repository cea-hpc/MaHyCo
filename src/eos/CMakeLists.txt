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

arcane_accelerator_enable()
arcane_accelerator_add_source_files(perfectgas/PerfectGasEOSService.cc)
arcane_accelerator_add_source_files(stiffenedgas/StiffenedGasEOSService.cc)
#arcane_accelerator_add_source_files(stdperfectgas/StdPerfectGasEOSService.cc)
#arcane_accelerator_add_source_files(stdperfectgas/PhyVarType.cc)
arcane_accelerator_add_to_target(PerfectGas)
arcane_accelerator_add_to_target(StiffenedGas)
#arcane_accelerator_add_to_target(StdPerfectGas)

arcane_target_add_axl(PerfectGas perfectgas/PerfectGasEOS)
arcane_target_add_axl(StiffenedGas stiffenedgas/StiffenedGasEOS)
arcane_target_add_axl(StdPerfectGas stdperfectgas/StdPerfectGasEOS)

add_library(eos INTERFACE)
target_link_libraries(eos INTERFACE PerfectGas StiffenedGas StdPerfectGas)
