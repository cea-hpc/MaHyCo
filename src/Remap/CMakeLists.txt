add_library(Remap
DualPhaseRemap.cc
DualPhaseRemapArcane.cc
IRemap.h
RemapArcane.axl
RemapArcaneFinal.cc
RemapArcaneService.cc
RemapArcaneService.h
RemapADI.axl
RemapADIFinal.cc
RemapADIService.cc
RemapADIService.h
RemapALE.axl
RemapALEService.cc
RemapALEService.h
Remap-ordre-3.cc
RemapArcane-ordre-3.cc
UtilesRemapALE.cc
UtilesRemap.cc
UtilesRemapArcane.cc
)

target_include_directories(Remap PRIVATE ..)
target_link_libraries(Remap PUBLIC Arcane::arcane_core)
target_link_libraries(Remap PRIVATE mahyco_config)
install(TARGETS Remap EXPORT ${MAHYCO_EXPORT_TARGET})

arcane_accelerator_enable()
arcane_accelerator_add_source_files(DualPhaseRemap.cc)
arcane_accelerator_add_source_files(DualPhaseRemapArcane.cc)
arcane_accelerator_add_source_files(RemapADIFinal.cc)
arcane_accelerator_add_source_files(RemapADIService.cc)
arcane_accelerator_add_source_files(RemapArcaneFinal.cc)
arcane_accelerator_add_source_files(RemapArcaneService.cc)
arcane_accelerator_add_source_files(RemapALEService.cc)
arcane_accelerator_add_source_files(Remap-ordre-3.cc)
arcane_accelerator_add_source_files(RemapArcane-ordre-3.cc)
arcane_accelerator_add_source_files(UtilesRemapALE.cc)
arcane_accelerator_add_source_files(UtilesRemap.cc)
arcane_accelerator_add_source_files(UtilesRemapArcane.cc)
arcane_accelerator_add_to_target(Remap)

arcane_target_add_axl(Remap RemapArcane)
arcane_target_add_axl(Remap RemapADI)
arcane_target_add_axl(Remap RemapALE)
