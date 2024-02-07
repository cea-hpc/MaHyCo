# Services communs pour les accélérateurs
add_library(accenv
    AcceleratorUtils.h
    AccEnv.axl
    AccEnvDefault.axl
    AccEnvDefaultOptions.h
    AccEnvDefaultService.cc
    AccEnvDefaultService.h
    AccEnvModule.cc
    AccEnvModule.h
    BufAddrMng.h
    IAccEnv.h
    IncidentHandler.h
    IndexSelecter.h
    MultiEnvMng.cc
    MultiEnvMng.h
    MultiEnvUtils.h
    ProfAcc.cc
    ProfAcc.h
    SingletonIAccEnv.h
    SortSmall.h
    UtilsAtomic.h)

target_include_directories(accenv PRIVATE ..)
target_link_libraries(accenv PUBLIC Arcane::arcane_core)
target_link_libraries(accenv PRIVATE mahyco_config)
install(TARGETS accenv EXPORT ${MAHYCO_EXPORT_TARGET})

arcane_accelerator_enable()
arcane_accelerator_add_source_files(AccEnvDefaultService.cc)
arcane_accelerator_add_source_files(AccEnvModule.cc)
arcane_accelerator_add_source_files(MultiEnvMng.cc)
arcane_accelerator_add_source_files(ProfAcc.cc)
arcane_accelerator_add_to_target(accenv)

arcane_target_add_axl(accenv AccEnvDefault)
arcane_target_add_axl(accenv AccEnv)

