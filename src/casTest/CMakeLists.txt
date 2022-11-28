add_library(RIDER RIDER/RIDERService.cc
  RIDER/RIDERService.h IInitialisations.h)
target_include_directories(RIDER PRIVATE .)
target_link_libraries(RIDER PUBLIC Arcane::arcane_core)
install(TARGETS RIDER EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(SOD SOD/SODService.cc
  SOD/SODService.h IInitialisations.h)
target_include_directories(SOD PRIVATE .)
target_link_libraries(SOD PUBLIC Arcane::arcane_core)
install(TARGETS SOD EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(SEDOV SEDOV/SEDOVService.cc
SEDOV/SEDOVService.h IInitialisations.h)
target_include_directories(SEDOV PRIVATE .)
target_link_libraries(SEDOV PUBLIC Arcane::arcane_core)
install(TARGETS SEDOV EXPORT ${MAHYCO_EXPORT_TARGET})

add_library(OTHER OTHER/OTHERService.cc
OTHER/OTHERService.h IInitialisations.h)
target_include_directories(OTHER PRIVATE .)
target_link_libraries(OTHER PUBLIC Arcane::arcane_core)
install(TARGETS OTHER EXPORT ${MAHYCO_EXPORT_TARGET})

arcane_target_add_axl(RIDER RIDER/RIDER)
arcane_target_add_axl(SOD SOD/SOD)
arcane_target_add_axl(SEDOV SEDOV/SEDOV)
arcane_target_add_axl(OTHER OTHER/OTHER)

add_library(casTest INTERFACE)
target_link_libraries(casTest INTERFACE RIDER SOD SEDOV OTHER)
