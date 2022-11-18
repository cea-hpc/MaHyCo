# -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
if (NOT EXAMPLE_NAME)
  message(FATAL_ERROR "Variable EXAMPLE_NAME not defined")
endif()

find_package(Arcane)
# Indique la version minimale de Arcane nécessaire
# À partir de la version 3.8 de Arcane, il sera possible
# d'utiliser directement 'find_package(Arcane 3.8 REQUIRED)'
if (Arcane_VERSION VERSION_LESS "3.6")
  message(FATAL_ERROR "Arcane version '${Arcane_VERSION}' is too old."
    " Version 3.6+ is required")
endif()

include(${Arcane_DIR}/ArcaneDotNet.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/mahyco.utils.cmake)

add_executable(${EXAMPLE_NAME} ${EXAMPLE_NAME}Module.cc main.cc ${EXAMPLE_NAME}_axl.h)

arcane_generate_axl(${EXAMPLE_NAME})
configure_file(${EXAMPLE_NAME}.config ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
configure_file(${EXAMPLE_NAME}.arc ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
arcane_add_arcane_libraries_to_target(${EXAMPLE_NAME})
target_include_directories(${EXAMPLE_NAME} PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
