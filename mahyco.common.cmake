﻿# -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
if (NOT EXAMPLE_NAME)
  message(FATAL_ERROR "Variable EXAMPLE_NAME not defined")
endif()

find_package(Arcane)

include(${Arcane_DIR}/ArcaneDotNet.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/mahyco.utils.cmake)

add_executable(${EXAMPLE_NAME} 
               main 
               ${EXAMPLE_NAME}Module.cc 
               ${EXAMPLE_NAME}_axl.h 
               time_history/TimeHistoryModule.cc 
               time_history/TimeHistory_axl.h)

arcane_generate_axl(${EXAMPLE_NAME})
arcane_generate_axl(time_history/TimeHistory)

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
