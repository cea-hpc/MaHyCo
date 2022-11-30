cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
# Minimum version for Arcane

project(Mahyco LANGUAGES C CXX VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_SHARED_LIBS TRUE)

# Pour utiliser l'implémentation Arcane du cartésien
option(WANT_IMPL_CART_ARCANE "Utilisation de l'implémentation Arcane:: dans arcane/cea/" FALSE)

# Pour la gestion du profiling sur accélérateur
option(WANT_PROF_ACC "Active le profiling des points d'entrées pour accélérateur" FALSE)

add_subdirectory(src)

enable_testing()
add_subdirectory(NONREGRESSION)

include(GNUInstallDirs)

install(DIRECTORY NONREGRESSION/
  DESTINATION share/cases
  PATTERN "output" EXCLUDE
  )
