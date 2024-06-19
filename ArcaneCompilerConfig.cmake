# Indique qu'on veut des bibliothèques dynamiques (indispensable avec les services)
set(BUILD_SHARED_LIBS TRUE)

# Si aucun chemin versArcane n'est renseigné pour la compilation, on essaie des chemins
# particuliers, sur PC portable et sur calculateur
if (Arcane_ROOT STREQUAL "")
    MESSAGE(STATUS "HOSTNAME = $ENV{HOSTNAME}")
    if ( $ENV{HOSTNAME} STREQUAL "login14")
        message("Passage pour INTI")
        set(Arcane_ROOT "$ENV{HOME}/local_arcane" CACHE PATH "Arcane cmake root path" FORCE)
    else()
        message("Passage pour Ubuntu")
        set(Arcane_ROOT "$ENV{HOME}/Workspace/arcane" CACHE PATH "Arcane cmake root path" FORCE)
    endif()
endif()

message("ArcaneRoot = ${Arcane_ROOT}")
set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Verbose makefile?" FORCE)
set(CMAKE_CXX_FLAGS_INIT "  " CACHE STRING "Default flags for C++ compiler" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "" CACHE STRING "Default flags for linker" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS "" CACHE STRING "Default flags for linker" FORCE)

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8
# End:
