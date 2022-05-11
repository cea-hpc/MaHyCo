# Indique qu'on veut des bibliothèques dynamiques (indispensable avec les services)
set(BUILD_SHARED_LIBS TRUE)
#set(Arcane_ROOT "/home/perlatj/ARCANE/arcane_build/arcframework/arcane/arcane/install/2.22.0.0" CACHE PATH "Arcane cmake root path" FORCE)
#set(Arcane_ROOT "/home/meltzb/workspace/arcane" CACHE PATH "Arcane cmake root path" FORCE)
set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Verbose makefile?" FORCE)
#set(CMAKE_C_COMPILER "/usr/bin/cc" CACHE PATH "C Compiler" FORCE)
#set(CMAKE_CXX_COMPILER "/usr/bin/c++" CACHE PATH "C++ Compiler" FORCE)
set(CMAKE_CXX_FLAGS_INIT "  " CACHE STRING "Default flags for C++ compiler" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "" CACHE STRING "Default flags for linker" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS "" CACHE STRING "Default flags for linker" FORCE)

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8
# End:
