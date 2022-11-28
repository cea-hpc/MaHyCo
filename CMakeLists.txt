cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
# Minimum version for Arcane

set(BUILD_SHARED_LIBS TRUE)

project(ArcaneSamples LANGUAGES CXX)

enable_testing()

add_subdirectory(src)
