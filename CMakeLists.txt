## CMake project file for liege

################################################################################
## PROJECT: liege
################################################################################
cmake_minimum_required (VERSION 3.8)
project (liege VERSION 3.0.0 LANGUAGES CXX )
# note: the following is old-style cmake. But
set (LIEGE_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set (LIEGE_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set (LIEGE_PATCH_VERSION ${PROJECT_VERSION_PATCH})
set (LIEGE_TWEAK_VERSION ${PROJECT_VERSION_TWEAK})
set (LIEGE_VERSION 
  ${LIEGE_MAJOR_VERSION}.${LIEGE_MINOR_VERSION}.${LIEGE_PATCH_VERSION})
set (LIEGE_SOVERSION ${LIEGE_MAJOR_VERSION})

option(COMPILE_FOR_BROADWELL "Enable the compiler flags for broadwell avx2 support" OFF)
option(COMPILE_FOR_HASWELL   "Enable the compiler flags for haswell avx2 support" OFF)
option(COMPILE_FOR_KNL       "Enable the compiler flags for KNL instruction set support" OFF)

################################################################################
## CMAKE Settings 
################################################################################
include(cmake/options.cmake)
include(cmake/compiler.cmake)
include(cmake/rpath.cmake)

################################################################################
## Add subdirectories
################################################################################
set (TARGETS "")

add_subdirectory(liege/core)
add_subdirectory(liege/physics)
add_subdirectory(liege/proc)
add_subdirectory(liege/gen)
add_subdirectory(program)

################################################################################
## Prepare the project to be EXPORTed, generate configuration files
################################################################################
include(cmake/export.cmake)
