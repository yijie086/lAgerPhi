## CMake project file for lager

################################################################################
## PROJECT: lager
################################################################################
cmake_minimum_required (VERSION 3.8)
project (lager VERSION 3.2.1 LANGUAGES CXX )
# note: the following is old-style cmake. But
set (LAGER_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set (LAGER_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set (LAGER_PATCH_VERSION ${PROJECT_VERSION_PATCH})
set (LAGER_TWEAK_VERSION ${PROJECT_VERSION_TWEAK})
set (LAGER_VERSION 
  ${LAGER_MAJOR_VERSION}.${LAGER_MINOR_VERSION}.${LAGER_PATCH_VERSION})
set (LAGER_SOVERSION ${LAGER_MAJOR_VERSION})

option(COMPILE_FOR_BROADWELL "Enable the compiler flags for broadwell avx2 support" OFF)
option(COMPILE_FOR_HASWELL   "Enable the compiler flags for haswell avx2 support" OFF)
option(COMPILE_FOR_KNL       "Enable the compiler flags for KNL instruction set support" OFF)

################################################################################
## CMAKE Settings 
################################################################################
include(cmake/options.cmake)
include(cmake/compiler.cmake)
include(cmake/rpath.cmake)

add_definitions(-DLAGER_VERSION=\"${LAGER_VERSION}\")

################################################################################
## Add subdirectories
################################################################################
set (TARGETS "")

add_subdirectory(src/lager/core)
add_subdirectory(src/lager/physics)
add_subdirectory(src/lager/proc)
add_subdirectory(src/lager/gen)
add_subdirectory(src/program)

################################################################################
## Prepare the project to be EXPORTed, generate configuration files
################################################################################
include(cmake/export.cmake)
