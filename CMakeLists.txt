## CMake project file for pcsim

################################################################################
## PROJECT: pcsim
################################################################################
cmake_minimum_required (VERSION 3.8)
project (pcsim VERSION 2.2.3 LANGUAGES CXX )
# note: the following is old-style cmake. But
set (PCSIM_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set (PCSIM_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set (PCSIM_PATCH_VERSION ${PROJECT_VERSION_PATCH})
set (PCSIM_TWEAK_VERSION ${PROJECT_VERSION_TWEAK})
set (PCSIM_VERSION 
  ${PCSIM_MAJOR_VERSION}.${PCSIM_MINOR_VERSION}.${PCSIM_PATCH_VERSION})
set (PCSIM_SOVERSION ${PCSIM_MAJOR_VERSION})

option(COMPILE_FOR_BROADWELL "Enable the compiler flags for broadwell avx2 support" OFF)
option(COMPILE_FOR_HASWELL   "Enable the compiler flags for haswell avx2 support" OFF)
option(COMPILE_FOR_KNL       "Enable the compiler flags for KNL instruction set support" OFF)


################################################################################
## CMAKE Settings 
################################################################################
include(cmake/options.cmake)
include(cmake/compiler.cmake)

################################################################################
## Add subdirectories
################################################################################
set (TARGETS "")

add_subdirectory(pcsim/core)
add_subdirectory(pcsim/physics)
add_subdirectory(pcsim/proc)
add_subdirectory(pcsim/gen)
add_subdirectory(pcsim/root)
add_subdirectory(program)
#add_subdirectory(python)

################################################################################
## Prepare the project to be EXPORTed, generate configuration files
################################################################################
include(cmake/export.cmake)
