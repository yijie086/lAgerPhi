## CMake project file for pcsim

################################################################################
## PROJECT: pcsim
################################################################################
cmake_minimum_required (VERSION 2.8)
project (pcsim)
set (PCSIM_MAJOR_VERSION 1)
set (PCSIM_MINOR_VERSION 1)
set (PCSIM_PATCH_VERSION 1)
set (PCSIM_VERSION 
  ${PCSIM_MAJOR_VERSION}.${PCSIM_MINOR_VERSION}.${PCSIM_PATCH_VERSION})
set (PCSIM_SOVERSION ${PCSIM_MAJOR_VERSION})

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
add_subdirectory(pcsim/gen)
add_subdirectory(pcsim/physics)
add_subdirectory(program)

################################################################################
## Prepare the project to be EXPORTed, generate configuration files
################################################################################
include(cmake/export.cmake)
