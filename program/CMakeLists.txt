################################################################################
## CMAKE Settings
################################################################################

## function to build our programs
function (BUILD_PROGRAM EXE)
  ## program name
  set (PROG "${EXE}")
  ## Sources and install headers
  set (SOURCES "${EXE}.cc")
  ## compile and link
  add_executable(${PROG} ${SOURCES})
  target_link_libraries(${PROG}
    pcsim_core pcsim_physics pcsim_gen_process pcsim_gen_beam pcsim_gen_decay
    pcsim_gen_detector
    pcsim_core pcsim_physics pcsim_gen_process pcsim_gen_beam pcsim_gen_decay
    pcsim_gen_detector
  ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
  set_target_properties(${PROG} PROPERTIES VERSION ${PCSIM_VERSION})
  ## install instructions
  install(TARGETS ${PROG}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin)
endfunction ()

################################################################################
## Include directories
################################################################################
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

################################################################################
## External Libraries
################################################################################
## ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
## boost
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

################################################################################
## Build all programs
################################################################################
BUILD_PROGRAM(pcsim-lp_gamma)
