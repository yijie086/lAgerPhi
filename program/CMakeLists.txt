################################################################################
## CMAKE Settings
################################################################################


## function to build our programs
function (BUILD_PROGRAM EXE PROG)
  ## Sources and install headers
  set (SOURCES "${EXE}.cc")
  ## compile and link
  add_executable(${PROG} ${SOURCES})
  target_link_libraries(${PROG}
    pcsim_core pcsim_physics pcsim_proc pcsim_gen
    pcsim_core pcsim_physics pcsim_proc pcsim_gen
  ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
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
## Build all programs
################################################################################
BUILD_PROGRAM(pcsim-lp_gamma pcsim)
