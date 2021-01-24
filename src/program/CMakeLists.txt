################################################################################
## CMAKE Settings
################################################################################


## function to build our programs
function (BUILD_PROGRAM PROG)
  ## Sources and install headers
  set (SOURCES "${PROG}.cc")
  ## compile and link
  add_executable(${PROG} ${SOURCES})
  target_link_libraries(${PROG}
    lager_core lager_physics lager_proc lager_gen
    lager_core lager_physics lager_proc lager_gen
    ${PHOTOSPP_LIBRARIES}
    ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${JPACPHOTO_LIBRARIES})
#target_compile_options(${PROG} PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")
  target_compile_options(${PROG} PUBLIC "$<$<CONFIG:RELEASE>:${PROJECT_EXTRA_CXX_FLAGS}>")
  ## install instructions
  install(TARGETS ${PROG}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT bin)
endfunction ()

################################################################################
## Include directories
################################################################################
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

################################################################################
## Build all programs
################################################################################
BUILD_PROGRAM(lager)
