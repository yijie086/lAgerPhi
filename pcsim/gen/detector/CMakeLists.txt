################################################################################
## CMAKE Settings
################################################################################
set (SUBDIR "gen_detector")
set (LIBRARY "${PROJECT_NAME}_${SUBDIR}")
set (TARGETS ${TARGETS} ${LIBRARY} PARENT_SCOPE)

################################################################################
## Sources and install headers
################################################################################
file (GLOB SOURCES "*.cc")
file (GLOB HEADERS "*.hh")

################################################################################
## Include directories
################################################################################
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

################################################################################
## Add definitions
## TODO: more robust installable data dir
################################################################################
add_definitions(-DEIC_ACCEPTANCE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/eic-fastmc/acceptance_files/\" -DJLEIC_IMPL)

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
## Compile and Link
################################################################################
add_library(${LIBRARY} ${SOURCES})
set_target_properties(${LIBRARY} PROPERTIES 
  VERSION ${PCSIM_VERSION} 
  SOVERSION ${PCSIM_SOVERSION}
  PUBLIC_HEADER "${HEADERS}")

################################################################################
## Export and Install
################################################################################
install(TARGETS ${LIBRARY}
  EXPORT ${PROJECT_NAME}-targets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}" COMPONENT dev)
