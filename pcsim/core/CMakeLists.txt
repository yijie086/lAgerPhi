cmake_minimum_required (VERSION 3.0)

################################################################################
## CMAKE Settings
################################################################################
set (SUBDIR "core")
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
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

################################################################################
## Configure file
################################################################################
#PCSIM commit version
execute_process(COMMAND git describe 
  OUTPUT_VARIABLE PCSIM_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file(version.hh.in ${CMAKE_CURRENT_BINARY_DIR}/version.hh)
set(HEADERS ${HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/version.hh)

################################################################################
## Compile and Link
################################################################################
add_library(${LIBRARY} SHARED ${SOURCES})
target_link_libraries(${LIBRARY} ${ROOT_EG_LIBRARY} ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
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
