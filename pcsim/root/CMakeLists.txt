################################################################################
## CMAKE Settings
################################################################################
set (SUBDIR "root")
set (LIBRARY "${PROJECT_NAME}_${SUBDIR}")
set (TARGETS ${TARGETS} ${LIBRARY} PARENT_SCOPE)

################################################################################
## Sources and install headers
################################################################################
file (GLOB SOURCES "dataframe/*.cc" "physics/*.cc")
file (GLOB DATAFRAME_HEADERS "dataframe/*.hh")
file (GLOB PHYSICS_HEADERS "physics/*.hh")

################################################################################
## Include directories
################################################################################
include_directories("${PROJECT_SOURCE_DIR}")
#include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_BINARY_DIR}")

################################################################################
## External Libraries
################################################################################

## ANL/LCIO2 (optional)
find_package(podio QUIET)
find_package(lcio2 QUIET)
if(podio_FOUND AND lcio2_FOUND)
  message("LCIO2/PODIO found, including ANL ROOT bindings")
  set(SUPPORT_ANL 1)
  include_directories(${LCIO2_INCLUDE_DIRS})
else()
  message("No LCIO2/PODIO found, not building ANL ROOT support")
  set(SUPPORT_ANL 0)
endif()

################################################################################
## Configure files
################################################################################
configure_file(setup_root.cc.in
    ${CMAKE_CURRENT_BINARY_DIR}/setup_root.cc)
set(CONFIG_ROOT ${CMAKE_CURRENT_BINARY_DIR}/setup_root.cc)

configure_file(pcsim-root.in 
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim-root)
set(CONFIG_EXE ${CMAKE_CURRENT_BINARY_DIR}/pcsim-root)

configure_file(config.hh.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hh)
set(CONFIG_HEADERS "${CMAKE_CURRENT_BINARY_DIR}/config.hh")

################################################################################
## Compile and Link
################################################################################
add_library(${LIBRARY} SHARED ${SOURCES})
target_link_libraries(${LIBRARY} 
                      ${PROJECT_NAME}_core ${PROJECT_NAME}_physics
                      ${PROJECT_NAME}_proc
                      ${ROOT_GenVector_LIBRARY} ${ROOT_EG_LIBRARY} ${ROOT_LIBRARIES} 
                      ${Boost_LIBRARIES})
set_target_properties(${LIBRARY} PROPERTIES 
  VERSION ${PCSIM_VERSION} 
  SOVERSION ${PCSIM_SOVERSION})

################################################################################
## Export and Install
################################################################################
install(TARGETS ${LIBRARY}
  EXPORT ${PROJECT_NAME}-targets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib)
  
install(PROGRAMS
  ${CONFIG_EXE}
  DESTINATION "${INSTALL_BIN_DIR}"
  COMPONENT bin)
install(FILES
  ${CONFIG_ROOT} 
  DESTINATION "${INSTALL_LIB_DIR}/${PROJECT_NAME}"
  COMPONENT lib)

install(FILES 
  ${CONFIG_HEADERS} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}" 
  COMPONENT dev)
install(FILES 
  ${DATAFRAME_HEADERS} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/dataframe" 
  COMPONENT dev)
install(FILES 
  ${PHYSICS_HEADERS} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/physics" 
  COMPONENT dev)
