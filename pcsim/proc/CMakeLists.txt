################################################################################
## CMAKE Settings
################################################################################
set (SUBDIR "proc")
set (LIBRARY "${PROJECT_NAME}_${SUBDIR}")
set (TARGETS ${TARGETS} ${LIBRARY} PARENT_SCOPE)

################################################################################
## Sources and install headers
################################################################################
file (GLOB SOURCES "decay/*.cc" "detector/*.cc" "reconstruction/*.cc")
file (GLOB HEADERS_DECAY "decay/*.hh")
file (GLOB HEADERS_DETECTORS "detectors/*.hh")
file (GLOB HEADERS_RC "reconstruction/*.hh")

set (HEADERS ${HEADERS_DECAY} ${HEADERS_DETECTORS})

################################################################################
## Include directories
################################################################################
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

################################################################################
## Add definitions
## TODO: more robust installable data dir
################################################################################
add_definitions(-DSOLID_ACCEPTANCE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/detector/solid-fastmc/\")

################################################################################
## Compile and Link
################################################################################
add_library(${LIBRARY} SHARED ${SOURCES})
target_compile_options(${LIBRARY} PUBLIC ${PROJECT_EXTRA_CXX_FLAGS})
target_link_libraries(${LIBRARY} 
                      ${PROJECT_NAME}_core ${PROJECT_NAME}_physics
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

install(FILES 
  ${HEADERS_DECAY} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/decay" 
  COMPONENT dev)
install(FILES 
  ${HEADERS_DETECTOR} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/detector" 
  COMPONENT dev)
install(FILES 
  ${HEADERS_RC} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/reconstruction" 
  COMPONENT dev)
