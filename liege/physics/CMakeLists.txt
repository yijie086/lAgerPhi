################################################################################
## CMAKE Settings
################################################################################
set (SUBDIR "physics")
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
## Compile and Link
################################################################################
add_library(${LIBRARY} SHARED ${SOURCES})
target_compile_options(${LIBRARY} PUBLIC ${PROJECT_EXTRA_CXX_FLAGS})
target_link_libraries(${LIBRARY} 
                      ${PROJECT_NAME}_core 
                      ${ROOT_GenVector_LIBRARY} ${ROOT_EG_LIBRARY} ${ROOT_LIBRARIES} 
                      ${Boost_LIBRARIES})
set_target_properties(${LIBRARY} PROPERTIES 
  VERSION ${LIEGE_VERSION} 
  SOVERSION ${LIEGE_SOVERSION}
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
