################################################################################
## CMAKE Settings
################################################################################
set (SUBDIR "gen")
set (LIBRARY "${PROJECT_NAME}_${SUBDIR}")
set (TARGETS ${TARGETS} ${LIBRARY} PARENT_SCOPE)

################################################################################
## Sources and install headers
################################################################################
file (GLOB SOURCES "*.cc" "initial/*.cc" "lA/*.cc")
file (GLOB HEADERS "*.hh")
file (GLOB INITIAL_HEADERS "initial/*.hh")
file (GLOB LA_HEADERS "lA/*.hh")

################################################################################
## Include directories
################################################################################
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

################################################################################
## Compile and Link
################################################################################
add_library(${LIBRARY} SHARED ${SOURCES})
target_link_libraries(${LIBRARY} 
                      ${PROJECT_NAME}_core ${PROJECT_NAME}_physics
                      ${PROJECT_NAME}_proc
                      ${ROOT_GenVector_LIBRARY} 
                      ${ROOT_EG_LIBRARY}
                      ${ROOT_MathMore_LIBRARY}
                      ${ROOT_Gsl_LIBRARY} ${ROOT_LIBRARIES}
                      ${GSL_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${JPACPHOTO_LIBRARIES})
target_compile_options(${LIBRARY} PUBLIC ${PROJECT_EXTRA_CXX_FLAGS})
set_target_properties(${LIBRARY} PROPERTIES 
  VERSION ${LAGER_VERSION} 
  SOVERSION ${LAGER_SOVERSION})

################################################################################
## Export and Install
################################################################################
install(TARGETS ${LIBRARY}
  EXPORT ${PROJECT_NAME}-targets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib)
  
install(FILES 
  ${HEADERS} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}" 
  COMPONENT dev)
install(FILES 
  ${INITIAL_HEADERS} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/initial" 
  COMPONENT dev)
install(FILES 
  ${LA_HEADERS} 
  DESTINATION "${INSTALL_INCLUDE_DIR}/${SUBDIR}/lA" 
  COMPONENT dev)
