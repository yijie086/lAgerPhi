################################################################################
## External Libraries
################################################################################

## LCIO2 (optional)
find_package(podio QUIET)
find_package(lcio2 QUIET)
if(podio_FOUND AND lcio2_FOUND)
  message("LCIO2/PODIO found")
  set(have_lcio2 1)
else()
  message("No LCIO2/PODIO found, not building ANL routines")
  set(have_lcio2 0)
endif()

################################################################################
## Copy over python library components
################################################################################
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pcsim 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

################################################################################
## Generate location-specific library components
################################################################################
configure_file(pcsim-jupyter.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim-jupyter)
configure_file(rootlogon.C.in
    ${CMAKE_CURRENT_BINARY_DIR}/rootlogon.C)
configure_file(pcsim/util/setup_root.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/setup_root.py)
configure_file(pcsim/util/setup_root.cc.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/setup_root.cc)
configure_file(pcsim/util/path.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/path.py)

################################################################################
## Cleanup and compile
################################################################################
add_custom_target(jupyter_cleanup ALL
  COMMAND find pcsim -name "*.in" -delete
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(jupyter_bindings ALL
  COMMAND python -c "import pcsim"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Building python/ROOT bindings")
