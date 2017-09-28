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
configure_file(pcsim-python.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim-python)
configure_file(pcsim-ipython.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim-ipython)
configure_file(pcsim/util/setup_root.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/setup_root.py)
configure_file(pcsim/util/path.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/path.py)

################################################################################
## Cleanup and compile
################################################################################
add_custom_target(python_cleanup ALL
  COMMAND find pcsim -name "*.in" -delete
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(python_bindings ALL
  COMMAND python -c "import pcsim"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Building python/ROOT bindings")
