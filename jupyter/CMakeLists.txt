file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pcsim 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


################################################################################
## External Libraries
################################################################################
## ROOT
find_package(ROOT REQUIRED)
## boost
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)

configure_file(pcsim/util/setup_root.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/setup_root.py)
configure_file(pcsim/util/path.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/pcsim/util/path.py)
