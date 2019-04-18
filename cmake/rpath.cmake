################################################################################
## CMAKE Settings
################################################################################

## function to set the RPATH settings for osx
function (init_rpath_macos)
  set(CMAKE_MACOSX_RPATH 1 PARENT_SCOPE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE PARENT_SCOPE)
  #set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE PARENT_SCOPE)
endfunction ()

## function to add a library directory to the RPATH if needed
function (add_to_rpath LIBDIR)
  ## only add if not a system dir
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIBDIR}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}" "${LIBDIR}" PARENT_SCOPE)
    message(STATUS "Added ${LIBDIR} to the CMAKE_INSTALL_RPATH for MacOS")
  else ()
    message(STATUS "${LIBDIR} is a system directory, not added to CMAKE_INSTALL_RPATH")
  endif()
endfunction ()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  init_rpath_macos()
  add_to_rpath(${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR})
endif()
