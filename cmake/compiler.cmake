################################################################################
## global defines
################################################################################
## set the RPATH for OsX
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_MACOSX_RPATH 1)
endif()

################################################################################
## CXX Compiler Settings 
################################################################################
enable_language (CXX)

## set special compiler flags
get_filename_component(CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
## gcc and clang
if (CXX_COMPILER_NAME MATCHES "c\\+\\+.*" OR 
    CXX_COMPILER_NAME MATCHES "g\\+\\+.*" OR 
    CXX_COMPILER_NAME MATCHES "clang\\+\\+.*")
  set (CXX_EXTRA_FLAGS "-fomit-frame-pointer -std=c++14 -fpermissive")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${CXX_EXTRA_FLAGS}")
## add additional compilers here
## other compilers use the defaults:
else ()
  message ("CMAKE_CXX_COMPILER full path: " ${CMAKE_CXX_COMPILER})
  message ("C++ compiler: " ${CXX_COMPILER_NAME})
  message ("No optimized C++ compiler flags are known, using the defaults...")
  message ("Add the correct rules to cmake/compiler.cmake if other behavior is"
           "required.")
endif ()
## do not check for missing symbols when creating shared libraries with clang
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

################################################################################
## External Libraries
################################################################################
## ROOT
if(EXISTS $ENV{ROOTSYS}/cmake/ROOTConfig.cmake)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else()
  list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
endif()
find_package(ROOT COMPONENTS GenVector EG REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

## boost
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

## RPATH on MAcOS
## Set correct rpath on MacOs
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH
    ${CMAKE_INSTALL_RPATH} 
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
message("Setting RPath to ${CMAKE_INSTALL_RPATH}")
endif()
