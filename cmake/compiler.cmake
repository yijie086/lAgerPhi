################################################################################
## global defines
################################################################################
## NOTHING HERE

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
  set (CXX_EXTRA_FLAGS "-fomit-frame-pointer -std=c++14")
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
