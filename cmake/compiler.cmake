################################################################################
## setup
################################################################################
include(${PROJECT_SOURCE_DIR}/cmake/rpath.cmake)


################################################################################
## global defines
################################################################################

## set the RPATH for OsX
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  init_rpath_macos()
  add_to_rpath_macos(${INSTALL_LIB_DIR})
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
  set (CXX_EXTRA_FLAGS "-fomit-frame-pointer -fpermissive")
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
find_package(ROOT COMPONENTS GenVector EG MathMore REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

## boost
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


## gsl
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_to_rpath_macos("${ROOT_LIBRARY_DIR}")
  add_to_rpath_macos("${Boost_LIBRARY_DIRS}")
  add_to_rpath_macos("${GSL_LIBDIR}")
endif()

#include (cmake/debug.cmake)
 
if( ${COMPILE_FOR_BROADWELL} ) 
 # http://www.prace-ri.eu/best-practice-guide-haswellbroadwell-january-2017/#id-1.4.2.6
 # From Table 13. found in the link above
 set(PROJECT_EXTRA_C_FLAGS       -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions)
 set(PROJECT_EXTRA_CXX_FLAGS     -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions)
 set(PROJECT_EXTRA_FORTRAN_FLAGS -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions)
 message("PROJECT_EXTRA_CXX_FLAGS : ${PROJECT_EXTRA_CXX_FLAGS}")
elseif( ${COMPILE_FOR_HASWELL} )
#same as broadwell
 set(PROJECT_EXTRA_C_FLAGS       -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions)
 set(PROJECT_EXTRA_CXX_FLAGS     -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions)
 set(PROJECT_EXTRA_FORTRAN_FLAGS -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions)

elseif( ${COMPILE_FOR_KNL} )
  # TODO 
  # http://www.prace-ri.eu/best-practice-guide-knights-landing-january-2017/
endif() 

