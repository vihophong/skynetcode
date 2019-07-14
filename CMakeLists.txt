# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

enable_testing()

# set Release build type as default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Configuration type (one of Debug, RelWithDebInfo, Release, MinSizeRel)"
    FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")


project(skynetuser1 CXX C Fortran)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# enable warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# enable C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# generate debugging symbols for release and debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")

# explicitly set DEBUG flag in Debug mode
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# to make shared libraries we need position independent code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

# linking the shared Intel libraries suppresses warnings about feupdateenv
# not being implemented, and to be sure that we're not missing any math functions
# we also link against libm
if(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") OR ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel"))
  set(SKYNET_EXTERNAL_LIBS "imf;m;${SKYNET_EXTERNAL_LIBS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -shared-intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -shared-intel")
endif ()

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{SKYNETSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
#find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Spectrum Net)
# find skynet
find_package(SkyNet REQUIRED)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Spectrum Net)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

find_package(Threads REQUIRED)
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

# we link MKL statically so that it will work from within Python
link_directories(${MKL_LIBRARY_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSKYNET_USE_MKL")


include_directories(${CMAKE_SOURCE_DIR} ${SKYNET_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

#---Create  a main program using the library
add_executable(r-process r-process.cpp)
target_link_libraries(r-process ${SKYNET_LIBRARIES})
target_link_libraries(r-process ${ROOT_LIBRARIES})   # Link2Dictionary!
