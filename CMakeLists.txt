# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.10)

#sets the project name
PROJECT(myproject)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(INCLUDES
  MLRooDataStore.h
  MLWCSHandler.h  
  )

set(SOURCES
  src/MLRooDataStore.cxx
  src/MLWCSHandler.cpp
  )

#################
# Configuration #
#################


## ============================================================================
##  ROOT
## ============================================================================

# Look for ROOT
message("Looking for ROOT > 6.0...")
find_package(ROOT 6 REQUIRED COMPONENTS
  RooFitCore
  TreePlayer
  MathCore
  Gui
  pthread
  MathMore
  )
list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})
message(STATUS "Found ROOT!")
include(${ROOT_USE_FILE})
message(STATUS "ROOT include directories: ${ROOT_INCLUDE_DIRS}")
include_directories("${ROOT_INCLUDE_DIRS}")
message(STATUS "Using ROOT libraries ${ROOT_LIBRARIES}")

## ============================================================================
##  PKG-CONFIG
## ============================================================================

##Find PkgConfig which is required for VEGAS to work
find_package(PkgConfig REQUIRED)
if(PKG_CONFIG_FOUND)
  message(STATUS "PkgConfig found.")
else(NOT PKG_CONFIG_FOUND)
  message(FATAL_ERROR "PkgConfig not found")
endif(PKG_CONFIG_FOUND)


## ============================================================================
##  WCSTOOLS
## ============================================================================

# Look for libwcs
message("Looking for wcstools...")
find_path(WCSTOOLS_INCLUDE_DIRS NAMES fitshead.h PATHS $ENV{WCSTOOLS})
if(WCSTOOLS_INCLUDE_DIRS)
  message(STATUS "Using WCSTOOLS : ${WCSTOOLS_INCLUDE_DIRS}")
  include_directories(${WCSTOOLS_INCLUDE_DIRS})
  set(WCSTOOLS_LIBRARIES "${WCSTOOLS_INCLUDE_DIRS}/libwcs.a")
else(NOT WCSTOOLS_INCLUDE_DIRS)
  message(STATUS "Could not find a usable WCSTOOLS installation!")
endif(WCSTOOLS_INCLUDE_DIRS)

## ============================================================================
##  VEGAS
## ============================================================================


# Look for VEGAS
message("Looking for VEGAS...")
find_path(VEGAS_BASE_DIR NAMES resultsExtractor/include/VAStage6Analysis.h PATHS ${VEGAS} $ENV{VEGAS})
if(VEGAS_BASE_DIR)
  message(STATUS "Found VEGAS install at ${VEGAS_BASE_DIR}" )
  # add includes and libraries according to the old convention prior to 2.5.6
  set(VEGAS_INCLUDE_DIRS
      "${VEGAS_BASE_DIR}/include/common/include"
      "${VEGAS_BASE_DIR}/include/resultsExtractor/include"
      "${VEGAS_BASE_DIR}/include/coordinates/include"
      "${VEGAS_BASE_DIR}/include/showerReconstruction2/include"
  )
  set(VEGAS_LIBRARIES "${VEGAS_BASE_DIR}/lib/libVEGAScommon.so" "${VEGAS_BASE_DIR}/lib/libVEGAScoord.so")

  #add_library(VEGAS_SP24_LIB SHARED IMPORTED GLOBAL)
  #set_target_properties(VEGAS_SP24_LIB PROPERTIES IMPORTED_LOCATION "${VEGAS_BASE_DIR}/common/lib/libSP24sharedLite.so")
  #list(APPEND VEGAS_LIBRARIES ${VEGAS_SP24_LIB})

  #add_library(VEGAS_S6_LIB SHARED IMPORTED GLOBAL)
  #set_target_properties(VEGAS_S6_LIB PROPERTIES IMPORTED_LOCATION "${VEGAS_BASE_DIR}/resultsExtractor/lib/libStage6shared.so")
  #list(APPEND VEGAS_LIBRARIES ${VEGAS_S6_LIB})

  #add_library(VEGAS_COORD_LIB SHARED IMPORTED GLOBAL)
  #set_target_properties(VEGAS_COORD_LIB PROPERTIES IMPORTED_LOCATION "${VEGAS_BASE_DIR}/coordinates/lib/vaCoordShared.so")
  #list(APPEND VEGAS_LIBRARIES ${VEGAS_COORD_LIB})

  #add_library(VEGAS_S4_LIB SHARED IMPORTED GLOBAL)
  #set_target_properties(VEGAS_S4_LIB PROPERTIES IMPORTED_LOCATION "${VEGAS_BASE_DIR}/showerReconstruction2/lib/libStage4.so")
  #list(APPEND VEGAS_LIBRARIES ${VEGAS_S4_LIB})

  include_directories(${VEGAS_INCLUDE_DIRS})

else()
  message(FATAL_ERROR "Could not find a usable VEGAS installation! Exiting.")
endif(VEGAS_BASE_DIR)


message(STATUS "VEGAS include directories: ${VEGAS_INCLUDE_DIRS}")
#include_directories("${VEGAS_INCLUDE_DIRS}")
message(STATUS "Using VEGAS libraries ${VEGAS_LIBRARIES}")


include_directories(/work/LAS/amandajw-lab/users/achrmy/macros/RooDataSetClassLoopFilesCmakeITMonly/include)


add_definitions(
  -D_OAWG
  -D__STDC_LIMIT_MACROS
  -D__STDC_CONSTANT_MACROS
  -DNOHDF5
  -DENABLE_ML_ROOTFILE_IO
  -D_VEGASMLM
  -Wno-deprecated
  )

#add_executable(MLRooDataStore src/main.cpp)
add_executable(MLRooDataStore src/main.cpp ${SOURCES} ${ROOT_INCLUDE_DIRS} ${VEGAS_INCLUDE_DIRS} ${WCSTOOLS_INCLUDE_DIRS})
target_link_libraries(MLRooDataStore ${ROOT_LIBRARIES} ${VEGAS_LIBRARIES} ${WCSTOOLS_LIBRARIES})
 # ${SOURCES}
 # ${INCLUDES}
 # ${ROOT_LIBRARIES}
 # ${WCSTOOLS_LIBRARIES}
 # ${VEGAS_LIBRARIES}
 # )

#install(TARGETS MLRooDataStore DESTINATION ${CMAKE_INSTALL_BINDIR})
#install(MLRooDataStore ${CMAKE_INSTALL_BINDIR})

#requests that the executable be built
#add_executable(MLRooDataStore ${SOURCES})
#cmake ./ -DVEGAS_BASE_DIR=$VEGAS
#cmake --build .
#rm -r CMakeCache.txt CMakeFiles cmake_install.cmake
