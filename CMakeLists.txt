# ------------------------------------------------------------------------------
# CMake file for TopoMS
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.9)

message(STATUS "> configuring TopoMS in (" ${CMAKE_CURRENT_SOURCE_DIR} ")")

project(TopoMS LANGUAGES CXX VERSION 1.1 DESCRIPTION "TopoMS")
set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_VERBOSE_MAKEFILE 1)

# ------------------------------------------------------------------------------
# compiler configuration

#if(NOT CMAKE_COMPILER_IS_GNUCXX)
#  message(FATAL_ERROR "TopoMS requires a GNU compiler. CMake will fail.")
#endif()

# Compiler flags
set(CMAKE_CXX_FLAGS         "-O3")                         ## Optimize
set(CMAKE_CXX_FLAGS         "-std=c++0x -fpermissive -w")

set(CMAKE_SHARED_LINKER_FLAGS "-s -lc++")                     ## Strip binary
set(CMAKE_EXE_LINKER_FLAGS    "-s -lc++")                     ## Strip binary
    #-lc++ for Apple

# ------------------------------------------------------------------------------
# subdirectories

set(PATH_EXT    ${CMAKE_CURRENT_SOURCE_DIR}/external)
set(PATH_MSC    ${CMAKE_CURRENT_SOURCE_DIR}/msc)
set(PATH_TOPOMS ${CMAKE_CURRENT_SOURCE_DIR}/topoms)
set(PATH_CLI    ${CMAKE_CURRENT_SOURCE_DIR}/topoms-cli)
set(PATH_UI     ${CMAKE_CURRENT_SOURCE_DIR}/topoms-ui)

#message(STATUS " > PATH_MSC: "    ${PATH_MSC})
#message(STATUS " > PATH_TOPOMS: " ${PATH_TOPOMS})
#message(STATUS " > PATH_CLI: "    ${PATH_CLI})
#message(STATUS " > PATH_UI: "     ${PATH_UI})

# target names
set(TRG_LIB topoms)
set(TRG_CLI topoms-cli)
set(TRG_UI  topoms-ui)

# ------------------------------------------------------------------------------
# OpenMP

find_package(OpenMP REQUIRED)

#message(STATUS "OpenMP_CXX_FLAGS: " ${OpenMP_CXX_FLAGS})
#message(STATUS "OpenMP_EXE_LINKER_FLAGS: " ${OpenMP_EXE_LINKER_FLAGS})

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}         ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}  ${OpenMP_EXE_LINKER_FLAGS}")

# ------------------------------------------------------------------------------
# VTK

find_package(VTK 7.1
             HINTS ${PATH_EXT}/lib/cmake/vtk-7.1)

if(VTK_FOUND)
  message(STATUS "Found VTK " ${VTK_VERSION} ": " ${VTK_INCLUDE_DIRS})
  include(${VTK_USE_FILE})
  add_definitions(-DUSE_VTK=TRUE)
endif(VTK_FOUND)

# ------------------------------------------------------------------------------
# TopoMS

add_subdirectory(topoms)
add_subdirectory(topoms-cli)

OPTION(TOPOMS_BUILD_UI "Build User Interface" OFF)
if(TOPOMS_BUILD_UI)
    add_subdirectory(topoms-ui)
endif()

# ------------------------------------------------------------------------------
