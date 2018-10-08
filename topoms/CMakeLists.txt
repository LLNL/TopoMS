# ------------------------------------------------------------------------------
# CMake file for TopoMS
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.6)

enable_language(CXX)
project(TopoMS CXX)

set(CMAKE_BUILD_TYPE Release)

# ------------------------------------------------------------------------------
# compiler configuration

if(NOT CMAKE_COMPILER_IS_GNUCXX)
  message(FATAL_ERROR "TopoMS requires a GNU compiler. CMake will fail.")
endif()

## Compiler flags
set(CMAKE_CXX_FLAGS         "-O3")                         ## Optimize
set(CMAKE_CXX_FLAGS         "-std=c++0x -fpermissive -w")
#set(CMAKE_EXE_LINKER_FLAGS "-s")                          ## Strip binary

# ------------------------------------------------------------------------------
# OpenMP

find_package(OpenMP REQUIRED)

message(STATUS "OpenMP_CXX_FLAGS: " ${OpenMP_CXX_FLAGS})
message(STATUS "OpenMP_EXE_LINKER_FLAGS: " ${OpenMP_EXE_LINKER_FLAGS})

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}         ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}  ${OpenMP_EXE_LINKER_FLAGS}")

# ------------------------------------------------------------------------------
# VTK

find_package(VTK 7.1
             HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../external/lib/cmake/vtk-7.1)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
  add_definitions(-DUSE_VTK=TRUE)
  message(STATUS "Found VTK " ${VTK_VERSION} ": " ${VTK_INCLUDE_DIRS})
  message(STATUS "Found VTK " ${VTK_VERSION} ": " ${VTK_LIBRARIES})
endif(VTK_FOUND)

# ------------------------------------------------------------------------------
# TopoMS Code

set(MSC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../msc)
set(TOPOMS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

message(STATUS "MSC_DIR: "    ${MSC_DIR} )
message(STATUS "TopoMS_DIR: " ${TOPOMS_DIR} )

include_directories(${MSC_DIR}/include ${TOPOMS_DIR}/include)

file(GLOB HDRS_MSC ${MSC_DIR}/include/*.h)
file(GLOB SRCS_MSC ${MSC_DIR}/src/*.cpp)

file(GLOB HDRS_TOPOMS ${TOPOMS_DIR}/include/*.h)
file(GLOB SRCS_TOPOMS ${TOPOMS_DIR}/src/*.cpp ${TOPOMS_DIR}/src/*.cxx)

# ------------------------------------------------------------------------------
# create executable

add_executable(topoms ${HDRS_MSC} ${HDRS_TOPOMS} ${SRCS_MSC} ${SRCS_TOPOMS})

if(VTK_FOUND)
  target_link_libraries(topoms LINK_PUBLIC ${VTK_LIBRARIES})
endif(VTK_FOUND)

# ------------------------------------------------------------------------------
