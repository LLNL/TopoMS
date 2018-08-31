# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6)

enable_language(CXX)
project(TopoMS CXX)
set(CMAKE_BUILD_TYPE Release)

# if(NOT CMAKE_COMPILER_IS_GNUCXX)
  # message(FATAL_ERROR "TopoMS requires a GNU compiler. CMake will fail.")
# endif()

## Compiler flags
set(CMAKE_CXX_FLAGS         "-O3")                         ## Optimize
set(CMAKE_CXX_FLAGS         "-std=c++0x -fpermissive -w")
#set(CMAKE_EXE_LINKER_FLAGS "-s")                          ## Strip binary

# other project definitions
add_definitions(-DUSE_KAHAN_SUM=TRUE)

# ------------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

message(STATUS "OpenMP_CXX_FLAGS: "        ${OpenMP_CXX_FLAGS})
message(STATUS "OpenMP_EXE_LINKER_FLAGS: " ${OpenMP_EXE_LINKER_FLAGS})
message(STATUS "CMAKE_EXE_LINKER_FLAGS: "  ${CMAKE_EXE_LINKER_FLAGS})

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}         ${OpenMP_CXX_FLAGS} -fpermissive")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}  ${OpenMP_EXE_LINKER_FLAGS}")

# ------------------------------------------------------------------------------
find_package(VTK 7.1
                 HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../external/lib/cmake/vtk-7.1)
if(VTK_FOUND)
  add_definitions(-DUSE_VTK=TRUE)
  message(STATUS "Found VTK " ${VTK_VERSION} ": " ${VTK_INCLUDE_DIRS})
  include(${VTK_USE_FILE})
  #set(VTK_LIBRARIES vtkCommonCore vtkCommonDataModel vtkIOCore vtkIOXML vtkIOLegacy)
endif(VTK_FOUND)

# ------------------------------------------------------------------------------
set(MSC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../msc)
set(TOPOMS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

message(STATUS "MSC_DIR: "    ${MSC_DIR} )
message(STATUS "TopoMS_DIR: " ${TOPOMS_DIR} )

include_directories(${MSC_DIR}/include ${TOPOMS_DIR}/include)

file(GLOB HDRS_MSC ${MSC_DIR}/include/*.h)
file(GLOB SRCS_MSC ${MSC_DIR}/src/*.cxx)

file(GLOB HDRS_TOPOMS ${TOPOMS_DIR}/include/*.h)
file(GLOB SRCS_TOPOMS ${TOPOMS_DIR}/src/*.cpp)

# ------------------------------------------------------------------------------
add_executable(TopoMS ${HDRS_MSC} ${HDRS_TOPOMS} ${SRCS_MSC} ${SRCS_TOPOMS})

if(VTK_FOUND)
  target_link_libraries(TopoMS LINK_PUBLIC ${VTK_LIBRARIES})
endif(VTK_FOUND)
# ------------------------------------------------------------------------------
