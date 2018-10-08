# ------------------------------------------------------------------------------
# CMake file for TopoMS-UI
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.6)

enable_language(CXX)
project(TopoMS-UI CXX)

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

message(STATUS "OpenMP_CXX_FLAGS: "        ${OpenMP_CXX_FLAGS})
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
# OpenGL and Qt

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})

find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets Xml PrintSupport REQUIRED)
set(QT_LIBRARIES Qt5::Core Qt5::Gui Qt5::OpenGL Qt5::Widgets Qt5::Xml Qt5::PrintSupport)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)

include_directories(${Qt5Core_INCLUDE_DIRS}   ${Qt5Gui_INCLUDE_DIRS}
                    ${Qt5OpenGL_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS}
                    ${Qt5Xml_INCLUDE_DIRS}    ${Qt5PrintSupport_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# libQGLViewer
# TODO: move this into a FindQGLViewer.cmake file

find_path(QGLVIEWER_INCLUDE_DIR
          NAMES QGLViewer/qglviewer.h
          PATHS
                /usr/include
                /opt/local/include
                /usr/local/include
                /sw/include
                ${QGL_DIR}/include
          HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../external/include
          PATH_SUFFIXES Headers)

find_library(QGLVIEWER_LIBRARY
             NAMES qglviewer QGLViewer qglviewer-qt5 QGLViewer-qt5
             PATHS
                   /usr/lib
                   /usr/local/lib
                   /opt/local/lib
                   /sw/lib
                   ${QGL_DIR}/lib
                   ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
             HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../external/lib
             PATH_SUFFIXES QGLViewer QGLViewer/release)

if (QGLVIEWER_INCLUDE_DIR)
  message(STATUS "Found QGLViewer include directory: " ${QGLVIEWER_INCLUDE_DIR})
  include_directories(${QGLVIEWER_INCLUDE_DIR})
else()
  message(FATAL_ERROR "QGLViewer include directory not found.")
endif()

if (QGLVIEWER_LIBRARY)
  message(STATUS "Found QGLViewer library: " ${QGLVIEWER_LIBRARY})
else()
  message(FATAL_ERROR "QGLViewer library not found.")
endif()

# ------------------------------------------------------------------------------
# TopoMS Code

set(MSC_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/../msc)
set(TOPOMS_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/../topoms)
set(TOPOMSUI_DIR ${CMAKE_CURRENT_SOURCE_DIR})

message(STATUS "MSC_DIR: "       ${MSC_DIR})
message(STATUS "TopoMS_DIR: "    ${TOPOMS_DIR})
message(STATUS "TopoMS-UI_DIR: " ${TOPOMSUI_DIR})

include_directories(${MSC_DIR}/include ${TOPOMS_DIR}/include ${TOPOMSUI_DIR}/include)

file(GLOB HDRS_MSC ${MSC_DIR}/include/*.h)
file(GLOB SRCS_MSC ${MSC_DIR}/src/*.cpp)

file(GLOB HDRS_TOPOMS ${TOPOMS_DIR}/include/*.h)
file(GLOB SRCS_TOPOMS ${TOPOMS_DIR}/src/*.cpp)

file(GLOB HDRS_TOPOMSUI ${TOPOMSUI_DIR}/include/*.h)
file(GLOB SRCS_TOPOMSUI ${TOPOMSUI_DIR}/src/*.cpp ${TOPOMSUI_DIR}/src/*.cxx)

# ------------------------------------------------------------------------------
# create executable

add_executable(topoms-ui ${HDRS_MSC} ${HDRS_TOPOMS} ${HDRS_TOPOMSUI}
                         ${SRCS_MSC} ${SRCS_TOPOMS} ${SRCS_TOPOMSUI})

if(VTK_FOUND)
  target_link_libraries(topoms-ui LINK_PUBLIC ${VTK_LIBRARIES})
endif(VTK_FOUND)

target_link_libraries(topoms-ui LINK_PUBLIC
                        ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${QGLVIEWER_LIBRARY})

# ------------------------------------------------------------------------------
