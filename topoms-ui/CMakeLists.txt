# ------------------------------------------------------------------------------
# CMake file for TopoMS-UI
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.9)

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
# Glew
# TODO: move this into a FindQGLViewer.cmake file
OPTION(USE_GLEW "Use GLEW for opengl wrangling (Windows)" OFF)
if (USE_GLEW)
	find_package(GLEW REQUIRED)
	include_directories(${GLEW_INCLUDE_DIRS})
	add_definitions(-DUSE_GLEW)
endif(USE_GLEW)
	

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
          HINTS ${PATH_EXT}/include
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
             HINTS ${PATH_EXT}/lib
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
# ui code

include_directories(${PATH_MSC}/include ${PATH_TOPOMS}/include)
include_directories(${PATH_UI}/include)

file(GLOB HDRS_TOPOMSUI ${PATH_UI}/include/*.h)
file(GLOB SRCS_TOPOMSUI ${PATH_UI}/src/*.cpp ${PATH_UI}/src/*.cxx)

# ------------------------------------------------------------------------------
# create the ui application

add_executable(${TRG_UI} ${HDRS_TOPOMSUI} ${SRCS_TOPOMSUI})

target_link_libraries(${TRG_UI} LINK_PUBLIC ${TRG_LIB})
target_link_libraries(${TRG_UI} LINK_PUBLIC
                        ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${QGLVIEWER_LIBRARY} )
if(USE_GLEW)
	target_link_libraries(${TRG_UI} LINK_PUBLIC ${GLEW_LIBRARIES})
endif(USE_GLEW)

install(TARGETS ${TRG_UI}
        BUNDLE DESTINATION bin
        RUNTIME DESTINATION bin)

# ------------------------------------------------------------------------------

