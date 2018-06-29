CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)
PROJECT(TopoMS-UI CXX)

IF(NOT CMAKE_COMPILER_IS_GNUCXX)
  MESSAGE(FATAL_ERROR "TopoMS requires a GNU compiler. CMake will fail.")
  exit()
ENDIF()

set(CMAKE_BUILD_TYPE Release)

## Compiler flags
set(CMAKE_CXX_FLAGS "-O3")          ## Optimize
set(CMAKE_CXX_FLAGS "-std=c++0x -fpermissive -w")
#set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

# other project definitions
ADD_DEFINITIONS(-DUSE_KAHAN_SUM=TRUE)

# ------------------------------------------------------------------------------
FIND_PACKAGE(OpenMP REQUIRED)

MESSAGE( STATUS "OpenMP_C_FLAGS: " ${OpenMP_C_FLAGS} )
MESSAGE( STATUS "OpenMP_CXX_FLAGS: " ${OpenMP_CXX_FLAGS} )
MESSAGE( STATUS "CMAKE_EXE_LINKER_FLAGS: " ${CMAKE_EXE_LINKER_FLAGS} )

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fpermissive")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fpermissive")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# ------------------------------------------------------------------------------
FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  MESSAGE(STATUS "Found VTK " ${VTK_VERSION})
  INCLUDE(${VTK_USE_FILE})
  ADD_DEFINITIONS(-DUSE_VTK=TRUE)
ENDIF(VTK_FOUND)

# ------------------------------------------------------------------------------
find_package(OpenGL REQUIRED)

find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets Xml PrintSupport REQUIRED)
set(QT_LIBRARIES Qt5::Core Qt5::Gui Qt5::OpenGL Qt5::Widgets Qt5::Xml Qt5::PrintSupport)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)

include_directories( ${Qt5Core_INCLUDE_DIRS}
                     ${Qt5Gui_INCLUDE_DIRS}
                     ${Qt5OpenGL_INCLUDE_DIRS}
                     ${Qt5Widgets_INCLUDE_DIRS}
                     ${Qt5Xml_INCLUDE_DIRS}
                     ${Qt5PrintSupport_INCLUDE_DIRS} )

MESSAGE( STATUS "QGLPATH: " ${QGLPATH} )
INCLUDE_DIRECTORIES(${QGLPATH}/include)
set(QGLVIEWER_LIBRARIES ${QGLPATH}/lib/libQGLViewer.a)

# ------------------------------------------------------------------------------
SET(MSC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../msc)
SET(TOPOMS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../topoms)
SET(TOPOMSUI_DIR ${CMAKE_CURRENT_SOURCE_DIR})

MESSAGE( STATUS "MSC_DIR: " ${MSC_DIR} )
MESSAGE( STATUS "TOPOMS_DIR: " ${TOPOMS_DIR} )
MESSAGE( STATUS "TOPOMS_DIR: " ${TOPOMSUI_DIR} )

INCLUDE_DIRECTORIES(${MSC_DIR}/include ${TOPOMS_DIR}/include ${TOPOMSUI_DIR}/include)

FILE(GLOB HEADERS_MSC ${MSC_DIR}/include/*.h)
FILE(GLOB SOURCES_MSC ${MSC_DIR}/src/*.cxx)

FILE(GLOB HEADERS_TOPOMS ${TOPOMS_DIR}/include/*.h)
FILE(GLOB SOURCES_TOPOMS ${TOPOMS_DIR}/src/TopoMS.cpp ${TOPOMS_DIR}/src/Utils.cpp)

FILE(GLOB HEADERS_TOPOMSUI ${TOPOMSUI_DIR}/include/*.h)
FILE(GLOB SOURCES_TOPOMSUI ${TOPOMSUI_DIR}/src/*.cpp)

# ------------------------------------------------------------------------------

ADD_EXECUTABLE(TopoMS-UI ${HEADERS_MSC} ${HEADERS_TOPOMS} ${HEADERS_TOPOMSUI}
                         ${SOURCES_MSC} ${SOURCES_TOPOMS} ${SOURCES_TOPOMSUI})

TARGET_LINK_LIBRARIES(TopoMS-UI LINK_PUBLIC
                      vtkCommonCore vtkCommonDataModel vtkIOCore vtkIOXML vtkIOLegacy
                      ${QT_LIBRARIES} ${OPENGL_LIBRARIES} ${QGLVIEWER_LIBRARIES})
# ------------------------------------------------------------------------------
