CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)
PROJECT(TopoMS CXX)

set(CMAKE_BUILD_TYPE Release)

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    message("found GNU compiler")
    set(CMAKE_CXX_FLAGS "-O3")          ## Optimize
    #set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()

SET(CMAKE_CXX_FLAGS "-std=c++0x -fpermissive")

# ------------------------------------------------------------------------------
FIND_PACKAGE(OpenMP REQUIRED )

if (OPENMP_FOUND)
    message(STATUS "found OMP")
    message( STATUS "OpenMP_C_FLAGS: " ${OpenMP_C_FLAGS} )
    message( STATUS "OpenMP_CXX_FLAGS: " ${OpenMP_CXX_FLAGS} )
    message( STATUS "CMAKE_EXE_LINKER_FLAGS: " ${CMAKE_EXE_LINKER_FLAGS} )

    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fpermissive")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fpermissive")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# ------------------------------------------------------------------------------
FIND_PACKAGE(VTK)

IF(VTK_FOUND)
  message(STATUS "found VTK. Version:" ${VTK_VERSION} ". VTK_DIR: " ${VTK_DIR})
  include(${VTK_USE_FILE})
  add_definitions (-DUSE_VTK=TRUE)
ENDIF(VTK_FOUND)

# ------------------------------------------------------------------------------
add_definitions (-DUSE_KAHAN_SUM=TRUE)

# ------------------------------------------------------------------------------

SET(MSC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../msc)
SET(TOPOMS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

message( STATUS "MSC_DIR: " ${MSC_DIR} )
message( STATUS "TOPOMS_DIR: " ${TOPOMS_DIR} )

INCLUDE_DIRECTORIES(${MSC_DIR}/include ${TOPOMS_DIR}/include)

FILE(GLOB HEADERS1 ${MSC_DIR}/include/*.h)
FILE(GLOB SOURCES1 ${MSC_DIR}/src/*.cxx)

FILE(GLOB HEADERS2 ${TOPOMS_DIR}/include/*.h)
FILE(GLOB SOURCES2 ${TOPOMS_DIR}/src/*.cpp)

ADD_EXECUTABLE(TopoMS ${HEADERS1} ${HEADERS2} ${SOURCES1} ${SOURCES2} ${SOURCES})
IF(VTK_FOUND)
target_link_libraries(TopoMS LINK_PUBLIC vtkCommonCore vtkCommonDataModel vtkIOCore vtkIOXML vtkIOLegacy vtkIOMPIParallel)
endif(VTK_FOUND)
