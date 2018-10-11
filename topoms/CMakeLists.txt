# ------------------------------------------------------------------------------
# CMake file for topoms (shared library)
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.9)

# ------------------------------------------------------------------------------
# TopoMS source

include_directories(${PATH_MSC}/include ${PATH_TOPOMS}/include)

file(GLOB HDRS_MSC ${PATH_MSC}/include/*.h)
file(GLOB SRCS_MSC ${PATH_MSC}/src/*.cpp)

file(GLOB HDRS_TOPOMS ${PATH_TOPOMS}/include/*.h)
file(GLOB SRCS_TOPOMS ${PATH_TOPOMS}/src/*.cpp)

# ------------------------------------------------------------------------------
# create and install target

add_library(${TRG_LIB} STATIC
            ${HDRS_MSC} ${HDRS_TOPOMS} ${SRCS_MSC} ${SRCS_TOPOMS})

if(VTK_FOUND)
  target_link_libraries(${TRG_LIB} LINK_PUBLIC ${VTK_LIBRARIES})
endif(VTK_FOUND)

install(TARGETS ${TRG_LIB} DESTINATION lib )

# ------------------------------------------------------------------------------
