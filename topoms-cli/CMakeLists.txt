# ------------------------------------------------------------------------------
# CMake file for topoms-cli
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.9)

include_directories(${PATH_MSC}/include ${PATH_TOPOMS}/include)

add_executable(${TRG_CLI} ${PATH_CLI}/src/main.cpp)
target_link_libraries(${TRG_CLI} LINK_PUBLIC ${TRG_LIB})

install(TARGETS ${TRG_CLI} DESTINATION bin)

# ------------------------------------------------------------------------------
