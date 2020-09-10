cmake_minimum_required(VERSION 2.6)

project(sep C)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

set(SOURCES
   ${CMAKE_SOURCE_DIR}/src/analyse.c
   ${CMAKE_SOURCE_DIR}/src/convolve.c
   ${CMAKE_SOURCE_DIR}/src/deblend.c
   ${CMAKE_SOURCE_DIR}/src/extract.c
   ${CMAKE_SOURCE_DIR}/src/lutz.c
   ${CMAKE_SOURCE_DIR}/src/aperture.c
   ${CMAKE_SOURCE_DIR}/src/background.c
   ${CMAKE_SOURCE_DIR}/src/util.c
   )

include_directories(${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/dsp/ ${CMAKE_SOURCE_DIR}/sep/src ${CFITSIO_INCLUDE_DIR})

link_directories(${CMAKE_LIBRARY_PATH} ${CMAKE_SOURCE_DIR}/sep/src)

add_library(sep SHARED ${SOURCES})
set_target_properties(sep PROPERTIES OUTPUT_NAME sep)

target_link_libraries(sep m)

install(TARGETS sep LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
