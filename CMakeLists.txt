cmake_minimum_required(VERSION 2.8.4)
project(advection)

set(SOURCE_FILES cfd_advection.cpp cfd.h cfd.cpp cfdUtility.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories("/usr/include")
find_library(OIIO "OpenImageIO" "/usr/lib")
find_library(GLUT "glut")
find_library(GL "GL")
find_library(GLU "GLU")

add_executable(cfd_advection ${SOURCE_FILES})

target_link_libraries(cfd_advection ${OIIO} ${GLUT} ${GL} ${GLU})
