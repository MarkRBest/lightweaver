cmake_minimum_required (VERSION 2.8)
project (LightWeaver)

add_library (CODE vec3f.cpp Utils.cpp mesh.cpp camera.cpp ray.cpp sceen.cpp triangle.cpp light.cpp Octree.cpp texture.cpp bsptree.cpp)

FIND_PACKAGE( Boost COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)
#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)


#include_directories(library ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
endif()

add_executable(LightWeaver LightWeaver.cpp)
TARGET_LINK_LIBRARIES( LightWeaver CODE jpeg GL ${GLUT_LIBRARY} ${Boost_LIBRARIES})

#add_executable(Test Test.cpp)
#TARGET_LINK_LIBRARIES(Test CODE ${Boost_LIBRARIES})
