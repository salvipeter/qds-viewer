cmake_minimum_required(VERSION 3.20)
project(qds-viewer)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(OpenGL_GL_PREFERENCE LEGACY)

find_package(OpenGL REQUIRED)
find_package(OpenMesh REQUIRED)
find_package(Qt6 REQUIRED COMPONENTS Core OpenGL OpenGLWidgets Widgets Xml)
find_package(QGLViewer REQUIRED)
find_package(LibGeom REQUIRED)
find_package(Triangle REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -fPIC -Wno-reorder")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG -fsanitize=address")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DNDEBUG")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(qds_HEADERS # headers processed by MOC
  MyViewer.h
  MyWindow.h
  )

set(CMAKE_AUTORCC ON)
qt6_wrap_cpp(qds_HEADERS_MOC ${qds_HEADERS})

add_definitions(${Qt6Widgets_DEFINITIONS} -D_USE_MATH_DEFINES)

include_directories(
  .
  ${LIBGEOM_INCLUDE_DIR}
  ${TRIANGLE_DIR}
  ${OPENGL_INCLUDE_DIR}
  ${OPENMESH_INCLUDE_DIR}
  ${Qt6OpenGL_INCLUDE_DIRS}
  ${Qt6Widgets_INCLUDE_DIRS}
  ${Qt6Xml_INCLUDE_DIRS}
  ${QGLVIEWER_INCLUDE_DIR}
  )

add_executable(qds-viewer WIN32
  ${qds_HEADERS_MOC}
  main.cpp
  qds-viewer.qrc
  MyViewer.cpp
  MyWindow.cpp
  ${TRIANGLE_DIR}/triangle.o
  )

if(UNIX)
  add_custom_command(OUTPUT tags
    COMMAND etags *.h *.cpp
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  target_sources(qds-viewer PRIVATE tags)
endif(UNIX)

target_link_libraries(qds-viewer PUBLIC
  ${LIBGEOM_LIBRARY}
  ${OPENGL_LIBRARIES}
  ${OPENMESH_LIBRARIES}
  Qt6::Core
  Qt6::OpenGL
  Qt6::OpenGLWidgets
  Qt6::Widgets
  Qt6::Xml
  ${QGLVIEWER_LIBRARY}
  )
