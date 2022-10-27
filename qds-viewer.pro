# -*- mode: Makefile -*-

TARGET = qds-viewer
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml
equals (QT_MAJOR_VERSION, 6) {
    QT += openglwidgets
}

INCLUDEPATH += ../libgeom

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

QMAKE_CXXFLAGS += -O3

unix:LIBS *= -lQGLViewer-qt5 -lOpenMeshCore -L../libgeom/release -lgeom -lGL -lGLU

RESOURCES = qds-viewer.qrc
