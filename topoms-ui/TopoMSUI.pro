TARGET = TopoMS-UI
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.12

CONFIG *= qt opengl thread release warn_off
QT *=  gui widgets opengl xml printsupport

QMAKE_CXXFLAGS += -fpermissive -w -fopenmp -std=c++11
LIBS += -fopenmp
CONFIG += c++11

macx {
    CONFIG -= app_bundle
    INCLUDEPATH *= /usr/include     # for Mac OS 10.13, this additional include is needed
}

## -------------------------------------------------
## get the cpp compiler
## -------------------------------------------------
QMAKE_LINK              = $$QMAKE_CXX
QMAKE_LINK_SHLIB        = $$QMAKE_CXX
message('QMAKE_CXX  : '$$QMAKE_CXX)

## -------------------------------------------------
## paths to temporary output files
## -------------------------------------------------

DESTDIR = $$OUT_PWD
MOC_DIR = $$OUT_PWD/obj/
OBJECTS_DIR = $$OUT_PWD/obj/
UI_DIR = $$OUT_PWD/obj/
RCC_DIR = $$OUT_PWD/obj/

## -------------------------------------------------
## specify external libraries and paths to be used
## -------------------------------------------------

isEmpty(VTKVERSION){   VTKVERSION = 7.1   }
isEmpty(VTKPATH){      VTKPATH = $$(HOME)/macports   }
isEmpty(QGLPATH){      QGLPATH = $$(HOME)/usr   }
isEmpty(GLEWPATH){     GLEWPATH = $$(HOME)/macports   }

message('QGL_PATH   : '$$QGLPATH)
message('VTK_PATH   : '$$VTKPATH)
message('VTK_VERSION: '$$VTKVERSION)

# QGLViewer
INCLUDEPATH *= $$QGLPATH/include
LIBS *= -L$$QGLPATH/lib
LIBS *= -lQGLViewer

# GLEW
DEFINES += USE_GLEW
INCLUDEPATH *= $$GLEWPATH/include
LIBS *= -L$$GLEWPATH/lib
LIBS *= -lGLEW

# vtk
DEFINES += USE_VTK
INCLUDEPATH += $$VTKPATH/include/vtk-$$VTKVERSION
LIBS *= -L$$VTKPATH/lib
LIBS += \
-lvtkCommonCore-$$VTKVERSION \
-lvtkCommonDataModel-$$VTKVERSION \
-lvtkIOCore-$$VTKVERSION \
-lvtkIOXML-$$VTKVERSION \
-lvtkIOLegacy-$$VTKVERSION

# Kahan sum for more accurate summation
DEFINES += USE_KAHAN_SUM

## -------------------------------------------------
## all headers and sources for this project
## -------------------------------------------------

MSC_PATH = $$PWD/../msc
TOPOMS_PATH = $$PWD/../topoms

FORMS *= $$PWD/include/TopoMSUI.ui

INCLUDEPATH *= $$MSC_PATH/include \
               $$TOPOMS_PATH/include \
               $$PWD/include

HEADERS = $$MSC_PATH/include/*.h \
          $$TOPOMS_PATH/include/*.h \
          $$PWD/include/*.h

SOURCES = $$MSC_PATH/src/*.cxx \
          $$TOPOMS_PATH/src/TopoMS.cpp \
          $$TOPOMS_PATH/src/Utils.cpp  \
          $$PWD/src/*.cpp

## -------------------------------------------------
## -------------------------------------------------
