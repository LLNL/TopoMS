TARGET = TopoMS-UI
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.11

CONFIG *= qt release thread warn_off
QT *= xml opengl printsupport widgets

QMAKE_CXXFLAGS += -fpermissive -fopenmp -std=c++11
LIBS += -fopenmp
CONFIG += c++11

macx {
    CONFIG -= app_bundle
}

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

QGLPATH = $$(HOME)/usr
VTKPATH = $$(HOME)/usr
VTKVERSION = 7.0



# QGLViewer
INCLUDEPATH *= $$QGLPATH/include
LIBS *= $$QGLPATH/lib/libQGLViewer.a

# vtk
INCLUDEPATH += $$VTKPATH/include/vtk-$$VTKVERSION
LIBS *= -L$$VTKPATH/lib
DEFINES += USE_VTK

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

FORMS *= $$PWD/TopoMSUI.ui

INCLUDEPATH *= $$MSC_PATH/include \
               $$TOPOMS_PATH/include \
               $$PWD/include

HEADERS = $$MSC_PATH/include/*.h \
          $$TOPOMS_PATH/include/*.h \
          $$PWD/include/*.h

SOURCES = $$MSC_PATH/src/*.cxx \
          $$TOPOMS_PATH/src/TopoMS.cpp $$TOPOMS_PATH/src/Utils.cpp  \
          $$PWD/src/*.cpp

## -------------------------------------------------
## -------------------------------------------------
