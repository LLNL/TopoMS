# ------------------------------------------------------------------------------
# Qt pro file for TopoMS-UI
# ------------------------------------------------------------------------------

TARGET = TopoMS-UI
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.12

CONFIG *= qt opengl thread release warn_off
QT *=  gui widgets opengl xml printsupport

QMAKE_CXXFLAGS += -std=c++11 -fpermissive -w -fopenmp
LIBS += -fopenmp
CONFIG += c++11

macx {
    CONFIG -= app_bundle
    INCLUDEPATH *= /usr/include     # for Mac OS 10.13, this additional include is needed
}

# ------------------------------------------------------------------------------
# get the cpp compiler
# ------------------------------------------------------------------------------
QMAKE_LINK              = $$QMAKE_CXX
QMAKE_LINK_SHLIB        = $$QMAKE_CXX
message('QMAKE_CXX  : '$$QMAKE_CXX)

# ------------------------------------------------------------------------------
# paths to temporary output files
# ------------------------------------------------------------------------------

DESTDIR = $$OUT_PWD
MOC_DIR = $$OUT_PWD/obj/
OBJECTS_DIR = $$OUT_PWD/obj/
UI_DIR = $$OUT_PWD/obj/
RCC_DIR = $$OUT_PWD/obj/

# ------------------------------------------------------------------------------
# specify external libraries and paths to be used
# ------------------------------------------------------------------------------

isEmpty(VTKVERSION){   VTKVERSION = 7.1   }
isEmpty(VTKPATH){      VTKPATH = $$(HOME)/macports   }
isEmpty(QGLPATH){      QGLPATH = $$(HOME)/usr   }

message('QGL_PATH: '$$QGLPATH)
message('VTK_PATH: '$$VTKPATH)
message('VTK_VERSION: '$$VTKVERSION)

# QGLViewer
INCLUDEPATH *= $$QGLPATH/include
LIBS *= -L$$QGLPATH/lib
LIBS *= -lQGLViewer

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

LIBS += \
-lvtkChartsCore-$$VTKVERSION                           -lvtkIOPLYPython27D-$$VTKVERSION \
-lvtkChartsCorePython27D-$$VTKVERSION                  -lvtkIOParallel-$$VTKVERSION \
-lvtkCommonColor-$$VTKVERSION                          -lvtkIOParallelPython27D-$$VTKVERSION \
-lvtkCommonColorPython27D-$$VTKVERSION                 -lvtkIOParallelXML-$$VTKVERSION \
-lvtkCommonComputationalGeometry-$$VTKVERSION          -lvtkIOParallelXMLPython27D-$$VTKVERSION \
-lvtkCommonComputationalGeometryPython27D-$$VTKVERSION -lvtkIOSQL-$$VTKVERSION \
-lvtkCommonCore-$$VTKVERSION                           -lvtkIOSQLPython27D-$$VTKVERSION \
-lvtkCommonCorePython27D-$$VTKVERSION                  -lvtkIOTecplotTable-$$VTKVERSION \
-lvtkCommonDataModel-$$VTKVERSION                      -lvtkIOTecplotTablePython27D-$$VTKVERSION \
-lvtkCommonDataModelPython27D-$$VTKVERSION             -lvtkIOVideo-$$VTKVERSION \
-lvtkCommonExecutionModel-$$VTKVERSION                 -lvtkIOVideoPython27D-$$VTKVERSION \
-lvtkCommonExecutionModelPython27D-$$VTKVERSION        -lvtkIOXML-$$VTKVERSION \
-lvtkCommonMath-$$VTKVERSION                           -lvtkIOXMLParser-$$VTKVERSION \
-lvtkCommonMathPython27D-$$VTKVERSION                  -lvtkIOXMLParserPython27D-$$VTKVERSION \
-lvtkCommonMisc-$$VTKVERSION                           -lvtkIOXMLPython27D-$$VTKVERSION \
-lvtkCommonMiscPython27D-$$VTKVERSION                  -lvtkImagingColor-$$VTKVERSION \
-lvtkCommonSystem-$$VTKVERSION                         -lvtkImagingColorPython27D-$$VTKVERSION \
-lvtkCommonSystemPython27D-$$VTKVERSION                -lvtkImagingCore-$$VTKVERSION \
-lvtkCommonTransforms-$$VTKVERSION                     -lvtkImagingCorePython27D-$$VTKVERSION \
-lvtkCommonTransformsPython27D-$$VTKVERSION            -lvtkImagingFourier-$$VTKVERSION \
-lvtkDICOMParser-$$VTKVERSION                          -lvtkImagingFourierPython27D-$$VTKVERSION \
-lvtkDomainsChemistry-$$VTKVERSION                     -lvtkImagingGeneral-$$VTKVERSION \
-lvtkDomainsChemistryOpenGL2-$$VTKVERSION              -lvtkImagingGeneralPython27D-$$VTKVERSION \
-lvtkDomainsChemistryOpenGL2Python27D-$$VTKVERSION     -lvtkImagingHybrid-$$VTKVERSION \
-lvtkDomainsChemistryPython27D-$$VTKVERSION            -lvtkImagingHybridPython27D-$$VTKVERSION \
-lvtkFiltersAMR-$$VTKVERSION                           -lvtkImagingMath-$$VTKVERSION \
-lvtkFiltersAMRPython27D-$$VTKVERSION                  -lvtkImagingMathPython27D-$$VTKVERSION \
-lvtkFiltersCore-$$VTKVERSION                          -lvtkImagingMorphological-$$VTKVERSION \
-lvtkFiltersCorePython27D-$$VTKVERSION                 -lvtkImagingMorphologicalPython27D-$$VTKVERSION \
-lvtkFiltersExtraction-$$VTKVERSION                    -lvtkImagingSources-$$VTKVERSION \
-lvtkFiltersExtractionPython27D-$$VTKVERSION           -lvtkImagingSourcesPython27D-$$VTKVERSION \
-lvtkFiltersFlowPaths-$$VTKVERSION                     -lvtkImagingStatistics-$$VTKVERSION \
-lvtkFiltersFlowPathsPython27D-$$VTKVERSION            -lvtkImagingStatisticsPython27D-$$VTKVERSION \
-lvtkFiltersGeneral-$$VTKVERSION                       -lvtkImagingStencil-$$VTKVERSION \
-lvtkFiltersGeneralPython27D-$$VTKVERSION              -lvtkImagingStencilPython27D-$$VTKVERSION \
-lvtkFiltersGeneric-$$VTKVERSION                       -lvtkInfovisCore-$$VTKVERSION \
-lvtkFiltersGenericPython27D-$$VTKVERSION              -lvtkInfovisCorePython27D-$$VTKVERSION \
-lvtkFiltersGeometry-$$VTKVERSION                      -lvtkInfovisLayout-$$VTKVERSION \
-lvtkFiltersGeometryPython27D-$$VTKVERSION             -lvtkInfovisLayoutPython27D-$$VTKVERSION \
-lvtkFiltersHybrid-$$VTKVERSION                        -lvtkInteractionImage-$$VTKVERSION \
-lvtkFiltersHybridPython27D-$$VTKVERSION               -lvtkInteractionImagePython27D-$$VTKVERSION \
-lvtkFiltersHyperTree-$$VTKVERSION                     -lvtkInteractionStyle-$$VTKVERSION \
-lvtkFiltersHyperTreePython27D-$$VTKVERSION            -lvtkInteractionStylePython27D-$$VTKVERSION \
-lvtkFiltersImaging-$$VTKVERSION                       -lvtkInteractionWidgets-$$VTKVERSION \
-lvtkFiltersImagingPython27D-$$VTKVERSION              -lvtkInteractionWidgetsPython27D-$$VTKVERSION \
-lvtkFiltersModeling-$$VTKVERSION                      -lvtkNetCDF-$$VTKVERSION \
-lvtkFiltersModelingPython27D-$$VTKVERSION             -lvtkNetCDF_cxx-$$VTKVERSION \
-lvtkFiltersParallel-$$VTKVERSION                      -lvtkParallelCore-$$VTKVERSION \
-lvtkFiltersParallelImaging-$$VTKVERSION               -lvtkParallelCorePython27D-$$VTKVERSION \
-lvtkFiltersParallelImagingPython27D-$$VTKVERSION      -lvtkRenderingAnnotation-$$VTKVERSION \
-lvtkFiltersParallelPython27D-$$VTKVERSION             -lvtkRenderingAnnotationPython27D-$$VTKVERSION \
-lvtkFiltersPoints-$$VTKVERSION                        -lvtkRenderingContext2D-$$VTKVERSION \
-lvtkFiltersPointsPython27D-$$VTKVERSION               -lvtkRenderingContext2DPython27D-$$VTKVERSION \
-lvtkFiltersProgrammable-$$VTKVERSION                  -lvtkRenderingContextOpenGL2-$$VTKVERSION \
-lvtkFiltersProgrammablePython27D-$$VTKVERSION         -lvtkRenderingContextOpenGL2Python27D-$$VTKVERSION \
-lvtkFiltersPython-$$VTKVERSION                        -lvtkRenderingCore-$$VTKVERSION \
-lvtkFiltersPythonPython27D-$$VTKVERSION               -lvtkRenderingCorePython27D-$$VTKVERSION \
-lvtkFiltersSMP-$$VTKVERSION                           -lvtkRenderingFreeType-$$VTKVERSION \
-lvtkFiltersSMPPython27D-$$VTKVERSION                  -lvtkRenderingFreeTypePython27D-$$VTKVERSION \
-lvtkFiltersSelection-$$VTKVERSION                     -lvtkRenderingGL2PSOpenGL2-$$VTKVERSION \
-lvtkFiltersSelectionPython27D-$$VTKVERSION            -lvtkRenderingGL2PSOpenGL2Python27D-$$VTKVERSION \
-lvtkFiltersSources-$$VTKVERSION                       -lvtkRenderingImage-$$VTKVERSION \
-lvtkFiltersSourcesPython27D-$$VTKVERSION              -lvtkRenderingImagePython27D-$$VTKVERSION \
-lvtkFiltersStatistics-$$VTKVERSION                    -lvtkRenderingLOD-$$VTKVERSION \
-lvtkFiltersStatisticsPython27D-$$VTKVERSION           -lvtkRenderingLODPython27D-$$VTKVERSION \
-lvtkFiltersTexture-$$VTKVERSION                       -lvtkRenderingLabel-$$VTKVERSION \
-lvtkFiltersTexturePython27D-$$VTKVERSION              -lvtkRenderingLabelPython27D-$$VTKVERSION \
-lvtkFiltersVerdict-$$VTKVERSION                       -lvtkRenderingOpenGL2-$$VTKVERSION \
-lvtkFiltersVerdictPython27D-$$VTKVERSION              -lvtkRenderingOpenGL2Python27D-$$VTKVERSION \
-lvtkGeovisCore-$$VTKVERSION                           -lvtkRenderingVolume-$$VTKVERSION \
-lvtkGeovisCorePython27D-$$VTKVERSION                  -lvtkRenderingVolumeOpenGL2-$$VTKVERSION \
-lvtkIOAMR-$$VTKVERSION                                -lvtkRenderingVolumeOpenGL2Python27D-$$VTKVERSION \
-lvtkIOAMRPython27D-$$VTKVERSION                       -lvtkRenderingVolumePython27D-$$VTKVERSION \
-lvtkIOCore-$$VTKVERSION                               -lvtkViewsContext2D-$$VTKVERSION \
-lvtkIOCorePython27D-$$VTKVERSION                      -lvtkViewsContext2DPython27D-$$VTKVERSION \
-lvtkIOEnSight-$$VTKVERSION                            -lvtkViewsCore-$$VTKVERSION \
-lvtkIOEnSightPython27D-$$VTKVERSION                   -lvtkViewsCorePython27D-$$VTKVERSION \
-lvtkIOExodus-$$VTKVERSION                             -lvtkViewsInfovis-$$VTKVERSION \
-lvtkIOExodusPython27D-$$VTKVERSION                    -lvtkViewsInfovisPython27D-$$VTKVERSION \
-lvtkIOExport-$$VTKVERSION                             -lvtkWrappingPython27Core-$$VTKVERSION \
-lvtkIOExportPython27D-$$VTKVERSION                    -lvtkalglib-$$VTKVERSION \
-lvtkIOGeometry-$$VTKVERSION                           -lvtkexoIIc-$$VTKVERSION \
-lvtkIOGeometryPython27D-$$VTKVERSION                  -lvtkexpat-$$VTKVERSION \
-lvtkIOImage-$$VTKVERSION                              -lvtkfreetype-$$VTKVERSION \
-lvtkIOImagePython27D-$$VTKVERSION                     -lvtkgl2ps-$$VTKVERSION \
-lvtkIOImport-$$VTKVERSION                             -lvtkglew-$$VTKVERSION \
-lvtkIOImportPython27D-$$VTKVERSION                    -lvtkhdf5-$$VTKVERSION \
-lvtkIOInfovis-$$VTKVERSION                            -lvtkhdf5_hl-$$VTKVERSION \
-lvtkIOInfovisPython27D-$$VTKVERSION                   -lvtkjpeg-$$VTKVERSION \
-lvtkIOLSDyna-$$VTKVERSION                             -lvtkjsoncpp-$$VTKVERSION \
-lvtkIOLSDynaPython27D-$$VTKVERSION                    -lvtklibxml2-$$VTKVERSION \
-lvtkIOLegacy-$$VTKVERSION                             -lvtkmetaio-$$VTKVERSION \
-lvtkIOLegacyPython27D-$$VTKVERSION                    -lvtkoggtheora-$$VTKVERSION \
-lvtkIOMINC-$$VTKVERSION                               -lvtkpng-$$VTKVERSION \
-lvtkIOMINCPython27D-$$VTKVERSION                      -lvtkproj4-$$VTKVERSION \
-lvtkIOMovie-$$VTKVERSION                              -lvtksqlite-$$VTKVERSION \
-lvtkIOMoviePython27D-$$VTKVERSION                     -lvtksys-$$VTKVERSION \
-lvtkIONetCDF-$$VTKVERSION                             -lvtktiff-$$VTKVERSION \
-lvtkIONetCDFPython27D-$$VTKVERSION                    -lvtkverdict-$$VTKVERSION \
-lvtkIOPLY-$$VTKVERSION                                -lvtkzlib-$$VTKVERSION


# ------------------------------------------------------------------------------
# all headers and sources for this project
# ------------------------------------------------------------------------------

MSC_PATH = $$PWD/../msc
TOPOMS_PATH = $$PWD/../topoms

FORMS *= $$PWD/include/TopoMSUI.ui

INCLUDEPATH *= $$MSC_PATH/include \
               $$TOPOMS_PATH/include \
               $$PWD/include

HEADERS = $$MSC_PATH/include/*.h \
          $$TOPOMS_PATH/include/*.h \
          $$PWD/include/*.h

SOURCES = $$MSC_PATH/src/*.cpp \
          $$TOPOMS_PATH/src/*.cpp \
          $$PWD/src/*.cpp \
          $$PWD/src/*.cxx

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
