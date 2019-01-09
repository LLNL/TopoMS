#!/bin/bash

# ------------------------------------------------------------------------------
echo ''
echo '> This script will install TopoMS and its dependencies.'
echo '  Do you also want to install TopoMS-UI? (y/n) [default: n]'
read BUILD_UI

# ------------------------------------------------------------------------------
# check commands

platform=`uname`

command -v wget >/dev/null 2>&1 || { echo >&2 "Cannot find command 'wget'. Aborting."; exit 1;  }
command -v tar >/dev/null 2>&1 || { echo >&2 "Cannot find command 'tar'. Aborting."; exit 1;  }
command -v cmake >/dev/null 2>&1 || { echo >&2 "Cannot find command 'cmake'. Aborting."; exit 1;  }

if [ -n "${GNU_CXX}" ] ; then
    command -v "${GNU_CXX}" >/dev/null 2>&1 || { echo >&2 "Cannot find command '${GNU_CXX}'. Aborting."; exit 1;  }
    echo ''
    echo "> Using" `"${GNU_CXX}"  --version | head -n1`:  `type gcc`
else
    echo ''
    echo "> Using" `gcc --version | head -n1`:  `type gcc`
fi

if [[ $BUILD_UI == y ]]; then
    command -v qmake >/dev/null 2>&1 || { echo >&2 "Cannot find command 'qmake'. Aborting."; exit 1;  }
    echo "> "`qmake  --version | tail -n1`

    if [ $platform == 'Darwin' ]; then
        alias qmake='qmake -spec macx-g++ INCLUDEPATH+=/usr/include'
    fi
fi

# ------------------------------------------------------------------------------
if false; then
echo
echo "# path to me --------------->  ${0}     "
echo "# parent path -------------->  ${0%/*}  "
echo "# my name ------------------>  ${0##*/} "
fi

# get the correct path of the TopoMS root
if [ ${0} != ${0##*/} ];  then
    PATH_TopoMS=${0%/*}
    PATH_TopoMS=$(cd "$(dirname "$PATH_TopoMS")"; pwd)/$(basename "$PATH_TopoMS")
else
    PATH_TopoMS=`pwd`
fi

NPROCS=20

# ------------------------------------------------------------------------------
# install script for TopoMS and its depenencies
# ------------------------------------------------------------------------------
echo ''
echo '> Installing TopoMS and its dependencies for ('`whoami`') on ('`hostname`'). platform = ('$platform')'
echo ''
echo '  > pwd                  : '`pwd`
echo '  > TopoMS Root Directory: '$PATH_TopoMS
echo ''

PATH_Ext=$PATH_TopoMS/external
mkdir -p $PATH_Ext/downloads
cd $PATH_Ext/downloads

# ------------------------------------------------------------------------------
# VTK 7.1
# ------------------------------------------------------------------------------
VTK_VERSION='7.1'
VTK_VERSION_BUILD='1'
VTK_NAME='VTK-'$VTK_VERSION'.'$VTK_VERSION_BUILD

testfile=$PATH_Ext/lib/libvtkCommonCore-$VTK_VERSION*
if ! ls $testfile 1> /dev/null 2>&1 ; then
    echo '  > '$VTK_NAME

    if ! ls $VTK_NAME.tar.gz 1> /dev/null 2>&1 ; then
        echo '     > Downloading '$VTK_NAME
        rm $PATH_Ext/vtk.download.log 2>/dev/null
        wget -o $PATH_Ext/vtk.download.log https://www.vtk.org/files/release/$VTK_VERSION/$VTK_NAME.tar.gz
    fi

    if ! ls $VTK_NAME 1> /dev/null 2>&1 ; then
        echo '     > Untarring the downloaded code'
        tar -xf $VTK_NAME.tar.gz
    fi

    echo '     > Setting out of source build ('$VTK_NAME'/build)'
    mkdir -p $VTK_NAME/build
    cd $VTK_NAME/build

    echo '     > Configuring VTK'
    rm $PATH_Ext/vtk.cmake.log 2>/dev/null
    cmake -DCMAKE_INSTALL_PREFIX=$PATH_Ext \
          -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_SHARED_LIBS:BOOL=OFF \
          -DVTK_WRAP_PYTHON:BOOL=ON \
          -DVTK_USE_CXX11_FEATURES:BOOL=ON \
          .. > $PATH_Ext/vtk.cmake.log

    echo '     > Building VTK'
    rm $PATH_Ext/vtk.make.log 2>/dev/null
    make -j$NPROCS > $PATH_Ext/vtk.make.log

    echo '     > Installing VTK'
    make install > $PATH_Ext/vtk.make.log

    if ls $testfile 1> /dev/null 2>&1 ; then
        echo '     > VTK successfully installed in '$PATH_Ext'.'
    else
        echo '     > VTK build failed. Please check the build logs in '$PATH_Ext'for more information.'
    fi
else
    echo '  > '$VTK_NAME 'is already installed ('$PATH_Ext'). If you need to reinstall, please remove the existing installation.'
fi

# ------------------------------------------------------------------------------
# QGLViewer 2.7.1
# ------------------------------------------------------------------------------
if [[ $BUILD_UI == y ]]; then

    QGL_VERSION='2.7'
    QGL_VERSION_BUILD='1'
    QGL_NAME='libQGLViewer-'$QGL_VERSION'.'$QGL_VERSION_BUILD

    testfile=$PATH_Ext/lib/libQGLViewer*
    if ! ls $testfile 1> /dev/null 2>&1 ; then
        echo '  > '$QGL_NAME

        if ! ls $QGL_NAME.tar.gz 1> /dev/null 2>&1 ; then
            echo '    > Downloading '$QGL_NAME
            rm $PATH_Ext/qgl.download.log 2>/dev/null
            wget -o $PATH_Ext/qgl.download.log http://www.libqglviewer.com/src/$QGL_NAME.tar.gz
        fi

        if ! ls $QGL_NAME 1> /dev/null 2>&1 ; then
            echo '    > Untarring the downloaded code'
            tar -xf $QGL_NAME.tar.gz
        fi

        echo '    > Setting out of source build ('$QGL_NAME'/build)'
        mkdir -p $QGL_NAME/build
        cd $QGL_NAME/build

        echo '    > Configuring QGLViewer'
        rm $PATH_Ext/qgl.cmake.log 2>/dev/null

        COMPILER=""
        if [ -n "${GNU_CXX}" ] ; then
            echo '     > CXX_COMPILER: '$GNU_CXX
            COMPILER="QMAKE_CXX=$GNU_CXX"
        fi

        qmake $COMPILER QMAKE_CXXFLAGS+='-w' \
              QGLVIEWER_STATIC=yes PREFIX=$PATH_Ext \
              ../QGLViewer > $PATH_Ext/qgl.qmake.log

        echo '    > Building QGLViewer'
        rm $PATH_Ext/qgl.make.log 2>/dev/null
        make -j$NPROCS > $PATH_Ext/qgl.make.log

        echo '    > Installing QGLViewer'
        make install > $PATH_Ext/qgl.make.log

        if ls $testfile 1> /dev/null 2>&1 ; then
            echo '    > QGLViewer successfully installed in '$PATH_Ext'.'
        else
            echo '    > QGLViewer build failed. Please check the build logs in '$PATH_Ext' for more information.'
        fi
    else
        echo '  > '$QGL_NAME 'is already installed ('$PATH_Ext'). If you need to reinstall, please remove the existing installation.'
    fi
fi

cd $PATH_TopoMS

# ------------------------------------------------------------------------------
# TopoMS
# ------------------------------------------------------------------------------
echo '  > TopoMS'

mkdir -p build
cd build

echo '     > Configuring TopoMS'

# parameters for cmake
COMPILER=""
if [ -n "${GNU_CXX}" ] ; then
    echo '     > CXX_COMPILER: '$GNU_CXX
    COMPILER="-DCMAKE_CXX_COMPILER=$GNU_CXX"
fi

UI=""
if [[ $BUILD_UI == y ]]; then
    UI="-DTOPOMS_BUILD_UI=True"
fi

cmake $COMPILER $UI \
      -DCMAKE_INSTALL_PREFIX=$PATH_TopoMS/install \
      $PATH_TopoMS > topoms.cmake.log

echo '     > Building TopoMS'

make -j$NPROCS > topoms.make.log
make install

if ls TopoMS 1> /dev/null 2>&1 ; then
    echo '    > TopoMS successfully built in '`pwd`
else
    echo '    > TopoMS build failed. Please check the build logs in '`pwd`' for more information.'
fi

# ------------------------------------------------------------------------------
# end of the install script
# ------------------------------------------------------------------------------
