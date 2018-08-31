
#!/bin/bash

# ------------------------------------------------------------------------------
# check commands

command -v wget >/dev/null 2>&1 || { echo >&2 "Cannot find command 'wget'. Aborting."; exit 1;  }
command -v tar >/dev/null 2>&1 || { echo >&2 "Cannot find command 'tar'. Aborting."; exit 1;  }
command -v cmake >/dev/null 2>&1 || { echo >&2 "Cannot find command 'cmake'. Aborting."; exit 1;  }

if [ -n "${GNU_CXX}" ] ; then
  command -v "${GNU_CXX}" >/dev/null 2>&1 || { echo >&2 "Cannot find command '${GNU_CXX}'. Aborting."; exit 1;  }
  echo ''
  echo "Using" `"${GNU_CXX}"  --version | head -n1`:  `type gcc`
else
  echo ''
  echo "Using" `gcc --version | head -n1`:  `type gcc`
fi

command -v qmake >/dev/null 2>&1 || { echo >&2 "Cannot find command 'qmake'. Aborting."; exit 1;  }
echo `qmake  --version | tail -n1`

platform=`uname`
if [ $platform == 'Darwin' ]; then
  alias qmake='qmake -spec macx-g++ INCLUDEPATH+=/usr/include'
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
  TopoMS_ROOT=${0%/*}
  TopoMS_ROOT=$(cd "$(dirname "$TopoMS_ROOT")"; pwd)/$(basename "$TopoMS_ROOT")
else
  TopoMS_ROOT=`pwd`
fi

# ------------------------------------------------------------------------------
# install script for TopoMS and its depenencies
# ------------------------------------------------------------------------------

echo '> Installing TopoMS-UI and its dependencies for ('`whoami`') on ('`hostname`'). platform = ('$platform')'
echo ''
echo '  pwd                   : '`pwd`
echo '  TopoMS Root Directory : '$TopoMS_ROOT
echo ''
cd $TopoMS_ROOT
mkdir -p external/downloads
cd external/downloads

# ------------------------------------------------------------------------------
# VTK 7.1
# ------------------------------------------------------------------------------
VTK_VERSION='7.1'
VTK_VERSION_BUILD='1'
VTK_NAME='VTK-'$VTK_VERSION'.'$VTK_VERSION_BUILD

testfile=$TopoMS_ROOT/external/lib/libvtkCommonCore-$VTK_VERSION*
if ! ls $testfile 1> /dev/null 2>&1 ; then
  echo '  > '$VTK_NAME

  if ! ls $VTK_NAME.tar.gz 1> /dev/null 2>&1 ; then
    echo '    > Downloading '$VTK_NAME
    rm $TopoMS_ROOT/external/vtk.download.log 2>/dev/null
    wget -o $TopoMS_ROOT/external/vtk.download.log https://www.vtk.org/files/release/$VTK_VERSION/$VTK_NAME.tar.gz
  fi

  if ! ls $VTK_NAME 1> /dev/null 2>&1 ; then
    echo '    > Untarring the downloaded code'
    tar -xf $VTK_NAME.tar.gz
  fi

  echo '    > Setting out of source build ('$VTK_NAME'/build)'
  mkdir -p $VTK_NAME/build
  cd $VTK_NAME/build

  echo '    > Configuring VTK'
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$TopoMS_ROOT/external \
        .. > $TopoMS_ROOT/external/vtk.cmake.log

  echo '    > Building VTK'
  rm $TopoMS_ROOT/external/vtk.make.log 2>/dev/null
  make -j6 > $TopoMS_ROOT/external/vtk.make.log

  echo '    > Installing VTK'
  make install > $TopoMS_ROOT/external/vtk.make.log

  if ls $testfile 1> /dev/null 2>&1 ; then
    echo '    > VTK successfully installed in '$TopoMS_ROOT'/external.'
  else
    echo '    > VTK build failed. Please check the build logs in '$TopoMS_ROOT'/external for more information.'
  fi
else
  echo '  > '$VTK_NAME 'is already installed ('$TopoMS_ROOT'/external). If you need to reinstall, please remove the existing installation.'
fi

# ------------------------------------------------------------------------------
# QGLViewer 2.7.1
# ------------------------------------------------------------------------------
QGL_VERSION='2.7'
QGL_VERSION_BUILD='1'
QGL_NAME='libQGLViewer-'$QGL_VERSION'.'$QGL_VERSION_BUILD

testfile=$TopoMS_ROOT/external/lib/libQGLViewer*
if ! ls $testfile 1> /dev/null 2>&1 ; then
  echo '  > '$QGL_NAME

  if ! ls $QGL_NAME.tar.gz 1> /dev/null 2>&1 ; then
    echo '    > Downloading '$QGL_NAME
    rm $TopoMS_ROOT/external/qgl.download.log 2>/dev/null
    wget -o $TopoMS_ROOT/external/qgl.download.log http://www.libqglviewer.com/src/$QGL_NAME.tar.gz
  fi

  if ! ls $QGL_NAME 1> /dev/null 2>&1 ; then
    echo '    > Untarring the downloaded code'
    tar -xf $QGL_NAME.tar.gz
  fi

  echo '    > Setting out of source build ('$QGL_NAME'/build)'
  mkdir -p $QGL_NAME/build
  cd $QGL_NAME/build

  echo '    > Configuring QGLViewer'
  rm $TopoMS_ROOT/external/qgl.cmake.log 2>/dev/null

  if [ -n "${GNU_CXX}" ] ; then
    echo '    > CXX_COMPILER: '$GNU_CXX
    qmake QMAKE_CXX=$GNU_CXX \
          QMAKE_CXXFLAGS+='-w' \
          QGLVIEWER_STATIC=yes PREFIX=$TopoMS_ROOT/external \
          ../QGLViewer  > $TopoMS_ROOT/external/qgl.qmake.log
  else
    qmake QMAKE_CXXFLAGS+='-w' \
          QGLVIEWER_STATIC=yes PREFIX=$TopoMS_ROOT/external \
          ../QGLViewer  > $TopoMS_ROOT/external/qgl.qmake.log
  fi

  echo '    > Building QGLViewer'
  rm $TopoMS_ROOT/external/qgl.make.log 2>/dev/null
  make -j6 > $TopoMS_ROOT/external/qgl.make.log

  echo '    > Installing QGLViewer'
  make install > $TopoMS_ROOT/external/qgl.make.log

  if ls $testfile 1> /dev/null 2>&1 ; then
    echo '    > QGLViewer successfully installed in '$TopoMS_ROOT'/external.'
  else
    echo '    > QGLViewer build failed. Please check the build logs in '$TopoMS_ROOT'/external for more information.'
  fi
else
  echo '  > '$QGL_NAME 'is already installed ('$TopoMS_ROOT'/external). If you need to reinstall, please remove the existing installation.'
fi

# ------------------------------------------------------------------------------
# TopoMS (with gui)
# ------------------------------------------------------------------------------
cd $TopoMS_ROOT

echo '  > TopoMS-UI'

mkdir -p build/topoms-ui
cd build/topoms-ui

echo '    > Configuring TopoMS-UI'

if [ -n "${GNU_CXX}" ] ; then
  echo '    > CXX_COMPILER: '$GNU_CXX
  cmake -DCMAKE_CXX_COMPILER=$GNU_CXX \
        ../../topoms-ui > topoms-ui.cmake.log
else
  cmake ../../topoms-ui > topoms-ui.cmake.log
fi

echo '    > Building TopoMS'
make -j6 > topoms-ui.make.log

if ls TopoMS-UI 1> /dev/null 2>&1 ; then
  echo '    > TopoMS-UI successfully built in '`pwd`
else
  echo '    > TopoMS-UI build failed. Please check the build logs in '`pwd`' for more information.'
fi

# ------------------------------------------------------------------------------
# end of the install script
# ------------------------------------------------------------------------------