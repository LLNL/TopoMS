## TopoMS-UI, Version 1.0
##### Author: Harsh Bhatia, hbhatia@llnl.gov

TopoMS is a computational tool for detailed topological analysis of molecular
and condensed matter systems, including the computation of atomic volumes and
charges through the quantum theory of atoms in molecules (also known as Bader
analysis), as well as the complete molecular graph.  With roots in techniques
from computational topology, and using a shared-memory parallel approach,
TopoMS provides scalable, numerically robust, and topologically consistent
analysis.

### Installation

In addition to `TopoMS`, `TopoMS-UI` depends upon `Qt 5.7`, `QGLViewer 2.7.1`
(http://libqglviewer.com/), and `GLEW 2.1.0`.
Once these depenencies are installed, you may use `cmake` or Qt's `qmake`
system to build the tool.

```
$ pwd
TopoMS/topoms-ui
$ mkdir build
$ cmake -DCMAKE_CXX_COMPILER=<path-to-gnu-c++> .. # compiler path needed only for Mac
$ make
```
If you are more comfortable using `qmake`, you may proceed as follows

```
$ pwd
TopoMS/topoms-ui
$ mkdir build
$ qmake -spec macx-g++ QMAKE_CXX=<path-to-gnu-c++> \  # this line needed only for Mac
        QGLPATH=/Users/bhatia4/usr \
        VTKPATH=/Users/bhatia4/usr VTKVERSION=7.1 \
        ..
$ make
```

Note: If you're using the install scripts, the external dependencies will be installed
in `TopoMS/external`, and the executables will be installed in `TopoMS/build`.

### Execution

The application can be executed by passing a configuration file from command line.

```
$ TopoMS-UI configuration_file.conf
```
