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

In addition to `TopoMS`, `TopoMS-UI` depends upon Qt and QGLViewer
(http://libqglviewer.com/). Once Qt and QGLViewer have been installed,
please edit the `TopoMSUI.pro` file to provide the paths to QGLViewer and VTK
(lines 29--31). TopoMS-UI can be installed using the `qmake` system.

```
$ pwd
TopoMS/topoms-ui
$ qmake
$ make
```

### Execution

The application can be executed by passing a configuration file from command line.

```
$ TopoMS-UI configuration_file.conf
```