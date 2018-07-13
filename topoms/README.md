## TopoMS, Version 1.0
##### Author: Harsh Bhatia, hbhatia@llnl.gov

TopoMS is a computational tool for detailed topological analysis of molecular
and condensed matter systems, including the computation of atomic volumes and
charges through the quantum theory of atoms in molecules (also known as Bader
analysis), as well as the complete molecular graph.  With roots in techniques
from computational topology, and using a shared-memory parallel approach,
TopoMS provides scalable, numerically robust, and topologically consistent
analysis.

### Installation

Update: Please see the script `TopoMS/install.topoms.sh`, which automates the
following steps.

The only external dependency of `TopoMS` is VTK (https://www.vtk.org/); recommended
version `VTK 7.1`. Once VTK has been installed, `TopoMS` can be installed using
the `cmake` system. `TopoMS`  requires a C++ compiler supporting OpenMP.

```
$ pwd
TopoMS/topoms
$ mkdir build
$ cd build
$ cmake -DCMAKE_CXX_COMPILER=<path-to-gnu-c++> ..   # path needed for Mac
$ make
```

Note: If you're using the install scripts, the external dependencies will be installed
in `TopoMS/external`, and the executables will be installed in `TopoMS/build`.

### Execution

The application can be executed by passing a configuration file from command line.

```
$ TopoMS configuration_file.conf
```
