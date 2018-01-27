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

The only external dependency of `TopoMS` is VTK (https://www.vtk.org/). Once VTK has been installed, `TopoMS` can be installed using the `cmake` system. `TopoMS` also requires a C++ compiler supporting OpenMP.

```
$ pwd
TopoMS/topoms
$ mkdir build
# cd build
$ cmake ../
$ make
```

### Execution

The application can be executed by passing a configuration file from command line.

```
$ TopoMS configuration_file.conf
```