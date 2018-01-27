## TopoMS, Version 1.0
##### Author: Harsh Bhatia, hbhatia@llnl.gov


TopoMS is a computational tool for detailed topological analysis of molecular
and condensed matter systems, including the computation of atomic volumes and
charges through the quantum theory of atoms in molecules (also known as Bader
analysis), as well as the complete molecular graph.  With roots in techniques
from computational topology, and using a shared-memory parallel approach,
TopoMS provides scalable, numerically robust, and topologically consistent
analysis.

TopoMS can be used as a command-line tool (`./topoms`) or with a GUI
(`./topoms-ui`), where the latter also enables an interactive exploration of
the molecular graph.  Installations instructions are provided separately in the
respective folders.

TopoMS requires a configuration file to provide several important parameters
for analysis. Both applications take a single command-line input -- the
configuration file name. A detailed example of configuration file is
provided in the `./config` folder.

### License

TopoMS is released under dual license.

- BSD
  - The core functionality, `./topoms`, is released under BSD license by LLNL.
  - The MSC library, `./msc`, is released under BSD license by University of Utah.
- GPL
  - The UI application, `./topoms-ui`, is released under GPL by LLNL.

Please see the respective directories on more details about corresponding licenses.
