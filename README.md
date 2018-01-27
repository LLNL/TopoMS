## TopoMS, Version 1.0

<img align="right" src="https://user-images.githubusercontent.com/10440378/35475248-ebbbba0c-034f-11e8-8980-1199a5c5fb56.png" width="50%">

##### Author: Harsh Bhatia, hbhatia@llnl.gov; Attila G Gyulassy

TopoMS is a computational tool for detailed topological analysis of molecular
and condensed matter systems, including the computation of atomic volumes and
charges through the quantum theory of atoms in molecules (also known as Bader
analysis), as well as the complete molecular graph.  With roots in techniques
from computational topology, and using a shared-memory parallel approach,
TopoMS provides scalable, numerically robust, and topologically consistent
analysis.

If you use TopoMS, please cite the following publication.
* H Bhatia, A G Gyulassy, V Lordi, J E Pask, V Pascucci, and P-T Bremer, "TopoMS: Comprehensive Topological Exploration for Molecular and Condensed-Matter Systems," Journal of Computational Chemistry, to appear. doi:10.1002/jcc.25181.

TopoMS can be used as a command-line tool (`./topoms`) or with a GUI
(`./topoms-ui`), where the latter also enables an interactive exploration of
the molecular graph.  Installations instructions are provided separately in the
respective folders.

TopoMS requires a configuration file to provide several important parameters
for analysis. Both applications take a single command-line input -- the
configuration file name. A detailed example of configuration file is
provided in the `./configs` folder.

### Examples

* Bader volumes of the atoms in Ethylene molecule

<img src="https://user-images.githubusercontent.com/10440378/35475245-d51a282e-034f-11e8-833b-77b05169c70c.png" width="25%">

* Complete molecular graph of a lithium salt (LiPF6) in Ethylene Carbonate solution

<img src="https://user-images.githubusercontent.com/10440378/35475242-cd6257b4-034f-11e8-9383-a40415ec5242.png" width="25%">

### License

TopoMS is released under dual license.

- BSD
  - The core functionality, `./topoms`, is released under BSD license by LLNL.
  - The MSC library, `./msc`, is released under BSD license by University of Utah.
- GPL
  - The UI application, `./topoms-ui`, is released under GPL by LLNL.

Please see the respective directories on more details about corresponding licenses.
