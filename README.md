## TopoMS, Version 1.1
#### Released: Nov 18, 2018

<img align="right" src="./docs/jcc25344-toc-0001-m.jpg" width="50%">
<!--<img align="right" src="https://user-images.githubusercontent.com/10440378/35475248-ebbbba0c-034f-11e8-8980-1199a5c5fb56.png" width="50%">-->

##### Author: Harsh Bhatia, hbhatia@llnl.gov; Attila G Gyulassy

TopoMS is a computational tool for detailed topological analysis of molecular
and condensed matter systems, including the computation of atomic volumes and
charges through the quantum theory of atoms in molecules (also known as Bader
analysis), as well as the complete molecular graph.  With roots in techniques
from computational topology, and using a shared-memory parallel approach,
TopoMS provides scalable, numerically robust, and topologically consistent
analysis.

If you use TopoMS, please cite the following publication.
* H Bhatia, A G Gyulassy, V Lordi, J E Pask, V Pascucci, and P-T Bremer,
"TopoMS: Comprehensive Topological Exploration for Molecular and Condensed-Matter Systems,"
Journal of Computational Chemistry, vol. 39, issue 16, pp 936–952, June 2018.
[doi:10.1002/jcc.25181](https://doi.org/10.1002/jcc.25181).

TopoMS can be used as a command-line tool (`./topoms`) or with a GUI
(`./topoms-ui`), where the latter also enables an interactive exploration of
the molecular graph.

TopoMS requires a configuration file to provide several important parameters
for analysis. Both applications take a single command-line input -- the
configuration file name. A detailed example of configuration file is
provided in the `./configs` folder.

### Installation

* TopoMS requires a `gnu` C++ compiler (preferably `gcc 7.3`). Note that for Mac,
the default `gcc` is in fact `clang`. Unless you have it configured otherwise,
you need to set the path of `gnu` compiler explicitly (see below).

* In addition, `TopoMS-UI` requires the following to support visualization and
user interface. TopoMS assumes the following are available on your machine.
  * `Qt 5.7`: Note that Qt often changes the API substantially between versions;
    please use the recommended version. Qt is available through a number of
    package managers. For example, for Mac, you can use `macports` to get Qt:
    `port install qt57`. Or you may install Qt from source. See [here](http://doc.qt.io/qt-5/linux.html#downloading-and-installing-qt) and
    [here](http://doc.qt.io/qt-5/linux-building.html).
  * `OpenGL`: Comes pre-installed on your machine.
  * `GLEW 2.1.0`: Extension Wrangler for OpenGL. Also available through common
  package managers. For Mac, `GLEW` is available via
  `macports`: `port install glew`. You may also install it from source. See [here](http://glew.sourceforge.net/build.html).

* Visualization Toolkit (VTK): `TopoMS`  requires VTK in order to output
  molecular graphs. Although optional, `VTK 7.1.1` is strongly suggested as a
  dependency.

* Finally, `TopoMS-UI` also uses `QGLViewer 2.7` for visualization. See [here](http://libqglviewer.com/download.html).

To simplify installation, scripts have been provided to build `TopoMS` as well
as `TopoMS-UI` along with their dependencies (`VTK` and `QGLViewer`). Please see
the usage below. Additional instructions for individual installation are
provided separately in the respective folders, `topoms` and `topoms-ui`.

```
$ pwd
<your-path>/TopoMS
$ export GNU_CXX=/opt/local/bin/g++  # example path;    needed only for Mac
$ sh install.topoms.sh
$ sh install.topoms-ui.sh
```

`TopoMS` executables will be installed in `<your-path>/TopoMS/build`, and
dependencies (`VTK` and `QGLViewer`) will be installed in
`<your-path>/TopoMS/external`.

**Note:** The provided scripts have been tested on Mac OS 10.13 and Red Hat linux.
Support for installation on Windows will be provided in the future.

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
