### Welcome to the URDME Webpage.
URDME is a general software framework for modeling and simulation of stochastic reaction-diffusion processes on unstructured, tetrahedral and triangular meshes. URDME emphasizes modularity in order to be useful both as a simulation tool and as a framework for development of stochastic simulation algorithms.

### Quick Start 

Download the latest stable release and unpack it, alternatively clone the repository (requires a git client). In a terminal type:

```
$ git clone https://github.com/URDME/urdme.git
```

To start using URDME under Matlab, navigate to the file `startup.m`,
```
>> startup
```

To run a simple 1D model,
```
>> cd examples/annihilation
>> annihilation_run
```

For further instructions, examples and tutorials, consult the software
manual,
```
doc/manual.pdf 
```

### System requirements and software dependencies
See the file VERSION for a fuller list.

* OS, tested with Ubuntu 14.04 LTS, 16.04 LTS, Mac OS X 10.9.5

* Matlab, tested with versions: 8.4, 9.1

* Comsol Multiphysics 5.2a with Matlab integration components

* PDE Toolbox, tested with versions: 1.5, 2.3

### Copyright
URDME 1.3. Copyright 2008,2009,2010,2011,2012.
See the file AUTHORS for complete list of authors. 

### Publications
URDME is a research software and correctly citing it is important to
the maintainers.

* The software package; B. Drawert, S. Engblom and A. Hellander:
"URDME: A modular framework for stochastic simulation of
reaction-transport processes in complex geometries", BMC Systems
Biology 6(76):1--17 (2012)
[(doi)](http://dx.doi.org/10.1186/1752-0509-6-76)

* The AEM solver; P. Bauer and S. Engblom: "Sensitivity estimation and
inverse problems in spatial stochastic models of chemical kinetics",
pp. 519--527 in A. Abdulle, S. Deparis, D. Kressner, F. Nobile and
M. Picasso (editors): "Numerical Mathematics and Advanced
Applications: ENUMATH 2013", vol 103 of Lecture Notes in Computational
Science and Engineering, Springer (2015).
[(doi)](http://dx.doi.org/10.1007/978-3-319-10705-9_51)

* The UDS solver implements the first order method described in
S. Engblom: "Computing the Moments of High Dimensional Solutions of
the Master Equation" Appl. Math. Comput. 18(2):498--515 (2006).
[(doi)](http://dx.doi.org/10.1016/j.amc.2005.12.032)
