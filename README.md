# URDME

Unstructured Reaction-Diffusion Master Equation

www.urdme.org

URDME is a general software framework for modeling and simulation of
stochastic reaction-diffusion processes on arbitrary meshes. URDME
emphasizes modularity in order to be useful both as a simulation tool
and as a framework for development of stochastic simulation
algorithms.

URDME 1.4. Copyright 2008--2020

See the file [AUTHORS](https://github.com/URDME/urdme/blob/master/AUTHORS) for a complete list of authors. 

### Quick Start

Download the latest stable release and unpack it, alternatively clone
the repository (requires a git client),
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

See the file [VERSION](https://github.com/URDME/urdme/blob/master/VERSION) for a complete list.

### Licence: GPL 3

See the file [LICENCE](https://github.com/URDME/urdme/blob/master/LICENCE) for the full statement.

### Citing URDME:

URDME is a research software and correctly citing it is important to
the maintainers.

* The software package; B. Drawert, S. Engblom and A. Hellander:
"URDME: A modular framework for stochastic simulation of
reaction-transport processes in complex geometries", BMC Systems
Biology 6(76):1--17 (2012)
[(doi)](http://dx.doi.org/10.1186/1752-0509-6-76)

* The underlying numerical modeling; S. Engblom, L. Ferm, A. Hellander
and P. Lötstedt: "Simulation of stochastic reaction--diffusion
processes on unstructured meshes", SIAM J. Scientific. Comp.
31(3):1774--1797 (2009) [(doi)](http://dx.doi.org/10.1137/080721388)

* The subdiffusion workflow; S. Engblom, P. Lötstedt and L. Meinecke:
"Mesoscopic Modeling of Random Walk and Reactions in Crowded Media,
Phys. Rev. E 98(3):033304 (2018)
[(doi)](http://dx.doi.org/10.1103/PhysRevE.98.033304)

* The neuron workflow; P. Bauer, S. Engblom, S. Mikulovic and
A. Senek: "Multiscale modelling via split-step methods in neural
firing", Math. Comput. Model. Dyn. Syst. 24(4):409--425 (2018)
[(doi)](http://dx.doi.org/10.1080/13873954.2018.1488740)

* The DLCM workflow; S. Engblom, D. B. Wilson and R. E. Baker:
"Scalable population-level modelling of biological cells incorporating
mechanics and kinetics in continuous time", Roy. Soc. Open Sci. 5(8)
(2018) [(doi)](http://dx.doi.org/10.1098/rsos.180379), and related,
S.Engblom: "Stochastic simulation of pattern formation in growing
tissue: a multilevel approach", Bull. Math. Biol. 81:3010--3023 (2019)
[(doi)](http://dx.doi.org/10.1007/s11538-018-0454-y)

* The AEM solver; P. Bauer and S. Engblom: "Sensitivity estimation and
inverse problems in spatial stochastic models of chemical kinetics",
pp. 519--527 in A. Abdulle, S. Deparis, D. Kressner, F. Nobile and
M. Picasso (editors): "Numerical Mathematics and Advanced
Applications: ENUMATH 2013", vol 103 of Lecture Notes in Computational
Science and Engineering, Springer (2015) 
[(doi)](http://dx.doi.org/10.1007/978-3-319-10705-9_51)

* The UDS solver implements the first order method described in
S. Engblom: "Computing the Moments of High Dimensional Solutions of
the Master Equation" Appl. Math. Comput. 18(2):498--515 (2006) 
[(doi)](http://dx.doi.org/10.1016/j.amc.2005.12.032)
