### Welcome to the URDME Webpage.
URDME is a general software framework for modeling and simulation of stochastic reaction-diffusion processes on unstructured, tetrahedral and triangular meshes. URDME emphasizes modularity in order to be useful both as a simulation tool and as a framework for developemnt of stochastic simulation algorithms.

### Quick Start 

Download the latest stable release and unpack it, alternatively clone the repository (requires a git client). In a terminal type:

```
$ git clone https://github.com/URDME/urdme.git
$ sudo ./install.sh
```
Consult the software manual in _doc/manual.pdf_ for further instructions, examples and tutorials. 

### System requirements and software dependencies
* Linux or Apple OSX operating system.
* Matlab
	- Tested versions: 2007a, 2008a, 2008b, 2009a, 2009b, 2011b, 2012a, 2014b. Command line interface must be installed. 
     
* Comsol multiphysics 3.5a (with appropriate patches) or Comsol multiphysics 4.x 
	- Tested versions: 3.5a, 4.0, 4.2, 4.3, 4.3b (recommended)
	- Must have Matlab integration components installed
* GCC (Xcode on Apple computers+command-line tools)
    - Executables gcc and make must be in the path
    - Standard libraries must be installed
* The optional SBML support requires additionally
	- Python runtime libraries 2.6 or higher
  	- SBML library for python (libsbml)

### Copyright
URDME 1.2. Copyright 2008,2009,2010,2011,2012.
See the file AUTHORS for complete list of authors. 

### Publications
URDME: A modular framework for stochastic simulation of reaction-transport processes in complex geometries, 
B. Drawert, S. Engblom and A. Hellander, BMC Systems Biology 2012 6:76