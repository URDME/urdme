## README for ./DLCM

Discrete Laplacian Cell Mechanics: population-level modeling of
biological cells in continuous time.

E. Blom 2025-04-16 (Major revision, integrated URDME solver) \
S. Engblom 2018-06-21 (Minor revision) \
S. Engblom and D. Wilson 2017-08-30

The DLCM workflow is implemented as an URDME solver which can be used
to model living cells in a population. The modeling physics is
detailed in [1,5] and consists of a grid-based stochastic method.

In essence, the framework represents cells as discrete agents that can
proliferate, move, switch phenotype, etc., with internal states
representing protein synthesis, degradation and other signaling
processes. The cells are embedded in a micro-enviroment of continuous
quantities, such as pressure and nutrients, each modeled by stationary
heat equations. Surface tension effects between populations can be
included through Young-Laplace drops in the pressure field.

PDE Toolbox is used to assemble the Laplacian operator over a Delaunay
triangulation required to solve for the micro-environment quantities.
Note that the 3D models use fegeometry (R2023b) to construct the mesh.

### Animations

animations/ \

* **animate_dlcm.m** Copy-paste examples of animating dlcm-simulations.
* This folder also contains animations from each of the experiments in [1].

### Sample models:

experiments/

A collection of various examples highlighting modeling use. Folders
often contain a sub-folder 'visual' which contains code and data to
reproduce the figures from the experiments in [1]. See the associated
[README](experiments/README.md).

### Utility functions:

utils/
* **basic_mesh.m** Basic regular mesh (Cartesian/hexagonal).
* **cmap_dlcm.mat** 'DLCM'-colormap: grey to bluish green to vermillion
  (based on colors from graphics_color.m)
* **dlcm2urdme.m** Organises solver object into umod.solverargs.
* **dt_neighe** Neighbor edge matrix with weights.
* **dt_operators.m** Operators over Delaunay triangulation.
* **graphics_color.m** Color selection in plots.
* **map.m** Mapping of indices.
* **mesh_layers** Mesh boundary layers.
* **mesh2dual.m** Dual mesh from primal mesh.

### References:
  [1] E. Blom, S. Engblom. "DLCM: a versatile multi-level solver for
  heterogeneous multicellular systems". ArXiv Preprint (2025) \
  [2] E. Blom, S. Engblom, and G. Menz: "Modeling the hallmarks of
  avascular tumors", ENUMATH23 proceedings (2025) \
  [3] E. Blom and S. Engblom: "Morphological stability for in silico models
  of avascular tumors", Bull. Math. Biol. 86 (2024) \
  [4] S. Engblom: "Stochastic simulation of pattern formation in
  growing tissue: a multilevel approach", Bull. Math. Biol. 81 (2019) \
  [5]  S. Engblom, D. B. Wilson, and R. E. Baker: "Scalable
  population-level modeling of biological cells incorporating
  mechanics and kinetics in continuous time", Roy. Soc. Open
  Sci. 5(8) (2018). \
