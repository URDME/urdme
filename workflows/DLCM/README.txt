README for workflows/DLCM/

Discrete Laplacian Cell Mechanics: population-level modeling of
biological cells in continuous time.

S. Engblom 2018-06-21 (Minor revision)
S. Engblom and D. Wilson 2017-08-30

The DLCM workflow is a package distributed with URDME which can be
used to model living cells in a population. The modeling physics is
detailed in [1] and consists of a grid-based stochastic method. The
workflow may be used as a stand-alone package, or together with
solvers and other features of URDME.

This is the first release intended for URDME 1.4. A tighter
integration with URDME is scheduled for URDME 1.5. There are a few
(very mild) dependencies on functions from STENGLIB, download freely
from www.stenglib.org.

PDE Toolbox is used to assemble the Laplacian operator over a Delaunay
triangulation.

Sample models:

experiments/

1_basic_test/
	Basic test of simulation of a population of cells (2D/3D).
	Results summarized in §3.1 of [1].

2_avascular_tumour/
	Avascular tumour model, results in §3.4 of [1].

3_gradient_growth/
	Chemotaxis slit-model, results in §3.2 of [1].

4_delta_notch/
	Delta-notch signalling in a growing population of cells.
	Results summarized in §3.3 of [1].

6_NDR/
	Advanced Delta-notch signalling modeling in a growing
	population of cells. This example is used as the running model
	in [2]. Refer to the README-file within this directory.

Utility functions:

utils/
	basic_mesh.m		Basic regular mesh (Cartesian/hex).
	dt_operators.m	Operators over Delaunay triangulation.
	graphics_color.m	Color selection in plots.
	map.m			Mapping of indices.
	mesh2dual.m		Dual mesh from primal mesh.

References:
  [1] S. Engblom, D. B. Wilson, and R. E. Baker: "Scalable
  population-level modeling of biological cells incorporating
  mechanics and kinetics in continuous time", Roy. Soc. Open
  Sci. 5(8) (2018).
  [2] S. Engblom: "Stochastic simulation of pattern formation in
  growing tissue: a multilevel approach", Bull. Math. Biol. 81 (2019).
