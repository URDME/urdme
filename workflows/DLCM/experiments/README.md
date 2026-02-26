## README for ./DLCM/experiments

E. Blom 2025-04-16

### Sample models:

* 1_basic_test/
	Basic test of relaxation of a population of cells.

* 2_avascular_tumour/
	Avascular tumour model, results in §3.3 in [1].
  * old_examples/
  Earlier models from [3,4].

* 3_gradient_growth/
	Chemotaxis models (2D/3D), results in §3.4 in [1].

* 4_delta_notch/
	Delta-notch signalling in a growing population of cells.
	Results in §3.2 in [1].

* 6_NDR/
	Advanced Delta-notch signalling modeling in a growing population
	of cells. This example is used as the running model in [6]. Note:
	this is an earlier example not using the dlcm-solver
	directly. Refer to the README-file within this directory.

* 7_hes1/
  Hes1-Notch pathway models from [1,2].

* 9_cell_sorting/
  Cell sorting experiments. Results in §3.1 in [1].

### References:
  [1] E. Blom, S. Engblom. "DLCM: a versatile multi-level solver for
  heterogeneous multicellular systems". ArXiv Preprint (2025) \
  [2] Menz, G., Engblom, S. Modelling Population-Level Hes1 Dynamics:
  Insights from a Multi-framework Approach. Bull Math Biol. 87 (2025) \
  [3] E. Blom, S. Engblom, and G. Menz: "Modeling the hallmarks of
  avascular tumors", ENUMATH23 proceedings (2025) \
  [4] E. Blom, S. Engblom: "Morphological stability for in silico models
  of avascular tumors", Bull. Math. Biol. 86 (2024) \
  [5] S. Engblom: "Stochastic simulation of pattern formation in
  growing tissue: a multilevel approach", Bull. Math. Biol. 81 (2019) \
  [6]  S. Engblom, D. B. Wilson, and R. E. Baker: "Scalable
  population-level modeling of biological cells incorporating
  mechanics and kinetics in continuous time", Roy. Soc. Open
  Sci. 5(8) (2018). \
