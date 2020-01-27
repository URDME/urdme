README for workflows/DLCM/experiments/6_NDR/

S. Engblom 2018-06-21

Advanced Delta-notch signalling modeling in a growing population of
cells. This example is used as the running model in [1]. The model
itself is described in [2]

There are a few (very mild) dependencies on functions from STENGLIB,
download freely from www.stenglib.org. PDE Toolbox is used to assemble
the Laplacian operator over a Delaunay triangulation.

Contents:

6_NDR/
	static_run.m		Calibration on a static population.
	growth.m		Growth process (dynamic_growth.mat).
	dynamic_run.m	        Process on a growing population.

	NDR_ODEvsSSA_static.m	self-explanatory calibration scripts
	NDR_ODEvsSSA_dynamic.m
	NDR_RDME_static.m
	NDR_RDME_dynamic.m

data/
	dynamic_growth.mat	Saved decoupled growth process.
	dynamicRDME1.mat	Main results, from [1].
	dynamicRDME2.mat
	dynamicRDME6.mat

source/
	NDR_ODEvsSSA.c		Right-hand side rates.
	NDR_RDME.c
	NDRrhs.m

utils/
	contact.m		Cell-contact logic.
	DLCMlayer.m		Discretization layers.
	RDMElayer.m
	NDRparam.m		Model parameters, mostly from [2].
	NDRplot.m		Plot utility.

References:
  [1] S. Engblom: "Stochastic simulation of pattern formation in
  growing tissue: a multilevel approach", Bull. Math. Biol. 81 (2019).
  [2] Z. Hadjivasiliou, G. L. Hunter, and B. Baum: "A new mechanism
  for spatial pattern formation via lateral and protrusion-mediated
  lateral signalling", J. R. Soc. Interface 13(124):1--10 (2016).
