README for workflows/subdiffusion/

S. Engblom 2018-02-16 (Minor revision)
S. Engblom 2017-05-16

The SUBDIFFUSION workflow enriches URDME with the modeling of
subdiffusive effects via an internal states model [1,2]. Starting from
an ordinary diffusive reaction process, utility functions are used to
expand the model with internal states such that the overall transport
behavior agrees with subdiffusiv scaling.

There are a few (very mild) dependencies on functions from STENGLIB,
download freely from www.stenglib.org.

Sample model files:

pure1D_run.m
	simple test of subdiffusive movement in 1D
pure2D_run.m
	similar test in 2D, model from [1]
bimolecular2D_run.m
	bimolecular reaction in 2D, model from [1]
mincde3D_run.m
	the Min-problem interpreted as a subdiffusive process, model from [1]

Utility functions:

subdmatrix.m
	expands a diffusion matrix with internal states, given a
	transfer matrix between the states and scalar diffusion
	coefficients for the states

subdata.m
	constructs an internal states transfer matrix from
	coarse-grained pre-computed data

Data:

data/
  CoarseGrainedSD.mat	coarse-grained pre-computed data, from [1]
  bimolecular2D.mat	results from bimolecular2D_run 
  originalrun.mat	results from mincde3D_run
  subdrun.mat		(with and without subdiffusive behavior)

References:
  [1] S. Engblom, P.Lötstedt, and L. Meinecke: "Mesoscopic Modeling of
  Random Walk and Reactions in Crowded Media" Phys. Rev. E
  98(3):033304 (2018).
  [2] E. Blanc, S. Engblom, A. Hellander, and P. Lötstedt: "Mesoscopic
  modeling of stochastic reaction-diffusion kinetics in the
  subdiffusive regime", Multiscale Model. Simul. 14(2):668--707
  (2016).
