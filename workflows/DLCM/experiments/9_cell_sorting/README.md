## README for ./9_cell_sorting

Simple cell sorting experiment where two different cell types have a
Young-Laplace pressure drop between them. Can achieve fragmentation
and engulfment by tweaking the sigma-parameter as in [1].  Runs using
the DLCM-solver in URDME.

E. Blom 2025-04-17

### Codes
* **cell_sorting.m** Young-Laplace cell sorting model.

visual/
* **visualise_cell_sorting.m** Visualises the experiments of Fig. 3.1 in [1]
* **cellsort1.mat** Symmetric experiment data (population
  fragmentation).
* **cellsort2.mat** Asymmetric experiment data (almost engulfment).
* **cellsort3.mat** Asymmetric experiment data (full engulfment!).

### References:
  [1] E. Blom, S. Engblom. "DLCM: a versatile multi-level solver for
  heterogeneous multicellular systems". ArXiv Preprint (2025)
