## README for ./7_hes1
Hes1-Notch pathway models from [1].

G. Menz and S. Engblom 2024-11-14

### Data
* **alpha_parametrization.m** - Result from running `parametrization.m` with `final = true`
  - used to determine final parameterisation in Tab. 2.1
* **hes1_bifurcation_Nvoxels_20_VOL_20.mat** - Bifurcation results for RDME model from running `hes1umod_bifurcation.m` on a 20x20 grid with volume 20 micrometer^3.
  - used in the generation of Figure(s): 3.7
* **hes1_grid2D.mat** - Result from running `hes1_grid2D.m` (full ODE model) on a 20x20 grid with random ics.
  - used in the generation of Figure(s): 2.2, 2.3
* **hes1red_grid2D.mat** - Result from running `hes1red_grid2D.m` (reduced ODE model, alternative 1) on a 20x20 grid with random ics.
  - used in the generation of Figure(s): 2.3
* **hes1umod_vol1_rand.mat** - Result from running `hes1umod2D_run.m` (RDME model) on a 20x20 grid using a cell volume of 1 micrometer^3 with random ics.
  - used in the generation of Figure(s): 2.5 
* **hes1umod_vol50_rand.mat** - Result from running `hes1umod2D_run.m` (RDME model) on a 20x20 grid using a cell volume of 50 micrometer^3 with random ics.
  - used in the generation of Figure(s): 2.5
* **stationary_3cell.mat** - Bifurcation results of 3 cell problem determined in `stationary_3cell.m`.
  - used in the generation of Figure(s): 3.5, 3.7
* **urdme_patterning.mat** - Connectivity in RDME model results over 14 different volumes and 10 iterations each.
  - used in the generation of Figure(s): 3.8

### Experiments
* **hes1_grid2D.m** - Full ODE Model implemented on a 2D grid using DLCM.
* **hes1red_grid2D.m** - Reduced ODE model implemented on a 2D grid using DLCM.
* **hes1umod2D_run.m** - Run RDME model in 2D using URDME solvers (both NSM and UDS).
* **hes1umod_bifurcation.** - Generate bifurcation graph for RDME model.
  - Figure(s) produced in [1]: 3.7
* **hes1umod_patterning.m** - Determine how "perfect" pattern is for different volumes of NSM compared to UDS results by calculating grid coupling W.
  - Figure(s) produced: 3.8
* **panel_plots.m** - Generates panel plots and constituent comparison plot in manuscript.
  - Figure(s) produced: 2.2, 2.3 and 2.5
* **stationary_3cell.m** - Investigate the stability of the reduced Hes1-model over three cells.
  - Figure(s) produced: 3.5 and 3.6
* **stationary_bifurcation.m** - Find bifurcation diagram by scaling parameters using the stationary solution(s) of the Hes1-model on a 1D grid with periodic boundary conditions and as a function of the parameters.
  - Figure(s) produced: 3.2
* **stationary_bounds.m** - Stationary analysis of Hes1 using the reduced ODE model.
  - Figure(s) produced: 3.1
* **stationary_stability.m** - Spectral analysis using spectra for small perturbations from stationary states.

### Models
* **final_parametrization** - Check parameterisations of perturbed results from `hes1_conc.m` and mus from `hes1_params.m` followed by parameterisation of alphas.
* **hes1_Jacobian.m** - Jacobian for the full Hes1 ODE model.
* **hes1_System.m** - RHS function for the Hes1 ODE system on a grid.
* **hes1_buildODE.m** - Full Hes1 ODE model on a grid.
* **hes1_conc.m** - Wanted concentrations for all constituents. Possibility to generate perturbed concentrations.
* **hes1_params.m** - Function which returns all Hes1 model parameters as a struct. Possibility to generate
perturbed parameters.
* **hes1red_params.m** - Function returning scaled Hes1 model parameters as a struct. Possibility to generate
perturbed parameters with same noise as used in `hes1_params.m`.
* **hes1umod.m** - URDME model file for Hes1@5. Run by `hes1umod2D_run.m`
* **parametrization.m** - Hes1 parametrization for given data.

### Utils
* **hes1red_homogeneous.m** - Homogeneous solution of the reduced ODE model.
* **hes1red_nonhomogeneous.m** - Non-homogeneous solution of the
  reduced ODE model.
* **hes1reduce.m** - Reduced Hes1-model given full parameters.
* **hes1unreduce.m** - Map to full Hes1-model from reduced case.

References: \
  [1] G. Menz and S. Engblom, Modelling Population-Level Hes1 Dynamics: Insights from a
  Multi-Framework Approach, 2024. arXiv Preprint: *link to be added*
