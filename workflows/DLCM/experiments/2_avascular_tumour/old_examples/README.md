## README for ./old_examples
Avascular tumor models, both PDE and DLCM models, without the DLCM-solver.

### Sample models
* **AvascularTumour.m** - DLCM avascular tumor model in [3]
* **AvascularTumor_updated.m** - DLCM avascular tumor model in [2]
* **AvascularTumorPDE.m** - PDE avascular tumor model in [2]
* **AvascularTumor_hallmarks.m** - DLCM avascular tumor model with explicit 
    hallmark capabilities from [1]. To run experiment in [1], run default 
    parameter setting, or use the corresponding Simulation Data as 
    explained below.
* **RegionalDynamics2D_run.m** - runs `RegionalDynamics2D` example

### Simulation Data
* **AThallmarks_longrun.mat** - Contains data from experiment in [1].
    Load and run section Plots in `AvascularTumor_hallmarks.m` for the
    figures in [1].
* **morphology_X#.mat** - Contains data from experiments in [2] for the
    X = DLCM or X = PDE model respectively, and # = 1 (rad. symmetry) 
    or # = 2 (splits in two). Load and run section Plots in corresponding 
    model code for figures.
* **growth_rate_comparison_exp.mat** - Contains morphological data from the
    growth rate experiments in [2]. 20 samples per mode for modes 
    k = 1, ..., 8. Analyzed in `get_growth_factors.m`.

### Analysis
utils/
* **stability_analysis.m** - symbolic morphological analysis as in [2]
* **rad_stability.m** - symbolic radial symmetry analysis as in [2]

### Utility functions
utils/
* **animate_growth.m** - animates tumor growth simulation from datafile
* **curvature2D.m** - evaluates curvature of contourline
* **compare_growth_rates.m** - 1st section evaluates analytical and 
    experimental growth rates from data in `morph_k#1_exp#2.mat`.
    2nd section plots the corresponding figure from [2].
* **connect_boundary.m** - connects FEM mesh boundary
* **errorshade.m** - plots moving std from data as shade
* **get_growth_factors.m** - evaluates Lambda(k) from tumor geometry [2]
* **meshdata.m** - maps cell data structure to one amenable to MATLAB'S
    `meshgrid.m` function
* **solve_oxygen.m** - solves Poisson problem for oxygen [2]
* **rad_stability.m** - see Analysis
* **RegionalDynamics2D.m** - solves 2D PDE solution under assumption of
    radial symmetry
* **set_init.m** - sets DLCM initial condition from a set of options
* **stability_analysis.m** - see Analysis

### Visualisation
visual/
* **RegionalCharacteristicsVisualiser.m** - plots tumor regional sizes
* **VisualiseSigma.m** - plots growth factors vs time

References: \
  [1] E. Blom, S. Engblom, and G. Menz: "Modeling the hallmarks of
  avascular tumors", ENUMATH23 proceedings (2025) \
  [2] E. Blom and S. Engblom: "Morphological stability for in silico models
  of avascular tumors", Bull. Math. Biol. 86 (2024) \
  [3]  S. Engblom, D. B. Wilson, and R. E. Baker: "Scalable
  population-level modeling of biological cells incorporating
  mechanics and kinetics in continuous time", Roy. Soc. Open
  Sci. 5(8) (2018).
