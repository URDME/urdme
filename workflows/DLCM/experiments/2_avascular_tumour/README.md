## README for workflows/DLCM/experiments/2_avascular_tumor
Avascular tumor models, both PDE and DLCM models.

E. Blom 23-09-22

### Sample models
* **AvascularTumour.m** - DLCM avascular tumor model in [2]
* **AvascularTumor_updated.m** - DLCM avascular tumor model in [1]
* **AvascularTumorPDE.m** - PDE avascular tumor model
* **RegionalDynamics2D_run.m** - runs `RegionalDynamics2D` example

### Analysis
utils/
* **stability_analysis.m** - symbolic morphological analysis as in [1]
* **rad_stability.m** - symbolic radial symmetry analysis as in [1]

### Utility functions
utils/
* **animate_growth.m** - animates tumor growth simulation from datafile
* **curvature2D.m** - evaluates curvature of contourline
* **connect_boundary.m** - connects FEM mesh boundary
* **errorshade.m** - plots moving std from data as shade
* **get_growth_factors.m** - evaluates Lambda(k) from tumor geometry [1]
* **meshdata.m** - maps cell data structure to one amenable to MATLAB'S
    `meshgrid.m` function
* **solve_oxygen.m** - solves Poisson problem for oxygen [1]
* **rad_stability.m** - see Analysis
* **RegionalDynamics2D.m** - solves 2D PDE solution under assumption of
    radial symmetry
* **set_init.m** - sets DLCM initial condition from a set of options
* **stability_analysis.m** - see Analysis

### Visualisation
visual/
* **RegionalCharacteristicsVisualiser.m** - plots tumor regional sizes
* **VisualiseSigma.m** - plots growth factors vs time

References:
  [1] E. Blom and S. Engblom, Morphological stability for in silico models
      of avascular tumors, 2023. arXiv Preprint:
      https://doi.org/10.48550/arXiv.2309.07889
  [2]  S. Engblom, D. B. Wilson, and R. E. Baker: "Scalable
  population-level modeling of biological cells incorporating
  mechanics and kinetics in continuous time", Roy. Soc. Open
  Sci. 5(8) (2018).
