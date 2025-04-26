## README for ./3_gradient_growth

Five models of chemotaxis, one using custom Laplacian and one in 3D,
and visualisations. Runs using the DLCM-solver in URDME.

E. Blom 2025-04-17

### Codes
* **chemotaxis.m** Three 2D chemotaxis models.
* **chemotaxis_modLa.m** Chemotaxis with custom Laplacian.
* **chemotaxis3D.m** 3D chemotaxis model.

visual/
* **visualise_chemotaxis3D.m** Visualising 3D experiments, Fig. 3.5 in [1].
* **visualise_chemotaxis.m** Visualising 2D experiments, Fig. 3.4 in [1].
* **chemotaxis_3D.mat** Model 5 data; cells exert chemotactic signal.
* **chemotaxis_cons.mat** Model 3 data; cells consume ambient.
  chemotactic signal
* **chemotaxis_diffusion.mat** Model 2 data; diffusive cell movement +
  chemotaxis.
* **chemotaxis_modLa.mat** Model 4 data; chemotaxis & pressure physics
  combined.
* **chemotaxis_pressure.mat** Model 1 data; Pressure-driven movement +
  chemotaxis.

### References:
  [1] E. Blom, S. Engblom. "DLCM: a versatile multi-level solver for
  heterogeneous multicellular systems". ArXiv Preprint (2025)
