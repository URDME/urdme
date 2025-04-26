## README for ./4_delta_notch

Simple Delta-Notch dynamics using the Collier model [2], using both \
discrete and continuous internal states. Runs using the DLCM-solver [1].

E. Blom 2025-04-17

### Codes
* **delta_notch.m** Delta-Notch simulation. \
  Set flag 'internal_state' to 'cont' to use continuous internal delta
  & notch values.

visual/
* **visualise_delta_notch.m** Visualise the experiment data shown in
  Fig. 3.2 in [1].
* **contDNgrowth.mat** Continuous states data.
* **discDNgrowth.mat** Discrete states data.
* **cmap.mat** 'Paraview' colormap.

### References:
  [1] E. Blom, S. Engblom. "DLCM: a versatile multi-level solver for
  heterogeneous multicellular systems". ArXiv Preprint (2025) \
  [2] J. Collier, et al. "Pattern formation by lateral inhibition
  with feedback: a mathematical model of delta-notch intercellular
  signalling." Journal of theoretical Biology. (1996)
