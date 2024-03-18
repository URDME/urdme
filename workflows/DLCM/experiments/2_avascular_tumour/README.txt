README for workflows/DLCM/experiments/2_avascular_tumor

Avascular tumor simulations.

E. Blom 24-01-26

Models
    1) AvascularTumor.m - DLCM avascular tumor model from [3].

    2) AvascularTumor_updated.m - DLCM avascular tumor model
    with surface tension and pressure sinks from [2].
    Plots volumetric growth and solution at 3 predefined times.
    Run final Section to save data and animate full simulation.
     - Run as is for initially circular tumor with no surface tension
        growing in a nutrient-rich region. Radial symmetry breaks by
        nutrient starvation and part of the tumor keeps growing.
     - Run with sigma = 32*1e-4 for stable circular growth with
        creeping effect (migration towards oxygen source)
     - Run with sigma = 5*1e-4 for tumor that begins to split in two
     - Load 'morphology_DLCM1.mat' or 'morphology_DLCM2.mat' and run the
        section 'Plots' to get figures from the first two experiments in 
        §4.2 in [2]
     
    3) AvascularTumorPDE.m - PDE avascular tumor model as
    effective model of the one in AvascularTumor_updated.m [2]
    Plots volumetric growth and solution at 3 predefined times.
    Run final Section to save data and animate full simulation.
     - Run as is for simulation corresponding to AvascularTumor_updated.m
        in the first example above. Also gives similar results for the 
        suggested sigma values.
     - Load 'morphology_PDE1.mat' or 'morphology_PDE2.mat' and run the
        section 'Plots' to get figures from the first two experiments in 
        §4.2 in [2]

    4) RegionalDynamics2D_run.m - Runs RegionalDynamics2D.m example
    of PDE model in [2] under assumption of radial symmetry and plots
    volumetric growth and perturbation growth factors at stationary state.

    5) AvascularTumor_hallmarks.m - DLCM avascular tumor model from [1]
    with explicit hallmark capabilites. Run as is to simulate the
    experiment using the same parameters as in [1], or load 
    AThallmarks_longrun.mat and run section 'Plots' in 
    AvascularTumor_hallmarks.m to get the figures from the experiment 
    performed in [1].

Analysis
utils/
    stability_analysis.m
	Symbolic evaluation of analysis in §3.2 & Lemma A.3 [2].
	
    rad_stability.m
	Symbolic evaluation and visualisation of analysis in §3.1 [2].

References:
  [1] E. Blom, S. Engblom, and G. Menz, Modeling the hallmarks of avascular
      tumors, 2024. ENUMATH23 proceedings (submitted).
  [2] E. Blom and S. Engblom, Morphological stability for in silico models
      of avascular tumors, 2023. arXiv Preprint:
      https://doi.org/10.48550/arXiv.2309.07889
  [3]  S. Engblom, D. B. Wilson, and R. E. Baker: "Scalable
       population-level modeling of biological cells incorporating
       mechanics and kinetics in continuous time", Roy. Soc. Open
       Sci. 5(8) (2018).
