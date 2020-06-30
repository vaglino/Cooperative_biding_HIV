# Cooperative_biding_HIV

## Chemical kinetics modeling of the trimolecular Env-CD4-CXCR4 reaction

Performs fitting of adhesion frequency data and bond lifetime distributions with any kinetics model scheme defined, determines kinetic rates and compares different models statistically.

The code is written in Matlab and requires the optimization toolbox

### Association fitting
adhesion frequency modeling: the function adhesion_curve_fit.m allows to fit given adhesion frequency vs. time data to a ODE model provided as argument.

An example on how to use the adhesion_curve_fit function as well as associated code is provided in Pa_fit_HIV.m.
Examples on how to construct an ODE model for in the appropriate format are shown in the model folder.

### Dissociation fitting
bond lifetime distribution modeling: the function survival_curve_fit.m allows to fit given lifetime distribution to a ODE model provided as argument.

An example on how to use the survival_curve_fit function as well as associated code is provided in dissociation_fit_HIV.m.
Examples on how to construct an ODE model for in the appropriate format are shown in the model folder.

### Algorithm description
an initial set of parameters is determined via hypercube scanning of the parameter set. The parameters are then optimized to fit the data, starting with the best hypercube parameters. Both scanning and optimization are done in log space to better explore the wide parameter space. Then standard deviation of the fitted parameters is calculated by bootstrapping.

[![DOI](https://zenodo.org/badge/275973111.svg)](https://zenodo.org/badge/latestdoi/275973111)
