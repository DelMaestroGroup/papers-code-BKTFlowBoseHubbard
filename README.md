[![Paper](https://img.shields.io/badge/paper-arXiv%3AXXXX.YYYYY-B31B1B.svg)](https://arxiv.org/abs/XXXX.YYYYY)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXYYYYY)

# Berezinskii-Kosterlitz-Thouless Renormalization Group Flow at a Quantum Phase Transition

Matthias Thamm, Harini Radhakrishnan, Hatem Barghathi, Chris Herdman, Arpan Biswas, Bernd Rosenow, and Adrian Del Maestro 

[arXiv:XXXX.YYYYY](https://arxiv.org/abs/XXXX.YYYYY)

### Abstract
We present a controlled numerical study of the Berezinskii-Kosterlitz-Thouless (BKT) transition in the one-dimensional Bose-Hubbard model at unit filling, providing evidence of the characteristic logarithmic finite-size scaling of the BKT transition. Employing density matrix renormalization group and quantum Monte Carlo simulations under periodic boundary conditions, together with a systematic finite-size scaling analysis of bipartite particle number fluctuations, we resolve boundary-induced complications that previously obscured critical scaling. We demonstrate that a suitably chosen central region under open boundaries reproduces universal RG signatures, reconciling earlier discrepancies. Finally, leveraging a non-parametric Bayesian analysis, we determine the critical interaction strength with high precision, establishing a benchmark for BKT physics in one-dimensional quantum models.

### Description
This repository includes links, code, scripts, and data to generate the figures in a paper.

### Requirements
The data in this project can be generated using the code in the following repositories:

1. DMRG simulations: [ExtendedBH_DMRG_Fluctuations_Julia](https://github.com/DelMaestroGroup/ExtendedBH_DMRG_Fluctuations_Julia)
2. QMC simulations: [pigsfli](https://github.com/DelMaestroGroup/pigsfli)

Data is included in the [data](https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/tree/main/data) directory.
QMC raw data is available via a Zenodo archive:  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14879281.svg)](https://zenodo.org/record/14879281)

This code requires Julia version 10.4 or higher and the IJulia package to run the Jupyter notebook. Required Julia packages can be installed by running the code in the `create_figures.jl` script:

```julia
using Pkg 
Pkg.activate(".")
Pkg.add(["Plots","PyFormattedStrings","NonlinearSolve","StaticArrays","Printf","Integrals","FastClosures","LaTeXStrings","DataFrames","NPZ","Measures","PyCall"])
```

The python code for postprocessing of the QMC data in `data/pbc/QMC/postprocess` requires the following packages:
- matplotlib
- numpy
- tqdm
- scipy
- zipfile-deflat64

The python code for the BO and GP analysis of the $\zeta(K)$ data in `src/GP_with_BoTorch.ipynb` requires the following packages:
- botorch=0.10.0
- torch
- gpytorch
- numpy
- matplotlib
- dgutils (`pip install git+https://github.com/DelMaestroGroup/dgutils.git#egg=dgutils`)
- scipy


### Support
This work was partially supported by the National Science Foundation Materials Research Science and Engineering Center program through the UT Knoxville Center for Advanced Materials and Manufacturing (DMR-2309083). H.R. acknowledges AITennessee for financial support. Computations were performed using resources provided by the Leipzig University Computing Center and University of Tennessee Infrastructure for Scientific Applications and Advanced Computing (ISAAC). 
 


### Figures

#### Figure 01: BKT RG flow for Bose-Hubbard model.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/001_BKT_RG_flow_Bose-Hubbard.png" width="600px">

#### Figure 02: Finite size scaling of Luttinger parameter according to BKT flow.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/002_finite_size_scaling.svg" width="400px">

#### Figure 03: Extracting Luttinger parameter from bipartite particle number fluctuations.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/003_fit_K_to_fluctuations.svg" width="400px">

#### Figure 04: Extracting critical point of BKT transition.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/final_zeta_var.svg" width="400px">

#### Figure S01: Periodic vs. open boundary conditions.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/S001_compare_pbc_obc.svg" width="400px">

#### Figure S02: Fitting method for OBC.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/S003_obc_fit_interval_data.svg" width="400px">

#### Figure S03: Reproduction of literature result.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/S002_comparison_obc_literature.svg" width="400px">

#### Figure S04: GP mean combined with the confidence interval of the posterior distribution..
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/GP0_acquisition_wzoom.svg" width="400px">


These figures are relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

