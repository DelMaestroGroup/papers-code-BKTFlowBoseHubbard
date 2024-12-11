[![Paper](https://img.shields.io/badge/paper-arXiv%3AXXXX.YYYYY-B31B1B.svg)](https://arxiv.org/abs/XXXX.YYYYY)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXYYYYY)

# Direct Observation of the Berezinskii-Kosterlitz-Thouless Renormalization Group Flow at a Quantum Phase Transition

Author 1, Author 2, Author 3

[arXiv:XXXX.YYYYY](https://arxiv.org/abs/XXXX.YYYYY)

### Abstract
Abstract Here

### Description
This repository includes links, code, scripts, and data to generate the figures in a paper.

### Requirements
The data in this project can be generated using the code in the following repositories:

1. DMRG simulations: [ExtendedBH_DMRG_Fluctuations_Julia](https://github.com/DelMaestroGroup/ExtendedBH_DMRG_Fluctuations_Julia)
2. QMC simulations: [pigsfli](https://github.com/DelMaestroGroup/pigsfli)

Data is included in the [data](https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/tree/main/data) directory.

This code requires Julia version 10.4 or higher and the IJulia package to run the Jupyter notebook. Required Julia packages can be installed by running the code in the `create_figures.jl` script:

```julia
using Pkg 
Pkg.activate(".")
Pkg.add(["Plots","PyFormattedStrings","NonlinearSolve","StaticArrays","Printf","Integrals","FastClosures","LaTeXStrings","DataFrames","NPZ","Measures","PyCall"])
```

### Support
The creation of these materials was supported in part by the {INSERT FUNDING AGENCY} under Award No. [{AWARD NUMBER}](https://www.nsf.gov/awardsearch/simpleSearchResult?queryText=delmaestro).

<img width="400px" src="https://new.nsf.gov/themes/custom/nsf_theme/components/images/logo/logo-desktop.svg">
<img width="400px" src="https://science.osti.gov/assets/img/doe-logos/logo.png">


### Figures

#### Figure 01: BKT RG flow for Bose-Hubbard model.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/001_BKT_RG_flow_Bose-Hubbard.png" width="600px">

#### Figure 02: Finite size scaling of Luttinger parameter according to BKT flow.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/002_finite_size_scaling.svg" width="400px">

#### Figure 03: Extracting Luttinger parameter from particle number fluctuations.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/003_fit_K_to_fluctuations.svg" width="400px">

#### Figure 04: Extracting critical point of BKT transition.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/004_zeta_of_U.svg" width="400px">

#### Figure S01: Periodic vs. open boundary conditions.
<img src="https://github.com/DelMaestroGroup/papers-code-BKTFlowBoseHubbard/blob/main/figures/S001_compare_pbc_obc.svg" width="400px">

These figures are relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

