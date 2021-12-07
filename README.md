[![Paper](https://img.shields.io/badge/paper-arXiv%3AXXXX.YYYYY-B31B1B.svg)](https://arxiv.org/abs/XXXX.YYYYY)
[![DOI](https://zenodo.org/badge/214220909.svg)](https://zenodo.org/badge/latestdoi/214220909)

# Entanglement of Lattice Bosons via Path-Integral Ground State Monte Carlo

Emanuel Casiano-Diaz, Chris M. Herdman, Adrian Del Maestro

[arXiv:XXXX.YYYYY](https://arxiv.org/abs/XXXX.YYYYY)

### Abstract
A ground state quantum Monte Carlo algorithm based on Path Integral Monte Carlo (PIMC) is introduced that allows for the simulation of lattice bosons at zero tempera- ture. The algorithm is validated by computing energy benchmarks in the Bose-Hubbard model that agree with the exact results expected from exact diagonalization. After a successful validation of the algorithm, an estimator is introduced to measure the Rényi entanglement entropy between spatial subregions. The resulting entanglement entropy is measured across the phase diagram of the one-dimensional Bose-Hubbard model for systems consisting of up to L = 256 sites at unit filling, far beyond the reach of exact di- agonalization, and at the insulating-superfluid critical point in two spatial dimensions, confirming the existence of an entanglement area (boundary) law in the ground state. By imposing a particle number superselection rule, the Rényi estimator is extended to measure the symmetry resolved entanglement that is operationally accessible as a re- source.

### Description
This repository includes links, code, scripts, and data to generate the plots in the above paper.

### Requirements
The raw data in this project was generated via a new Path Integral Monte Carlo simulation for the ground state of bosonic lattice models. The source code for this Monte Carlo simulation can be found in the following [public repository](https://github.com/ecasiano/pimc/tree/master/pimc), alongside the Python scripts used to obtain the data in the [ProcessedData](https://github.com/DelMaestroGroup/papers-code-latticepigs/tree/main/ProcessedData) directory.

### Support
The creation of these materials was supported in part by the National Science Foundation under Award No. DMR-1553991.

[<img width="100px" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)

### Figures

#### Figure 08: Relative Error of Kinetic & Potential Energies
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/relativeErrorsVK_N8.svg" width="400px">

#### Figure 12: Relative Error of Second Rényi Entanglement Entropy
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/relativeErrorsS2_N8.svg" width="400px">

#### Figure 13: Entanglement as a Function of Interaction Strength
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/interactionStrengthSweep_N16N256.svg" width="400px">

#### Figure 14: Entanglement Boundary Law in the Square Bose-Hubbard Lattice
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/boundaryLaw_N1024.svg" width="400px">

#### Figure 15: Relative Error of Accessible Second Rényi Entanglement Entropy
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/relativeErrorsS2acc_N8.svg" width="400px">

This figure is relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

