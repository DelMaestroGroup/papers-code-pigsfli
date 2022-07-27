[![Paper](https://img.shields.io/badge/paper-arXiv%3A2207.11301-B31B1B.svg)](https://arxiv.org/abs/2207.11301)
[![DOI](https://zenodo.org/badge/355554146.svg)](https://zenodo.org/badge/latestdoi/355554146)

#  ``PIGSFLI``: A Path-Integral Ground State Monte Carlo Algorithm for Entanglement of Lattice Bosons

Emanuel Casiano-Diaz, Chris M. Herdman, Adrian Del Maestro

[arXiv:2207.11301](https://arxiv.org/abs/2207.11301)

### Abstract
 A ground state path integral quantum Monte Carlo algorithm is introduced that allows for the simulation of lattice bosons at zero temperature. The method is successfully benchmarked against the one dimensional Bose-Hubbard model through comparison with the potential and kinetic energy computed via exact diagonalization. After successful validation, an estimator is introduced to measure the Rényi entanglement entropy between spatial subregions which is explored across the phase diagram of the one dimensional Bose-Hubbard model for systems consisting of up to 256 sites at unit-filling, far beyond the reach of exact diagonalization. The favorable scaling of the algorithm is demonstrated through a further measurement of the Rényi entanglement entropy at the two dimensional superfluid-insulator critical point for large system sizes, confirming the existence of the expected entanglement boundary law in the ground state. The Rényi entanglement estimator is extended to measure the symmetry resolved entanglement that is operationally accessible as a resource. 

### Description
This repository includes links, code, scripts, and data to generate the plots in the above paper.

### Requirements
The data in this project was generated via a new path integral Monte Carlo algorithm for the ground state of bosonic lattice models: [pigsfli](https://github.com/DelMaestroGroup/pigsfli).

The data contained in the [ProcessedData](https://github.com/DelMaestroGroup/papers-code-latticepigs/tree/main/ProcessedData) directory of this repository was obtained by processing raw data [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6827186.svg)](https://doi.org/10.5281/zenodo.6827186) with the various `.py` scripts found in the [src](https://github.com/DelMaestroGroup/papers-code-latticepigs/tree/main/src) directory. Figures were generated using the `.ipynb` notebook files contained there.

### Support
The creation of these materials was supported in part by the National Science Foundation under Award Nos. [DMR-1553991](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991&HistoricalAwards=false) and [DMR-2041995](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2041995&HistoricalAwards=false).

[<img width="100px" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)

### Figures

#### Figure 10: Relative Error of Kinetic & Potential Energies
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/relativeErrorsVK_N8.svg" width="400px">

#### Figure 14: Relative Error of Second Rényi Entanglement Entropy
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/relativeErrorsS2_N8.svg" width="400px">

#### Figure 15: Entanglement as a Function of Interaction Strength
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/interactionStrengthSweep_N16N256.svg" width="400px">

#### Figure 16: Entanglement Boundary Law in the Square Bose-Hubbard Lattice
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/boundaryLaw_N1024.svg" width="400px">

#### Figure 17: $\beta$-scaling of Accessible Second Rényi Entanglement Entropy
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/relativeErrorsS2acc_N8.svg" width="400px">

#### Figure 18: $\beta$-scaling of symmetry-resolved entanglement
<img src="https://github.com/DelMaestroGroup/papers-code-latticepigs/blob/main/figures/symmetry_resolved_s2_vs_beta.svg" width="400px">

Figures are relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

