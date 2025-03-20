# Imperfect phase randomization and generalized decoy-state quantum key distribution

This is a public version of the code used in [*Imperfect phase randomization and generalized decoy-state quantum key distribution*](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.20.064031) \[[arXiv](https://arxiv.org/abs/2304.09401)\]. This repository pre-dates the open QKD security software.

This code computes key rates for an imperfectly phase randomized three state protocol as described in the above paper. It is intended to reproduce Fig. 6. Note that slight differences might be present due to parts of the code being updated to make it compatible with MATLAB 2024b.

> [!CAUTION]
> This code uses parallel computing to reduce run time, however CVX is known to be unstable when run in parallel.

## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes.

This code was designed for a previous version of Matlab, but was updated and tested on 2024b, though any other recent edition should work. In addition to Matlab, this code requires the following additional resources:
 - [CVX](https://cvxr.com/cvx/download/) v2.2, a library for solving convex optimization problems in MATLAB.
 - [QETLAB](https://github.com/nathanieljohnston/QETLAB) *above* v0.9, a MATLAB toolbox for operations involving quantum channels and entanglement. Note that you cannot use the version from their website as it has a bugs associated with implementing Choi matrices. *You must use their latest copy on Github*. At the time of writing, this means downloading their code with the green "code" button and *not* the v0.9 release.
 - [MOSEK](https://www.mosek.com/) a more advanced semidefinite programming (SDP) solver than the default (SDPT3) used in CVX. Note that the MOSEK solver can be downloaded together with CVX, but requires a separate license to use. See [this page](https://cvxr.com/cvx/doc/mosek.html) for more information.
 - The following Mathworks's toolboxes that can be installed with Matlab:
   - Parallel Computing
   - Statistics and Machine Learning

1. Download the latest release on the side bar or clone with git
```
git clone https://github.com/Optical-Quantum-Communication-Theory/Imperfect-phase-randomization-and-generalized-decoy-state-QKD
```
2. Unzip in your preferred directory and add the folder and all subfolders to your Matlab path.
3. To generate the data, run `SampleScript_ThreeState_decoy.m` then plot with `simplePlot.m`.
