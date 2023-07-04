# LFMS_Stiefel
Leapfrog Multiple Shooting on the Stiefel Manifold.

This is a collection of MATLAB files associated with the preprint
[The leapfrog algorithm as nonlinear Gauss-Seidel](https://arxiv.org/abs/2010.14137v3)
and [my Ph.D. thesis](https://archive-ouverte.unige.ch/unige:146438).
It collects the main script to run single shooting, multiple shooting, and the leapfrog algorithms to compute geodesics on the Stiefel manifold.

If there are any problems, bugs, or suggestions, feel free to [contact me](mailto:msutti@ncts.tw).


## I) Version History

- July 2023: changed formulation in `GetKA.m`.
- 8 Feb 2023: updated README file.
- Ver 1, 12 Nov 2021: initial release.


## II) Contents

The main folder `LFMS` contains the following subfolders:

- `Plots`: the folder where the plots are saved.
- `Postprocessing`: contains utilities for plotting and visualizing data.
- `Results`: contains the matfiles generated by `Driver_Figure_6_Boxplot.m`.
- `Utilities`: contains all the core functions for leapfrog, single shooting, and multiple shooting.

The additional folder named `Shared_Utilities` contains other core functions shared by other thesis parts.


## III) Installation and Usage

No installation is required. Use the `Driver_*.m` files to run the calculations or generate the figures in the preprint/Ph.D. thesis.
The code saves the results of the simulations as matfiles in the `Results` folder. 
The figures are saved in the `Plots` folder.


## IV) License

The code in this repository is GPL licensed.


## V) References
- [1] M. Sutti, "Riemannian Algorithms on the Stiefel and the Fixed-Rank Manifold", Ph.D. thesis, University of Geneva, December 2020.
- [2] M. Sutti and B. Vandereycken, "The leapfrog algorithm as nonlinear Gauss-Seidel", arXiv preprint, January 2023.
