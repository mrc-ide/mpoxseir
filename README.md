# mpoxseir

This package contains a compartmental SEIRD model of mpox disease transmission and progression as illustrated in the following schematic diagram:

<img width="351" height="344" alt="image" src="https://github.com/user-attachments/assets/f79b3f36-a1c2-41b5-854d-823d250e6f3e" />

(A) The population moves through the different disease states susceptible (S), exposed (Ea,Eb), infected (IR, ID) and recovered (R) or dead (D). The disease states follow an identical flow regardless of vaccination stratum j (although the likelihoods of becoming infected change). (B) The model captures key populations including sex workers and people who buy sex, drawn from the general population. The coloured arrows below indicate the three routes of transmission that have been modelled (sexual, nonsexual and zoonotic). (C) Vaccine stratum j = 1 accounts for individuals who have had previous smallpox
vaccination, with no one able to move in or out of j = 1. The remaining strata correspond to being unvaccinated (j = 2), vaccinated with a single dose of either MVA-BN or LC16m8 (j = 3) or vaccinated with two doses of MVA-BN (j = 4). The population may only move through the vaccine strata sequentially and cannot move backwards (e.g. going from being vaccinated to unvaccinated).

## Installation
This package relies heavily upon odin (a high-level language for implementing mathematical models) and monty (a package used to fit odin models to data). These packages require a compiler to install dependencies for the package and to build any models with odin. Windows users should install Rtools. Be sure to select the “edit PATH” checkbox during installation or the tools will not be found.

After installation of odin, ensure you have the devtools package installed by running the following:

install.packages("devtools")

Then install the mpoxseir package (https://github.com/mrc-ide/mpoxseir) directly from GitHub by running:

devtools::install_github("mrc-ide/mpoxseir")

If you have any problems installing then please raise an issue on the mpoxseir GitHub.

If everything has installed correctly, you then need to load the package:

library(mpoxseir)
