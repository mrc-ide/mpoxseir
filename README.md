# mpoxseir

This package contains a compartmental SEIRD model of mpox disease transmission and progression. A schematic diagram of the model can be found below:

<img width="351" height="344" alt="image" src="https://github.com/user-attachments/assets/f79b3f36-a1c2-41b5-854d-823d250e6f3e" />

This package relies heavily upon odin (a high-level language for implementing mathematical models) and monty (a package used to fit odin models to data). These packages require a compiler to install dependencies for the package and to build any models with odin. Windows users should install Rtools. Be sure to select the “edit PATH” checkbox during installation or the tools will not be found.

After installation of odin, ensure you have the devtools package installed by running the following:

install.packages("devtools")

Then install the mpoxseir package (https://github.com/mrc-ide/mpoxseir) directly from GitHub by running:

devtools::install_github("mrc-ide/mpoxseir")

If you have any problems installing then please raise an issue on the mpoxseir GitHub.

If everything has installed correctly, you then need to load the package:

library(mpoxseir)
