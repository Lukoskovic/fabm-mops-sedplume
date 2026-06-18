This repository is a fork ( *commit*  03fd736382e1b4a0b10a2193e7e48145a81383d8) from the FABM implementation of MOPS.  

Original `README.md`
> This is an implementation of the biogeochemical model [Model of Oceanic Pelagic Stoichiometry (MOPS)](https://doi.org/10.5194/gmd-8-2929-2015) for [FABM](https://fabm.net).
> 
> This work was funded by the European Union under grant agreement no. 101083922 (OceanICU) and UK Research and Innovation (UKRI) under the UK government’s Horizon Europe funding guarantee [grant number 10054454, 10063673, 10064020, 10059241, 10079684, 10059012, 10048179]. The views, opinions and practices used to produce this dataset/software are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

## Differences to the previous Implementation

This repository contains additional implementations, which are not per default part of MOPS. 

### Sediment Module

The main addition is the implementation of the Sediment module `sediment.F90`.