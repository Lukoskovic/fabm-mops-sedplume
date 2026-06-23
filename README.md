This repository is a fork (*commit*  03fd736382e1b4a0b10a2193e7e48145a81383d8) from the FABM implementation of MOPS [^1].  

Original `README.md`
> This is an implementation of the biogeochemical model [Model of Oceanic Pelagic Stoichiometry (MOPS)](https://doi.org/10.5194/gmd-8-2929-2015) for [FABM](https://fabm.net).
> 
> This work was funded by the European Union under grant agreement no. 101083922 (OceanICU) and UK Research and Innovation (UKRI) under the UK government’s Horizon Europe funding guarantee [grant number 10054454, 10063673, 10064020, 10059241, 10079684, 10059012, 10048179]. The views, opinions and practices used to produce this dataset/software are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

FABM-MOPS model implementation is based on the MOPS code that was published [here](https://hdl.handle.net/20.500.12085/b174de1c-0bed-47f5-9718-7a8d44d1d2d1)[^2], which is the supplement by Kriest et al. (2023)[^3].

## Differences to the previous Implementation

This repository contains additional implementations, which are not default in FABM-MOPS. 

### Sediment Module

The main addition is the implementation of the Sediment module `sediment.F90`.

A new variable "sediment" was created in order to simulate sediment plume that is transported as a tracer and contains a sinking velocity. By taking the inspiration from the detritus, the module for the sediment is structured accordingly. However, the main difference is here, that sinking of the sediment is not following a Martin b curve, which could be still activated with the assumption of sediment remineralization (e.g. by contents of organic matter). Instead, there is a constant depth-independent sinking velocity.

An inherent optical property of absorption and scattering of sediment particles leads to an increase of light attenuation. Therefore, a light attenuation coefficient `ACksed` was added. 

### Sediment Impact Processes

Sensitivity experiments of sediment impacts are possible to carry out by using 4 implemented factors of MOPS Zooplankton parameters (default are unchanged 1). By changing the factors \< 1 \< one can test the system responses of potential impacts.  
Following processes are included:
- changed (reduced) grazing rate as in `ACmuzoofac`
- changed (reduced) assimilation rate `ACefffac`
- changed (reduced) densitiy-dependent loss rate (`AComnizfac`)
- changed (increased) zooplankton mortality (`zlambdafac`)

*Remark*: `...fac` points out to the additional factor of the corresponding parameter that is given in the name before the suffix.


# Sources

[^1]: FABM-MOPS repository is [https://github.com/BoldingBruggeman/fabm-mops](https://github.com/BoldingBruggeman/fabm-mops) 
[^2]: Kriest, Iris, Getzlaff, Julia, Landolfi, Angela, Sauerland, Volkmar, Schartau, Markus, and Oschlies, Andreas (2023). Supplemental dataset to Kriest et al. (2023): Exploring the role of different data types and timescales for the quality of marine biogeochemical model calibration [dataset]. GEOMAR Helmholtz Centre for Ocean Research Kiel [distributor]. [hdl:20.500.12085/b174de1c-0bed-47f5-9718-7a8d44d1d2d1](https://hdl.handle.net/20.500.12085/b174de1c-0bed-47f5-9718-7a8d44d1d2d1)
[^3]: Kriest, I., Getzlaff, J., Landolfi, A., Sauerland, V., Schartau, M., and Oschlies, A.: Exploring the role of different data types and timescales in the quality of marine biogeochemical model calibration, Biogeosciences, 20, 2645–2669, [https://doi.org/10.5194/bg-20-2645-2023](https://doi.org/10.5194/bg-20-2645-2023), 2023. 
