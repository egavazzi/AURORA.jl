# AURORA.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://egavazzi.github.io/AURORA.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://egavazzi.github.io/AURORA.jl/dev)
[![Test workflow status](https://github.com/egavazzi/AURORA.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/egavazzi/AURORA.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/egavazzi/AURORA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/egavazzi/AURORA.jl)
[![Lint workflow Status](https://github.com/egavazzi/AURORA.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/egavazzi/AURORA.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/egavazzi/AURORA.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/egavazzi/AURORA.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11238620.svg)](https://doi.org/10.5281/zenodo.1123862)


AURORA is a time-dependent multi-stream electron transport code.

It is suitable for modeling the ionospheric electron fluxes when the electron precipitation varies rapidly, typically on sub-second timescales.

Below is an example of modeled ionospheric electron fluxes in response to a suprathermal electron burst measured by a rocket.

https://github.com/user-attachments/assets/e67fc4d3-1fe1-4c0b-a8f6-42275f6e4659



## Installation / Documentation
Instructions are available in the [**documentation**](https://egavazzi.github.io/AURORA.jl/dev/).

Descriptions of how the code works are available in Gustavsson (2022), in section 3 of Gavazzi (2022), and in [this document](https://github.com/egavazzi/AURORA.jl/blob/main/docs/other/AURORA_Documentation.pdf).



## Citation
If you use AURORA.jl in your work, please cite using the reference given in the *Cite this repository* button in the **About** section of this page. Note that the DOI ![10.5281/zenodo.11238620](https://doi.org/10.5281/zenodo.11238620) refers to the general project and all its versions. Version-specific DOIs can be found on Zenodo.



## References
Gavazzi, E. (2022). The effects of time-variation of electron fluxes from the auroral ionosphere on M-I coupling [Master thesis, UiT Norges arktiske universitet]. https://munin.uit.no/handle/10037/25897

Gustavsson, B. (2022). Time-Dependent Electron Transport I: Modelling of Supra-Thermal Electron Bursts Modulated at 5–10 Hz With Implications for Flickering Aurora. Journal of Geophysical Research: Space Physics, 127(6), e2019JA027608. https://doi.org/10.1029/2019JA027608
