# AURORA.jl

| **Documentation**                 | **DOI**                 |
|:---------------------------------:|:-----------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][doi-img]][doi-url] |

AURORA is a time-dependent multi-stream electron transport code, suitable for modeling ionospheric electron-fluxes during periods of rapidly varying electron-precipitation.


This is a Julia implementation of the original AURORA code written in MATLAB and available at https://github.com/egavazzi/AURORA. The version present here is the one we recommend to use. It is in active development, is faster (~1000x), and produces more accurate results (i.e. bugs have been fixed).

Descriptions of the code are available in Gustavsson (2022) and Gavazzi (2022).

## Installation
Instructions are available in the documentation.

## Documentation
The documentation is available [**here**](https://egavazzi.github.io/AURORA.jl/dev/).

## References
Gavazzi, E. (2022). The effects of time-variation of electron fluxes from the auroral ionosphere on M-I coupling [Master thesis, UiT Norges arktiske universitet]. https://munin.uit.no/handle/10037/25897

Gustavsson, B. (2022). Time-Dependent Electron Transport I: Modelling of Supra-Thermal Electron Bursts Modulated at 5–10 Hz With Implications for Flickering Aurora. Journal of Geophysical Research: Space Physics, 127(6), e2019JA027608. https://doi.org/10.1029/2019JA027608




[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://egavazzi.github.io/AURORA.jl/dev/
[doi-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.11238620.svg
[doi-url]: https://doi.org/10.5281/zenodo.11238620
