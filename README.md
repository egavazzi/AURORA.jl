# AURORA.jl

| **Documentation**                 | **DOI**                 |
|:---------------------------------:|:-----------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][doi-img]][doi-url] |

This is a julia implementation of the following time-dependent model of electron transport in the ionosphere : 
https://github.com/egavazzi/AURORA



## Installation

1. Clone the repository (e.g. with Git) or download and extract the .zip (available under the green _**code**_ button)

2. Open a terminal, **go in the folder where the code is installed**, and start Julia with the command:
```
$> julia
```

3. Then, **activate** the repository and download the packages required by *AURORA.jl* using the commands:
```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

4. *AURORA.jl* is now ready to use!



## Documentation
The documentation is available [**here**](https://egavazzi.github.io/AURORA.jl/dev/).





[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://egavazzi.github.io/AURORA.jl/dev/
[doi-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.11238620.svg
[doi-url]: https://doi.org/10.5281/zenodo.11238620
