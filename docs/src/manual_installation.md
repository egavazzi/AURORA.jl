# Installation

1. Clone the repository (e.g. with Git) or download and extract the .zip (available under the green _**code**_ button on the GitHub page)

2. Open a terminal, move into the AURORA.jl, and start Julia with the command
```
$> julia
```

3. Then, **activate** the repository and install the packages required by *AURORA.jl* using the commands
```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate() # this might take a while...
```

4. *AURORA.jl* is now ready to use!