## Installation

1. Clone the repository (e.g. with Git) or download and extract the .zip (available under the green _**code**_ button on the GitHub page)

2. Open a terminal, go in the folder where the code is installed, and start Julia with the command:
```
$> julia
```

3. Then, **activate** the repository and download the packages required by *Aurora.jl* using the commands:
```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

4. *Aurora.jl* is now ready to use!

## Folder structure
The code is structured as follow
```
Aurora/
├── data/
│   └── 20220926/...
│   └── ...
│
├── docs/...
│ 
├── e_cascading_data/
│   └── N2/...
│   └── O2/...
│   └── O/...
│
├── scripts/
│   └── Control_script.jl
│   └── ...
│ 
└── src/
    └── main.jl
    └── cascading.jl
    └── setup.jl
    └── ...

```