## Installation

1. Clone the repository (e.g. with Git) or download and extract the .zip (available under the green _**code**_ button on the GitHub page)

2. Open a terminal, go in the folder where the code is installed, and start Julia with the command:
```
$> julia
```

3. Then, **activate** the repository and download the packages required by *AURORA.jl* using the commands:
```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate() # this might take a while...
```

4. *AURORA.jl* is now ready to use!

## Folder structure
The code is structured as follow
```
AURORA/
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
The folder `data/` contains the subfolders where simulation results are saved.

The folder `docs/` contains all the necessary scripts to power this documentation.

The folder `e_cascading_data/` is where the cascading data produced by the simulations are saved for future use by the program itself. The cascading data are saved by species in the subfolders `N2/`, `O2/` and `O/`. You should not need to venture into this folder.

The folder `scripts/` contains the scripts for the user to start the simulations and plot some of the results.

The folder `src/` contains the source code of the model.

## Get started
!!! warning "Activating the AURORA environment"
    To be able to use AURORA.jl, the repository environment needs to be activated. This can be done for example by starting Julia from the *AURORA/* folder using the command
    ```
    $> julia --project=.
    ```

The idea is to use  using `Control_script.jl` to control the simulations. You can open the script, change the parameters of the simulation, and save it. Then, you can start the simulation with the command:
```julia-repl
julia> include("scripts/Control_script.jl")
```
The results will be saved in a folder under `data/` along with the parameters used to run the simulation.


