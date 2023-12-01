## Installation

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

Simulations are started using the function `calculate_e_transport(...)`. The function takes in many parameters, so it can be easier to use the script named `Control_script.jl` situated in the `scripts/` folder.
The script is pre-filled and you just need to modify the values of the parameters.
You can also use the script as a template to make your own control scripts.

After you have modified the `Control_script.jl` and saved it, you just have to 
execute it. This can be done from the Julia REPL or from the command line. 

### Starting simulation from the Julia REPL

!!! warning "Activating the AURORA environment"
    To be able to use AURORA.jl, the repository environment needs to be activated. This can be done for example by starting Julia from the *AURORA.jl/* folder using the command
    ```
    $> julia --project=.
    ```
    Or by entering the Pkg REPL using the `]` key and typing
    ```julia-repl
    pkg> activate .
    ```

!!! info "Using VS Code"
    If you are using VS Code with the Julia extension, the local environment should
    be automatically activated when you open the AURORA.jl/ folder.

Once the `AURORA.jl` environment activated, you can start simulations from the Julia REPL with the command
```julia-repl
julia> include("scripts/Control_script.jl")
```
If you are using VS Code, you can also use the "Execute active File in REPL" button.

The results will be saved in a folder under `data/` along with the parameters used to run the simulation.


### Starting simulation from the command line

Move to the *AURORA.jl/* folder. Then, execute the `Control_script.jl` using the
command 
```
$> julia --project=@. scripts/Control_script.jl 
```

The results will be saved in a folder under `data/` along with the parameters used to run the simulation.