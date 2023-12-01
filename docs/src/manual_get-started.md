# Get started

Simulations are started using the function `calculate_e_transport(...)`. The function takes in many parameters, so it can be easier to use the script named `Control_script.jl` situated in the `scripts/` folder.
The script is pre-filled and you just need to modify the values of the parameters.
You can also use the script as a template to make your own control scripts.

After you have modified the `Control_script.jl` and saved it, you just have to 
execute it. This can be done from the Julia REPL or from the command line. 

## Starting simulation from the Julia REPL

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


## Starting simulation from the command line

Move to the *AURORA.jl/* folder. Then, execute the `Control_script.jl` using the
command 
```
$> julia --project=@. scripts/Control_script.jl 
```

The results will be saved in a folder under `data/` along with the parameters used to run the simulation.