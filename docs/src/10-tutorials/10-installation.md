# Installation

1. Clone the repository (e.g. with Git) or download and extract the .zip (available under the green _**code**_ button on the GitHub page).

2. Open a terminal, move into the folder where AURORA was just installed, and start Julia with the command
   ```
   $> julia --project=.
   ```

3. We now need to install the packages required by AURORA. To do this, we need to **instantiate** the repository environment using the following command
   ```julia-repl
   julia> using Pkg
   julia> Pkg.instantiate() # this might take a while the first time...
   ```

4. AURORA is now ready to use! See the [Getting started](@ref) tutorial.


!!! note
    - "Instantiating" the environment will download and install dependencies in their correct version. This needs to be done only the first time you use AURORA.
    - The `--project=.` flag **activates** the local environment. This has to be done every time you start a new Julia session and want to use AURORA.

!!! tip "Using VS Code"
    If you are using VS Code with the Julia extension, the local environment should be automatically activated when you open the `AURORA.jl/` folder.

!!! warning "Supported OS"
	AURORA is developed on Linux but should also work on macOS and Windows platforms. If you encounter problems, open an issue on the [GitHub repository](https://github.com/egavazzi/AURORA.jl/issues) or contact [etienne.gavazzi@uit.no](mailto:etienne.gavazzi@uit.no).

!!! info "Python dependencies (MSIS & IRI)"
	AURORA uses `pymsis` and `iri2020` Python packages (via CondaPkg) for atmospheric background models. These packages depend on external libraries, so on some systems they can occasionally run into build or compilation issues. 
	If that happens, see the [Troubleshooting](@ref) page.