# Installation

1. Clone the repository (e.g. with Git) or download and extract the .zip (available under the green _**code**_ button on the GitHub page).

2. Open a terminal, move into the folder where AURORA was just installed, and start Julia with the command
```
$> julia
```

3. We now need to install the packages required by AURORA. To do this, we need to **activate** and **instantiate** the repository environment using the following commands
```julia-repl
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate() # this might take a while...
```

4. AURORA is now ready to use!

<br></br>
<br></br>

**Notes:**
- "Instantiating" the environment will download and install dependencies in their correct version. This needs to be done only the first time you use AURORA.
- "Activating" the local environment will however have to be done everytime you start a new Julia session and want to use AURORA.
