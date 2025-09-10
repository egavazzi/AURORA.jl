#=
This script is to create the documentation locally on the user's machine.
To use it, move into the `docs/` folder, activate the local environment, instantiate it
(if it has never been done before) and run the script.

Example of commands to run in the `docs/` folder:
```
using Pkg
Pkg.activate(".")
Pkg.instantiate() # you will have to run this if you haven't done it before
include("make_local.jl")
```
=#

using AURORA
using Documenter

makedocs(sitename = "AURORA.jl",
         modules = [AURORA],
         pages = [
             "Home" => "index.md",
             "Manual" => [
                 "Installation" => "manual_installation.md",
                 "Get started" => "manual_get-started.md",
                 "Folder structure" => "manual_folder-structure.md"
             ],
             "Troubleshooting" => "troubleshooting.md",
             "API" => "api.md"
         ],
         format = Documenter.HTML(prettyurls = false),
         warnonly = true)
