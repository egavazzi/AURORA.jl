#=
This script is to create the documentation locally on the user's machine.
To use it, move into the docs/ folder, and type the command > include("make_local.jl")
=#

# push!(LOAD_PATH, "../src/")
# using Pkg
# Pkg.activate("../")

using AURORA
using Documenter

makedocs(
        sitename = "AURORA.jl",
        modules = [AURORA],
        pages = [
                "Home" => "index.md",
                "Manual" => [
                    "Installation" => "manual_installation.md",
                    "Get started" => "manual_get-started.md",
                    "Folder structure" => "manual_folder-structure.md"
                    ],
                ],
        format = Documenter.HTML(prettyurls = false),
        warnonly = true)
