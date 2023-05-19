#=
This script is to create the documentation locally on the user's machine.
To use it, move into the docs/ folder, and type the command > include("make_local.jl")
=#

push!(LOAD_PATH, "../src/")
using Pkg
Pkg.activate("../")

using AURORA
using Documenter

makedocs(
        sitename="AURORA.jl",
        modules =[AURORA],
        pages = [
                "Home" => "index.md",
                "Manual" => "manual.md"
                ],
        format = Documenter.HTML(prettyurls = false))
