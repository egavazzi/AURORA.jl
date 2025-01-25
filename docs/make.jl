# push!(LOAD_PATH, "../src/")

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
                "Troubleshooting" => "troubleshooting.md",
                ],
        warnonly = true)

deploydocs(;
    repo="github.com/egavazzi/AURORA.jl.git",
    )
