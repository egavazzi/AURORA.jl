#=
To run locally:
- Install LiveServer.jl in your global environment
- Move to the AURORA folder
- Run
    ```
    $> julia --project=docs -e 'using AURORA, LiveServer; servedocs()`
    ```
=#

using AURORA
using Documenter

DocMeta.setdocmeta!(
    AURORA,
    :DocTestSetup,
    :(using AURORA);
    recursive=true,
    )

makedocs(
        sitename = "AURORA.jl",
        modules = [AURORA],
        pages = [
                "Home" => "index.md",
                "Manual" => [
                    "Installation" => "manual_installation.md",
                    "Get started" => "manual_get-started.md",
                    "Folder structure" => "manual_folder-structure.md",
                    "Visualization" => "manual_visualization.md"
                    ],
                "Troubleshooting" => "troubleshooting.md",
                "API" => "api.md"
                ],
        warnonly = false,
        checkdocs = :none,
        )

deploydocs(;
    repo="github.com/egavazzi/AURORA.jl.git",
    )
