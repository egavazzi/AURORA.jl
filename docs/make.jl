push!(LOAD_PATH, "../src/")

using Aurora
using Documenter

makedocs(
        sitename="Aurora.jl",
        modules =[Aurora],
        pages = [
                "Home" => "index.md"
                ])

# deploydocs(;
#     repo="github.com/egavazzi/Aurora.jl",
#     )