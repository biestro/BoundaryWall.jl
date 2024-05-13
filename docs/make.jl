using Documenter
using BoundaryWall
using StaticArrays
# using DocumenterTools: Themes

# Themes.compile(joinpath(@__DIR__, "src/assets/dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

makedocs(

sitename = "BoundaryWall.jl",
authors  = "Alberto Ruiz-Biestro",
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine=KaTeX(
                Dict(:macros => Dict("\\x" => "\\boldsymbol{x}"))
                ),
        assets = ["assets/custom.css",],
        ),

pages = Any[
        "Home" => "index.md",
        "Reference guide" => [
                "Geometry" => "geometry.md",
                "Incident waves" => "incident.md",
                "Lattices" => "grids.md"
                ],
        "Tutorials" => [
                "Beam splitters" =>"beamsplitter.md", 
                "Quantum Dots" => "billiards.md",
                "Vector fields" => "vectorial.md"
                ],
        "API" => "api.md"

                ]
)

