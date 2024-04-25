using Documenter
using BoundaryWall
using StaticArrays
# using DocumenterTools: Themes

makedocs(
         sitename = "BoundaryWall.jl",
         authors  = "Alberto Ruiz-Biestro",
         format = Documenter.HTML(
                    prettyurls = get(ENV, "CI", nothing) == "true",
                    mathengine=KaTeX(
                         Dict(:macros => Dict(
                            "\\x" => "\\boldsymbol{x}"
                                             )
                             )
                            )
                    ),

         pages = Any[
                     "Home" => "index.md",
                     "Reference guide" => [
                                           "Geometry" => "geometry.md",
                                           "Incident waves" => "incident.md",
                                           "Lattices" => "grids.md"
                                          ],
                     "API" => "api.md"

                    ]
        )

