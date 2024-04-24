using Documenter,BoundaryWall
# using DocumenterTools: Themes

makedocs(
         sitename = "BoundaryWall",
         authors  = "Alberto Ruiz-Biestro",
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true"
                                 ),
         pages = Any[
                     "Home" => "index.md"
                     "Reference guide" => [
                                           "Incident waves" => "incident.md"
                                          ]
                    ]
        )

