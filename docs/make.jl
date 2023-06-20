using Documenter, RigidBodyTools

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

makedocs(
    sitename = "RigidBodyTools.jl",
    doctest = true,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Bodies and Transforms" => ["manual/shapes.md",
                                    "manual/transforms.md",
                                    "manual/bodylists.md"],
        "Motions" => ["manual/dofmotions.md",
                      "manual/joints.md",
                      "manual/exogenous.md"]
        #"Internals" => [ "internals/properties.md"]
    ],
    #format = Documenter.HTML(assets = ["assets/custom.css"])
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
    #assets = ["assets/custom.css"],
    #strict = true
)


#if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
     repo = "github.com/JuliaIBPM/RigidBodyTools.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
#end
