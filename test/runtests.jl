using RigidBodyTools
using Literate
using Test

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

const GROUP = get(ENV, "GROUP", "All")

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"

for (root, dirs, files) in walkdir(litdir)
    if splitpath(root)[end] == "assets"
        for file in files
            cp(joinpath(root, file),joinpath(notebookdir,file),force=true)
            cp(joinpath(root, file),joinpath(docdir,file),force=true)
            cp(joinpath(root, file),joinpath(".",file),force=true)
        end
    end
end

if GROUP == "Transforms"
  include("transforms.jl")
end

if GROUP == "Bodies"
  include("bodies.jl")
end

if GROUP == "All" || GROUP == "Auxiliary"
  include("transforms.jl")
  include("bodies.jl")
end

if GROUP == "All" || GROUP == "Literate"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && @testset "$file" begin include(joinpath(root,file)) end
    end
  end
end

if GROUP == "Notebooks"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      #endswith(file,".jl") && startswith(file,"defor") && Literate.notebook(joinpath(root, file),notebookdir)
      endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir)
    end
  end
end

if GROUP == "Documentation"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir)
    end
  end
end
