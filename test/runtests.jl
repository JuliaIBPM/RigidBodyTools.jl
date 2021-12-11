using RigidBodyTools
using Literate
using Test


const GROUP = get(ENV, "GROUP", "All")

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"

if GROUP == "All" || GROUP == "Auxiliary"
  include("motions.jl")
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
