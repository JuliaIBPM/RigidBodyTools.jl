using RigidBodyTools
using Literate
using Test

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

const GROUP = get(ENV, "GROUP", "All")

macro mysafetestset(args...)
    name, expr = args
    quote
        ex = quote
          name_str = $$(QuoteNode(name))
          expr_str = $$(QuoteNode(expr))
          mod = gensym(name_str)
          ex2 = quote
              @eval module $mod
                      using Test
                      @testset $name_str $expr_str
                    end
              nothing
          end
          eval(ex2)
        end
        eval(ex)
    end
end

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
      global file_str = "$file"
      global body = :(begin include(joinpath($root,$file)) end)
      endswith(file,".jl") && @mysafetestset file_str body
      #endswith(file,".jl") && @testset "$file" begin include(joinpath(root,file)) end
    end
  end
end

if GROUP == "Notebooks"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      #endswith(file,".jl") && startswith(file,"exog") && Literate.notebook(joinpath(root, file),notebookdir)
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
