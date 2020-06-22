export Body

abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type Body{N,C<:BodyClosureType} end

include("bodylist.jl")
include("rigidtransform.jl")
include("tools.jl")
include("shapes.jl")
