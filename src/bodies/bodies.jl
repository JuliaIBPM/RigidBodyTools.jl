export Body

abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type Body{N,C<:BodyClosureType} end

include("rigidbodymotions.jl")
include("rigidtransform.jl")
include("bodylist.jl")

include("tools.jl")
include("shapes.jl")
