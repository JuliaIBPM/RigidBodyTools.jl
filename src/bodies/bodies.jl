export Body

abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type Body{N,C<:BodyClosureType} end

numpts(::Body{N}) where {N} = N
numpts(::Nothing) = 0

include("rigidbodymotions.jl")
include("rigidtransform.jl")
include("bodylist.jl")


include("tools.jl")
include("assignvelocity.jl")
include("shapes.jl")
