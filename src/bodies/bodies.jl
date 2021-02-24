export Body

const NDIM = 2
const CHUNK = 3*(NDIM-1)


abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type Body{N,C<:BodyClosureType} end

numpts(::Body{N}) where {N} = N
numpts(::Nothing) = 0

include("rigidbodymotions.jl")
include("directmotions.jl")
include("motionprofiles.jl")
include("kinematics.jl")

include("rigidtransform.jl")
include("bodylist.jl")


include("tools.jl")
include("assignvelocity.jl")
include("shapes.jl")
