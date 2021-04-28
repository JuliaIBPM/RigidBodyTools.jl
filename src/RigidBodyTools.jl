module RigidBodyTools


export Body

const NDIM = 2
const CHUNK = 3*(NDIM-1)


abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type Body{N,C<:BodyClosureType} end

numpts(::Body{N}) where {N} = N
numpts(::Nothing) = 0

include("bodies/rigidbodymotions.jl")
include("bodies/directmotions.jl")
include("bodies/kinematics.jl")

include("bodies/rigidtransform.jl")
include("bodies/bodylist.jl")


include("bodies/tools.jl")
include("bodies/assignvelocity.jl")
include("bodies/shapes.jl")

include("plot_recipes.jl")

end
