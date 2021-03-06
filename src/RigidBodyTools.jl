module RigidBodyTools


export Body

const NDIM = 2
const CHUNK = 3*(NDIM-1)


abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type PointShiftType end
abstract type Unshifted <: PointShiftType end
abstract type Shifted <: PointShiftType end

abstract type Body{N,C<:BodyClosureType} end

numpts(::Body{N}) where {N} = N
numpts(::Nothing) = 0

include("rigidbodymotions.jl")
include("directmotions.jl")
include("kinematics.jl")

include("rigidtransform.jl")
include("bodylist.jl")


include("tools.jl")
include("assignvelocity.jl")
include("shapes.jl")

include("plot_recipes.jl")

end
