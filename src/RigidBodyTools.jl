module RigidBodyTools

import Base: vec


export Body
export RigidBodyMotion, Kinematics, d_dt, motion_velocity, motion_state,
          surface_velocity!, surface_velocity
export Oscillation, OscillationX, OscillationY, OscillationXY, RotationalOscillation,
        PitchHeave, Pitchup, EldredgeRamp, ColoniusRamp


const NDIM = 2
const CHUNK = 3*(NDIM-1)



abstract type BodyClosureType end
abstract type OpenBody <: BodyClosureType end
abstract type ClosedBody <: BodyClosureType end

abstract type PointShiftType end
abstract type Unshifted <: PointShiftType end
abstract type Shifted <: PointShiftType end

abstract type Body{N,C<:BodyClosureType} end

abstract type AbstractMotion end


numpts(::Body{N}) where {N} = N
numpts(::Nothing) = 0

include("kinematics.jl")
include("rigidbodymotions.jl")
include("directmotions.jl")

include("rigidtransform.jl")
include("lists.jl")


include("tools.jl")
include("assignvelocity.jl")
include("shapes.jl")

include("plot_recipes.jl")

end
