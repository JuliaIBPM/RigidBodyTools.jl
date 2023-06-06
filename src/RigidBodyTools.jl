module RigidBodyTools

using StaticArrays
using LinearAlgebra
using UnPack

import Base: *, inv, transpose, vec


export Body
export RigidBodyMotion, AbstractKinematics, d_dt, motion_velocity, motion_state,
          surface_velocity!, surface_velocity, update_body!,
          AbstractDeformationMotion, ConstantDeformationMotion, DeformationMotion,
          RigidAndDeformingMotion,maxvelocity, maxlistvelocity

export AbstractDOFKinematics, AbstractPrescribedDOFKinematics, DOFKinematicData, SmoothRampDOF,
        OscillatoryDOF, ConstantVelocityDOF, CustomDOF, ExogenousDOF, ConstantStateDOF,
        UnconstrainedDOF, Kinematics

export RigidTransform, rotation_about_x, rotation_about_y, rotation_about_z, rotation_from_quaternion,
          quaternion, rotation_about_axis, rotation_identity,
          MotionTransform, ForceTransform, AbstractTransformOperator,
          cross_matrix, cross_vector, translation, rotation

export Joint, joint_transform, parent_to_child_transform

export Oscillation, OscillationX, OscillationY, OscillationXY, RotationalOscillation,
        PitchHeave, Pitchup, EldredgeRamp, ColoniusRamp, SwitchedKinematics,
        complex_translational_position, complex_translational_velocity, complex_translational_acceleration,
        angular_position, angular_velocity, angular_acceleration,
        translational_position, translational_velocity, translational_acceleration,
        KinematicData



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
include("joints.jl")

include("lists.jl")

include("tools.jl")
include("assignvelocity.jl")
include("shapes.jl")

include("plot_recipes.jl")

end
