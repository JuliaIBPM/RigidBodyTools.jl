module RigidBodyTools

using StaticArrays
using LinearAlgebra
using UnPack

import Base: *, +, -, inv, transpose, vec

import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect,view,findall

import LinearAlgebra: dot


export Body
export RigidBodyMotion, AbstractKinematics, d_dt, motion_velocity, motion_state,
          surface_velocity!, surface_velocity, update_body!,
          AbstractDeformationMotion, ConstantDeformationMotion, DeformationMotion,
          NullDeformationMotion,maxvelocity, zero_body, update_exogenous!,
          zero_exogenous, body_velocities, velocity_in_body_coordinates_2d

export parent_body_of_joint, child_body_of_joint, parent_joint_of_body, child_joints_of_body,
        position_vector,velocity_vector,deformation_vector, exogenous_position_vector,
        exogenous_velocity_vector, unconstrained_position_vector, unconstrained_velocity_vector

export AbstractDOFKinematics, AbstractPrescribedDOFKinematics, DOFKinematicData, SmoothRampDOF,
        OscillatoryDOF, ConstantVelocityDOF, CustomDOF, ExogenousDOF, ConstantPositionDOF,
        UnconstrainedDOF, Kinematics, dof_position, dof_velocity, dof_acceleration, ismoving,
        is_system_in_relative_motion

export RigidTransform, rotation_about_x, rotation_about_y, rotation_about_z, rotation_from_quaternion,
          quaternion, rotation_about_axis, rotation_identity,
          MotionTransform, ForceTransform, AbstractTransformOperator, rotation_transform,
          cross_matrix, cross_vector, translation, rotation, motion_subspace, joint_velocity

export PluckerForce, PluckerMotion, motion_rhs!, zero_joint, zero_motion_state, init_motion_state,
        angular_only, linear_only

export Joint, joint_transform, parent_to_child_transform, LinkedSystem, position_dimension,
        exogenous_dimension, constrained_dimension, unconstrained_dimension, position_and_vel_dimension,
        body_transforms, number_of_dofs

export RevoluteJoint, PrismaticJoint, CylindricalJoint, SphericalJoint, FreeJoint, FreeJoint2d,
        FixedJoint



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
abstract type AbstractDeformationMotion <: AbstractMotion end


numpts(::Body{N}) where {N} = N
numpts(::Nothing) = 0

include("kinematics.jl")
include("rigidtransform.jl")
include("joints.jl")

include("lists.jl")

include("rigidbodymotions.jl")
include("directmotions.jl")

include("tools.jl")
include("assignvelocity.jl")
include("shapes.jl")

include("plot_recipes.jl")

end
