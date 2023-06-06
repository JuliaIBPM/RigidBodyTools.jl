### JOINTS ###


abstract type AbstractJointType end


# Joint structure

struct Joint{ND,JT}

    "Parent body ID (0 for inertial)"
    parent_id :: Int

    "Child body ID"
    child_id :: Int

    "Constrained DOFs with prescribed behaviors"
    cdofs :: Vector{Int}

    "Constrained DOFs with exogenously-specified behaviors"
    edofs :: Vector{Int}

    "Unconstrained DOFS"
    udofs :: Vector{Int}

    "Kinematics for prescribed DOFs"
    kins :: Vector{AbstractPrescribedDOFKinematics}

    "Initial values for DOFs"
    xinit_dofs :: Vector{Float64}

    "Joint parameters"
    params :: Dict

    "Transform matrix from parent to joint"
    Xp_to_j :: MotionTransform

    "Transform matrix from child to joint"
    Xch_to_j :: MotionTransform
end



function Joint(::Type{JT},parent_id::Int,child_id::Int,Xp_to_j::MotionTransform{ND},Xch_to_j::MotionTransform{ND},dofs::Vector{T};
                  params = Dict(), xinit = zeros(Float64,3(ND-1))) where {JT<:AbstractJointType,ND,T<:AbstractDOFKinematics}
    @assert length(dofs)==number_of_dofs(JT) "input dofs must have correct length"
    dof_lists = Dict("constrained" => Int[], "exogenous" => Int[], "unconstrained" => Int[])
    kins = AbstractDOFKinematics[]

    # Note - should only need to provide enough dof kinematics for the *actual* dofs
    for (i,dof) in enumerate(dofs)
        _classify_joint!(dof_lists,kins,i,dof)
    end
    Joint{ND,JT}(parent_id,child_id,dof_lists["constrained"],dof_lists["exogenous"],dof_lists["unconstrained"],
            kins,xinit,params,Xp_to_j,Xch_to_j)
end


function check_q_dimension(q,C::Type{T}) where T<:AbstractJointType
  @assert length(q) == state_dimension(C) "Incorrect length of q"
end

function joint_transform(q::AbstractVector,joint::Joint{ND,JT}) where {ND,JT<:AbstractJointType}
  @unpack params = joint
  check_q_dimension(q,JT)
  joint_transform(q,JT,params,Val(ND))
end

function parent_to_child_transform(q::AbstractVector,joint::Joint)
  @unpack Xp_to_j, Xch_to_j = joint
  Xj = joint_transform(q,joint)
  Xp_to_ch = inv(Xch_to_j)*Xj*Xp_to_j
end


function _classify_joint!(dof_lists,kins,index,kin::AbstractDOFKinematics)
    push!(dof_lists["constrained"],index)
    push!(kins,kin)
    nothing
end

_classify_joint!(dof_lists,kins,index,::ExogenousDOF) = (push!(dof_lists["exogenous"],index); nothing)
_classify_joint!(dof_lists,kins,index,::UnconstrainedDOF) = (push!(dof_lists["unconstrained"],index); nothing)



# Joint types


abstract type RevoluteJoint <: AbstractJointType end

number_of_dofs(::Type{RevoluteJoint}) = 1
state_dimension(::Type{RevoluteJoint}) = 1

function joint_transform(q::AbstractVector,::Type{RevoluteJoint},p::Dict,::Val{ND}) where {ND}
  R = rotation_about_z(-q[1])
  x = SVector{3}([0.0,0.0,0.0])
  return MotionTransform{ND}(x,R)
end

abstract type PrismaticJoint <: AbstractJointType end

number_of_dofs(::Type{PrismaticJoint}) = 1
state_dimension(::Type{PrismaticJoint}) = 1

function joint_transform(q::AbstractVector,::Type{PrismaticJoint},p::Dict,::Val{3})
  R = rotation_identity()
  x = SVector{3}([0.0,0.0,q[1]])
  return MotionTransform{3}(x,R)
end

abstract type HelicalJoint <: AbstractJointType end

number_of_dofs(::Type{HelicalJoint}) = 1
state_dimension(::Type{HelicalJoint}) = 1

function joint_transform(q::AbstractVector,C::Type{HelicalJoint},p::Dict,::Val{3})
  R = rotation_about_z(-q[1])
  x = SVector{3}([0.0,0.0,p["pitch"]*q[1]])
  return MotionTransform{3}(x,R)
end

abstract type CylindricalJoint <: AbstractJointType end

number_of_dofs(::Type{CylindricalJoint}) = 2
state_dimension(::Type{CylindricalJoint}) = 2

function joint_transform(q::AbstractVector,::Type{CylindricalJoint},p::Dict,::Val{3})
  R = rotation_about_z(-q[1])
  x = SVector{3}([0.0,0.0,q[2]])
  return MotionTransform{3}(x,R)
end

abstract type SphericalJoint <: AbstractJointType end

number_of_dofs(::Type{SphericalJoint}) = 3
state_dimension(::Type{SphericalJoint}) = 4

function joint_transform(q::AbstractVector,::Type{SphericalJoint},p::Dict,::Val{3})
  R = rotation_from_quaternion(quaternion(q))
  x = SVector{3}([0.0,0.0,0.0])
  return MotionTransform{3}(x,R)
end

abstract type FreeJoint <: AbstractJointType end

number_of_dofs(::Type{FreeJoint}) = 6
state_dimension(::Type{FreeJoint}) = 7

function joint_transform(q::AbstractVector,::Type{FreeJoint},p::Dict,::Val{3})
  R = rotation_from_quaternion(quaternion(q[1:4]))
  x = SVector{3}(q[5:7])
  return MotionTransform{3}(R'*x,R)
end

abstract type FreeJoint2d <: AbstractJointType end

number_of_dofs(::Type{FreeJoint2d}) = 3
state_dimension(::Type{FreeJoint2d}) = 3

function joint_transform(q::AbstractVector,::Type{FreeJoint2d},p::Dict,::Val{2})
  R = rotation_about_z(-q[1])
  x = SVector{3}(q[2],q[3],0.0)
  return MotionTransform{2}(x,R)
end

###
