### JOINTS ###
#=
Note that all 3-d joints with Revolute, Prismatic, Helical, Cylindrical
joints all need the joint to be along its z axis, so this dictates
how the parent-to-joint and child-to-joint transforms are to be created.
=#

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

    #"Initial values for DOFs"
    #xinit_dofs :: Vector{Float64}

    "Joint parameters"
    params :: Dict

    "Buffer for holding joint dof velocity"
    vbuf :: Vector{Float64}

    "Transform matrix from parent to joint"
    Xp_to_j :: MotionTransform

    "Transform matrix from child to joint"
    Xch_to_j :: MotionTransform
end

"""
    Joint(jtype::AbstractJointType,parent_id,Xp_to_j::MotionTransform,child_id,Xch_to_j::MotionTransform,
            dofs::Vector{AbstractDOFKinematics};[params=Dict()])

Construct a joint of type `jtype`, connecting parent body of id `parent_id` to child body
of id `child_id`. The placement of the joint on the bodies is given by `Xp_to_j`
and `Xch_to_j`, the transforms from the parent and child coordinate systems to the joint
system, respectively. Each of the degrees of freedom for the joint is specified with
`dofs`, a vector of `AbstractDOFKinematics`. Any parameters required by the joint
can be passed along in the optional `params` dictionary.
"""
function Joint(::Type{JT},parent_id::Int,Xp_to_j::MotionTransform{ND},child_id::Int,Xch_to_j::MotionTransform{ND},dofs::Vector{T};
                  params = Dict()) where {JT<:AbstractJointType,ND,T<:AbstractDOFKinematics}
    @assert length(dofs)==number_of_dofs(JT) "input dofs must have correct length"
    dof_lists = Dict("constrained" => Int[], "exogenous" => Int[], "unconstrained" => Int[])
    kins = AbstractDOFKinematics[]

    for (i,dof) in enumerate(dofs)
        _classify_joint!(dof_lists,kins,i,dof)
    end

    vbuf = zeros(Float64,number_of_dofs(JT))

    Joint{ND,JT}(parent_id,child_id,dof_lists["constrained"],dof_lists["exogenous"],dof_lists["unconstrained"],
            kins,params,vbuf,Xp_to_j,Xch_to_j)
end



function Base.show(io::IO, joint::Joint{ND,JT}) where {ND,JT}
    println(io, "Joint of dimension $ND and type $JT")
    println(io, "   Constrained dofs = $(joint.cdofs)")
    println(io, "   Exogenous dofs = $(joint.edofs)")
    println(io, "   Unconstrained dofs = $(joint.udofs)")
end

physical_dimension(j::Joint{ND}) where {ND} = ND

function physical_dimension(jvec::Vector{<:Joint{ND}}) where {ND}
  #@assert allequal(physical_dimension.(jvec)) "Not all joints are same physical dimension"
  return physical_dimension(first(jvec))
end


position_dimension(j::Joint{ND,JT}) where {ND,JT} = position_dimension(JT)
number_of_dofs(j::Joint{ND,JT}) where {ND,JT} = number_of_dofs(JT)

function ismoving(joint::Joint)
    @unpack kins = joint
    any(map(dof -> ismoving(dof),kins))
end


"""
    motion_subspace(j::Joint)

Return the `6 x ndof` (3d) or `3 x ndof` (2d) matrix providing the mapping from joint dof velocities
to the full Plucker velocity vector of the joint. This matrix represents
the subspace of free motion in the full space. It is orthogonal to the constrained subspace.
"""
motion_subspace(j::Joint{ND,JT}) where {ND,JT} = motion_subspace(JT,j.params,Val(ND))

constrained_dimension(j::Joint) = length(j.cdofs)
exogenous_dimension(j::Joint) = length(j.edofs)
unconstrained_dimension(j::Joint) = length(j.udofs)
position_and_vel_dimension(j::Joint) = position_dimension(j) + exogenous_dimension(j) + unconstrained_dimension(j)

function check_q_dimension(q,C::Type{T}) where T<:AbstractJointType
  @assert length(q) == position_dimension(C) "Incorrect length of q"
end

"""
    joint_transform(q::AbstractVector,joint::Joint)

Use the joint dof entries in `q` to create a joint transform operator.
The number of entries in `q` must be apprpriate for the type of joint `joint`.
"""
function joint_transform(q::AbstractVector,joint::Joint{ND,JT}) where {ND,JT<:AbstractJointType}
  @unpack params = joint
  check_q_dimension(q,JT)
  joint_transform(q,JT,params,Val(ND))
end

"""
    parent_to_child_transform(q::AbstractVector,joint::Joint)

Use the joint dof entries in `q` to create a transform operator from the
parent of joint `joint` to the child of joint `joint`.
The number of entries in `q` must be apprpriate for the type of joint `joint`.
"""
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


### Joint evolution equations ###

"""
    zero_joint(joint::Joint[;dimfcn=position_and_vel_dimension])

Create a vector of zeros for different aspects of the joint state, based
on the argument `dimfcn`. By default, it uses `position_and_vel_dimension` and creates a zero vector sized
according to the the position of the joint and the parts of the
joint velocity that must be advanced (from acceleration). Alternatively, one
can use `position_dimension`, `constrained_dimension`, `unconstrained_dimension`,
or `exogenous_dimension`.
"""
function zero_joint(joint::Joint{ND,JT};dimfcn=position_and_vel_dimension) where {ND,JT}
    return zeros(Float64,dimfcn(joint))
end

"""
    init_joint(joint::Joint[;tinit=0.0])

Create an initial state vector for a joint. It initializes the joint's constrained
degrees of freedom with their kinematics, evaluated at time `tinit` (equal to zero, by default).
Other degrees of freedom (exogenous, unconstrained) are initialized to zero, so
they can be set manually.
"""
function init_joint(joint::Joint{ND,JT};tinit = 0.0) where {ND,JT}
    @unpack kins, cdofs = joint
    x = zero_joint(joint)
    q = view(x,1:position_dimension(joint))
    for (i,jdof) in enumerate(cdofs)
      kd = kins[i](tinit)
      q[jdof] = dof_position(kd)
    end
    return x
end

function joint_rhs!(dxdt::AbstractVector,x::AbstractVector,t::Real,a_edof::AbstractVector,a_udof::AbstractVector,joint::Joint{ND,JT}) where {ND,JT}
   @unpack kins, cdofs, edofs, udofs, vbuf = joint

   # Evaluate the velocity of the dofs and place the result in joint.vbuf
   _joint_velocity!(x,t,joint)

   # x holds both q and the parts of qdot that must be advanced
   q = view(x,1:position_dimension(joint))

   qdot = view(dxdt,1:position_dimension(joint))
   aeu = view(dxdt,position_dimension(joint)+1:position_and_vel_dimension(joint))

   # parse the exogenous accelerations into their entries
   for (i,jdof) in enumerate(edofs)
     aeu[i] = a_edof[i]
   end

   # parse the unconstrained accelerations into their entries
   for (i,jdof) in enumerate(udofs)
     aeu[exogenous_dimension(joint)+i] = a_udof[i]
   end

   joint_dqdt!(qdot,q,vbuf,JT)
end

"""
    joint_velocity(x::AbstractVector,t,joint::Joint) -> PluckerMotion

Given a joint `joint`'s full state vector `x`, compute the joint Plucker velocity at time `t`.
"""
function joint_velocity(x::AbstractVector,t::Real,joint::Joint)
    @unpack vbuf = joint

    _joint_velocity!(x,t,joint)

    return PluckerMotion(motion_subspace(joint)*vbuf)

end

function _joint_velocity!(x::AbstractVector,t::Real,joint::Joint)
  @unpack kins, cdofs, edofs, udofs, vbuf = joint

  veu = view(x,position_dimension(joint)+1:position_and_vel_dimension(joint))

  vbuf .= 0.0

  # evaluate the prescribed kinematics at the current time
  for (i,jdof) in enumerate(cdofs)
    kd = kins[i](t)
    vbuf[jdof] = dof_velocity(kd)
  end

  # parse the exogenous velocities into their entries
  for (i,jdof) in enumerate(edofs)
    vbuf[jdof] = veu[i]
  end

  # parse the unconstrained velocities into their entries
  for (i,jdof) in enumerate(udofs)
    vbuf[jdof] = veu[exogenous_dimension(joint)+i]
  end

end


function _joint_dqdt_standard!(dqdt,q,v)
    dqdt .= v
end

function _joint_dqdt_quaternion!(dqdt,q,w)
    dqdt[1] = -q[2]*w[1] - q[3]*w[2] - q[4]*w[3]
    dqdt[2] =  q[1]*w[1] - q[4]*w[2] + q[3]*w[3]
    dqdt[3] =  q[4]*w[1] + q[1]*w[2] - q[2]*w[3]
    dqdt[4] = -q[3]*w[1] + q[2]*w[2] + q[1]*w[3]
    dqdt .*= 0.5
    return dqdt
end


#### Joint types ####

## Default joint velocity calculation

function joint_dqdt!(dqdt,q,v::Vector,::Type{T}) where T<:AbstractJointType
  _joint_dqdt_standard!(dqdt,q,v)
  nothing
end

function motion_subspace(JT::Type{<:AbstractJointType},p::Dict,::Val{2})
  S3 = motion_subspace(JT,p,Val(3))
  m, n = size(S3)
  SMatrix{3,n}(S3[3:5,:])
end

## Revolute joint ##

abstract type RevoluteJoint <: AbstractJointType end

number_of_dofs(::Type{RevoluteJoint}) = 1
position_dimension(::Type{RevoluteJoint}) = 1

function joint_transform(q::AbstractVector,::Type{RevoluteJoint},p::Dict,::Val{ND}) where {ND}
  R = rotation_about_z(-q[1])
  x = SVector{3}([0.0,0.0,0.0])
  return MotionTransform{ND}(x,R)
end

motion_subspace(::Type{RevoluteJoint},p::Dict,::Val{3}) = SMatrix{6,1,Float64}([0 0 1 0 0 0]')

## Prismatic joint ##

abstract type PrismaticJoint <: AbstractJointType end

number_of_dofs(::Type{PrismaticJoint}) = 1
position_dimension(::Type{PrismaticJoint}) = 1

function joint_transform(q::AbstractVector,::Type{PrismaticJoint},p::Dict,::Val{3})
  R = rotation_identity()
  x = SVector{3}([0.0,0.0,q[1]])
  return MotionTransform{3}(x,R)
end

motion_subspace(::Type{PrismaticJoint},p::Dict,::Val{3}) = SMatrix{6,1,Float64}([0 0 0 0 0 1]')


## Helical joint ##

abstract type HelicalJoint <: AbstractJointType end

number_of_dofs(::Type{HelicalJoint}) = 1
position_dimension(::Type{HelicalJoint}) = 1

function joint_transform(q::AbstractVector,C::Type{HelicalJoint},p::Dict,::Val{3})
  R = rotation_about_z(-q[1])
  x = SVector{3}([0.0,0.0,p["pitch"]*q[1]])
  return MotionTransform{3}(x,R)
end

motion_subspace(::Type{HelicalJoint},p::Dict,::Val{3}) = SMatrix{6,1,Float64}([0 0 1 0 0 p["pitch"]]')

## Cylindrical joint ##

abstract type CylindricalJoint <: AbstractJointType end

number_of_dofs(::Type{CylindricalJoint}) = 2
position_dimension(::Type{CylindricalJoint}) = 2

function joint_transform(q::AbstractVector,::Type{CylindricalJoint},p::Dict,::Val{3})
  R = rotation_about_z(-q[1])
  x = SVector{3}([0.0,0.0,q[2]])
  return MotionTransform{3}(x,R)
end

motion_subspace(::Type{CylindricalJoint},p::Dict,::Val{3}) = SMatrix{6,2,Float64}([0 0 1 0 0 0; 0 0 0 0 0 1]')

## Spherical joint ##

abstract type SphericalJoint <: AbstractJointType end

number_of_dofs(::Type{SphericalJoint}) = 3
position_dimension(::Type{SphericalJoint}) = 4

function joint_transform(q::AbstractVector,::Type{SphericalJoint},p::Dict,::Val{3})
  R = rotation_from_quaternion(quaternion(q))
  x = SVector{3}([0.0,0.0,0.0])
  return MotionTransform{3}(x,R)
end

motion_subspace(::Type{SphericalJoint},p::Dict,::Val{3}) = SMatrix{6,3,Float64}([1 0 0 0 0 0;
                                                                                 0 1 0 0 0 0;
                                                                                 0 0 1 0 0 0]')

function joint_dqdt!(dqdt,q,v,::Type{SphericalJoint})
  _joint_dqdt_quaternion!(dqdt,q,v)
  nothing
end

## Free joint (3d) ##

abstract type FreeJoint <: AbstractJointType end

number_of_dofs(::Type{FreeJoint}) = 6
position_dimension(::Type{FreeJoint}) = 7

function joint_transform(q::AbstractVector,::Type{FreeJoint},p::Dict,::Val{3})
  R = rotation_from_quaternion(quaternion(q[1:4]))
  x = SVector{3}(q[5:7])
  return MotionTransform{3}(R'*x,R)
end

motion_subspace(::Type{FreeJoint},p::Dict,::Val{3}) = SMatrix{6,6,Float64}(I)


function joint_dqdt!(dqdt,q,v,::Type{FreeJoint})
  qrdot, qtdot = view(dqdt,1:4), view(dqdt,5:7)
  qr, qt = view(q,1:4), view(q,5:7)
  vr, vt = view(v,1:3), view(v,4:6)
  _joint_dqdt_quaternion!(qrdot,qr,vr)
  _joint_dqdt_standard!(qtdot,qt,vt)
  nothing
end

## Free joint (2d) ##

abstract type FreeJoint2d <: AbstractJointType end

number_of_dofs(::Type{FreeJoint2d}) = 3
position_dimension(::Type{FreeJoint2d}) = 3

function joint_transform(q::AbstractVector,::Type{FreeJoint2d},p::Dict,::Val{2})
  R = rotation_about_z(-q[1])
  x = SVector{3}(q[2],q[3],0.0)
  return MotionTransform{2}(x,R)
end

motion_subspace(::Type{FreeJoint2d},p::Dict,::Val{2}) = SMatrix{3,3,Float64}(I)

## Fixed joint ##

abstract type FixedJoint <: AbstractJointType end

"""
    Joint(X::MotionTransform,bid::Int)

Construct a joint that simply places the body with ID `bid` rigidly in the configuration
given by `X`.
"""
Joint(Xp_to_j::MotionTransform{ND},child_id::Int) where {ND} = Joint(FixedJoint,0,Xp_to_j,child_id,MotionTransform{ND}())

"""
    Joint(X::MotionTransform)

For a problem with a single body, construct a joint that simply places the body rigidly in the configuration
given by `X`.
"""
Joint(Xp_to_j::MotionTransform{ND}) where {ND} = Joint(FixedJoint,0,Xp_to_j,1,MotionTransform{ND}())


function Joint(::Type{FixedJoint},parent_id::Int,Xp_to_j::MotionTransform{2},child_id::Int,Xch_to_j::MotionTransform{2})
    dofs = [ConstantVelocityDOF(0.0) for i = 1:3]
    Joint(FreeJoint2d,parent_id,Xp_to_j,child_id,Xch_to_j,dofs)
end

function Joint(::Type{FixedJoint},parent_id::Int,Xp_to_j::MotionTransform{3},child_id::Int,Xch_to_j::MotionTransform{3})
    dofs = [ConstantVelocityDOF(0.0) for i = 1:6]
    Joint(FreeJoint,parent_id,Xp_to_j,child_id,Xch_to_j,dofs)
end
