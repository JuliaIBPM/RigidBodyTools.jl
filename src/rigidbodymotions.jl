# Linked systems #

"""
    RigidBodyMotion

Type containing the connectivities and motions of rigid bodies, linked to the
inertial system and possibly to each other, via joints. The basic constructor is
`RigidBodyMotion(joints::Vector{Joint},nbody::Int)`, in which `joints` contains
a vector of joints of [`Joint`](@ref) type, each specifying the connection between
a parent and a child body. (The parent may be the inertial coordinate system.)
"""
struct RigidBodyMotion{ND}
    nls :: Int
    nbody :: Int

    "List of linked bodies, starting with the one connected to inertial system"
    lslists :: Vector{Vector{Int}}
    joints :: Vector{Joint}
    parent_body :: Vector{Int}
    parent_joint :: Vector{Int}
    child_bodies :: Vector{Vector{Int}}
    child_joints :: Vector{Vector{Int}}

    state_indices :: Vector{Int}
    vel_indices :: Vector{Int}

end


function RigidBodyMotion(joints::Vector{<:Joint},nbody::Int)
    parent_body = zeros(Int,nbody)
    parent_joint = zeros(Int,nbody)
    child_bodies = [Int[] for i in 1:nbody]
    child_joints = [Int[] for i in 1:nbody]

    lscnt = 0
    lsroots = Int[]
    for bodyid in 1:nbody
        parent_joint[bodyid] = pjid = _parent_joint_of_body(bodyid,joints)
        parent_body[bodyid] = pbid = joints[pjid].parent_id
        if pbid == 0
            lscnt += 1
            push!(lsroots,bodyid)
        end
        child_joints[bodyid] = cjids = _child_joints_of_body(bodyid,joints)
        child_bodies[bodyid] = [joints[cjid].child_id for cjid in cjids]
    end
    @assert lscnt > 0 "There is no body connected to the inertial system"

    lslists = [Int[] for lsid in 1:lscnt]
    for (lsid,lsroot) in enumerate(lsroots)
        lslist = lslists[lsid]
        _list_of_linked_bodies!(lslist,lsroot,child_bodies)
    end
    total_members = mapreduce(x -> length(x),+,lslists)
    @assert total_members == nbody "Not all bodies have been added as members of linked lists"

    state_indices = mapreduce(jid -> _getrange(joints,state_and_vel_dimension,jid)[1:state_dimension(joints[jid])],
                              vcat,eachindex(joints))
    vel_indices = mapreduce(jid -> _getrange(joints,state_and_vel_dimension,jid)[state_dimension(joints[jid])+1:state_and_vel_dimension(joints[jid])],
                              vcat,eachindex(joints))

    RigidBodyMotion{physical_dimension(joints)}(lscnt,nbody,lslists,joints,parent_body,parent_joint,
                                              child_bodies,child_joints,state_indices,vel_indices)
end

function Base.show(io::IO, ls::RigidBodyMotion)
    println(io, "$(ls.nls) linked system(s) of bodies")
end

function _list_of_linked_bodies!(lslist,bodyid,child_bodies)
    push!(lslist,bodyid)
    for childid in child_bodies[bodyid]
        _list_of_linked_bodies!(lslist,childid,child_bodies)
    end
    nothing
end


# If this is empty, then there should be an error, since all bodies
# need a parent joint. If it has more than one result, also an error,
# since a body can only have one parent joint
function _parent_joint_of_body(bodyid,joints::Vector{<:Joint})
    idx = findall(x -> x.child_id == bodyid, joints)
    @assert length(idx) > 0 "Body "*string(bodyid)*" has no parent joint"
    @assert length(idx) == 1 "Body "*string(bodyid)*" has too many parent joints"
    return first(idx)
end

function _child_joints_of_body(bodyid,joints::Vector{<:Joint})
    idx = findall(x -> x.parent_id == bodyid, joints)
    return idx
end

number_of_linked_systems(ls::RigidBodyMotion) = ls.nls

"""
    first_body(lsid::Int,ls::RigidBodyMotion)

Return the index of the first body in linked system `lsid` in
the overall set of linked systems `ls`.  This body's parent is the
inertial coordinate system.
"""
function first_body(lsid::Int,ls::RigidBodyMotion)
  @unpack lslists, parent_joint = ls
  first(lslists[lsid])
end

"""
    first_joint(lsid::Int,ls::RigidBodyMotion)

Return the index of the first joint in linked system `lsid` in
the overall set of linked systems `ls`. This joint's parent is the
inertial coordinate system.
"""
function first_joint(lsid::Int,ls::RigidBodyMotion)
  @unpack lslists, parent_joint = ls
  parent_joint[first_body(lsid,ls)]
end


for f in [:number_of_dofs,:state_dimension,:state_and_vel_dimension,:constrained_dimension,
           :exogenous_dimension,:unconstrained_dimension]
   @eval $f(ls::RigidBodyMotion) = mapreduce(x -> $f(x),+,ls.joints)
end


"""
    getrange(ls::RigidBodyMotion,dimfcn::Function,jid::Int) -> Range

Return the subrange of indices in the global state vector
for the state corresponding to joint `jid` in linked system `ls`.
"""
function getrange(ls::RigidBodyMotion,dimfcn::Function,jid::Int)
    @unpack joints = ls
    _getrange(joints,dimfcn,jid)
end

function _getrange(joints::Vector{<:Joint},dimfcn::Function,jid::Int)
    0 < jid <= length(joints) || error("Unavailable joint")
    first = 1
    j = 1
    while j < jid
        first += dimfcn(joints[j])
        j += 1
    end
    last = first+dimfcn(joints[jid])-1
    return first:last
end


"""
    view(q::AbstractVector,ls::RigidBodyMotion,jid::Int[;dimfcn=state_dimension]) -> SubArray

Provide a view of the range of values in vector `q` corresponding to the state
of the joint with index `jid` in a RigidBodyMotion `ls`. The optional argument `dimfcn`
can be set to `state_dimension`, `constrained_dimension`, `unconstrained_dimension`,
or `exogenous_dimension`.
"""
function Base.view(q::AbstractVector,ls::RigidBodyMotion,jid::Int;dimfcn::Function=state_dimension)
    length(q) == dimfcn(ls) || error("Inconsistent size of data for viewing")
    return view(q,getrange(ls,dimfcn,jid))
end

"""
    linked_system_transform(q::AbstractVector,ls::RigidBodyMotion) -> MotionTransformList

Parse the overall state vector `q` into the individual joints and construct
the inertial system-to-body transforms for every body. Return these
transforms in a `MotionTransformList`.
"""
function linked_system_transform(q::AbstractVector,ls::RigidBodyMotion{ND}) where {ND}
    X = MotionTransform{ND}()
    ml = [deepcopy(X) for jb in 1:ls.nbody]
    for lsid in 1:number_of_linked_systems(ls)
      jid = first_joint(lsid,ls)
      _joint_descendants_transform!(ml,X,jid,q,ls)
    end
    return MotionTransformList(ml)
end

function _joint_descendants_transform!(ml,Xp::MotionTransform,jid::Int,q::AbstractVector,ls::RigidBodyMotion)
    @unpack joints, child_joints = ls
    joint = joints[jid]
    Xch = parent_to_child_transform(view(q,ls,jid),joint)*Xp
    bid = joint.child_id
    ml[bid] = deepcopy(Xch)
    for jcid in child_joints[bid]
        _joint_descendants_transform!(ml,Xch,jcid,q,ls)
    end
    nothing
end



"""
    body_velocities(x::AbstractVector,t::Real,ls::RigidBodyMotion) ->

Compute the velocities of bodies for the rigid body system `ls`, expressing each in its own coordinate system.
To carry this out, the function evaluates velocities of dofs with prescribed kinematics at time `t` and
and obtains the remaining free dofs (exogenous and unconstrained) from the state/velocity vector `x`.
"""
function body_velocities(x::AbstractVector,t::Real,ls::RigidBodyMotion{ND}) where {ND}
    v0 = PluckerMotion{2}()
    vl = [v0 for jb in 1:ls.nbody]
    for lsid in 1:number_of_linked_systems(ls)
      jid = first_joint(lsid,ls)
      _child_velocity_from_parent!(vl,v0,jid,x,t,ls)
    end
    return PluckerMotionList(vl)
end

function _child_velocity_from_parent!(vl,vp::PluckerMotion,jid::Int,x::AbstractVector,t::Real,ls::RigidBodyMotion)
    @unpack joints, child_joints = ls

    joint = joints[jid]
    xJ = view(x,ls,jid;dimfcn=state_and_vel_dimension)
    q = statevector(x,ls)
    vJ = joint_velocity(xJ,t,joint)
    qJ = view(q,ls,jid)

    Xp_to_ch = parent_to_child_transform(qJ,joint)
    xJ_to_ch = inv(joint.Xch_to_j)

    vch = _child_velocity_from_parent(Xp_to_ch,vp,xJ_to_ch,vJ)
    bid = joint.child_id
    vl[bid] = deepcopy(vch)

    for jcid in child_joints[bid]
        _child_velocity_from_parent!(vl,vch,jcid,x,t,ls)
    end
    nothing

end


function _child_velocity_from_parent(Xp_to_ch::MotionTransform,vp::PluckerMotion,XJ_to_ch::MotionTransform,vJ::PluckerMotion)
    vch = Xp_to_ch*vp + XJ_to_ch*vJ
end


"""
    zero_joint(ls::RigidBodyMotion[;dimfcn=state_and_vel_dimension])

Create a vector of zeros for the some aspect of the state of the linked system(s) `ls`, based
on the argument `dimfcn`. By default, it uses `state_and_vel_dimension` and creates a zero vector sized
according to the the state of the joint and the parts of the
joint velocity that must be advanced (from acceleration). Alternatively, one
can use `state_dimension`, `constrained_dimension`, `unconstrained_dimension`,
or `exogenous_dimension`.
"""
function zero_joint(ls::RigidBodyMotion;dimfcn=state_and_vel_dimension)
    mapreduce(joint -> zero_joint(joint;dimfcn=dimfcn),vcat,ls.joints)
end

"""
    init_joint(ls::RigidBodyMotion[;tinit = 0.0])

Initialize the global linked system state and velocity vector, using
the prescribed motions for constrained degrees of freedom to initialize
those components (evaluated at `tinit`, which by default is 0).
"""
function init_joint(ls::RigidBodyMotion;kwargs...)
    mapreduce(joint -> init_joint(joint;kwargs...),vcat,ls.joints)
end

"""
    statevector(x::AbstractVector,ls::RigidBodyMotion)

Returns a view of the global state/velocity vector for a linked system containing only the state.
"""
function statevector(x::AbstractVector,ls::RigidBodyMotion)
   @unpack state_indices = ls
   return view(x,state_indices)
end

"""
    velvector(x::AbstractVector,ls::RigidBodyMotion)

Returns a view of the global state/velocity vector for a linked system containing only the velocity.
"""
function velvector(x::AbstractVector,ls::RigidBodyMotion)
   @unpack vel_indices = ls
   return view(x,vel_indices)
end


"""
    joint_rhs!(dxdt::AbstractVector,x::AbstractVector,t::Real,a_edof,a_udof,ls::RigidBodyMotion)

Sets the right-hand side vector `dxdt` (mutating) for linked system `ls`, using the current state/velocity vector `x`,
the current time `t`, exogenous accelerations `a_edof` and unconstrained accelerations `a_udof`.
"""
function joint_rhs!(dxdt::AbstractVector,x::AbstractVector,t::Real,a_edof::AbstractVector,a_udof::AbstractVector,ls::RigidBodyMotion)
    @unpack joints = ls
    for (jid,joint) in enumerate(joints)
        dxdt_j = view(dxdt,ls,jid;dimfcn=state_and_vel_dimension)
        x_j = view(x,ls,jid;dimfcn=state_and_vel_dimension)
        a_edof_j = view(a_edof,ls,jid;dimfcn=exogenous_dimension)
        a_udof_j = view(a_udof,ls,jid;dimfcn=unconstrained_dimension)
        joint_rhs!(dxdt_j,x_j,t,a_edof_j,a_udof_j,joint)
    end
    nothing
end


function velocity_in_body_coordinates_2d(x̃,ỹ,vb::PluckerMotion)
    Xb_to_p = MotionTransform(x̃,ỹ,0.0)
    Xb_to_p*vb
end

function velocity_in_inertial_coordinates_2d(x̃,ỹ,vb::PluckerMotion,Xb_to_0::MotionTransform)
    rotation_transform(Xb_to_0)*velocity_in_body_coordinates_2d(x̃,ỹ,vb)
end

function velocity_in_body_coordinates_2d!(u::AbstractVector,v::AbstractVector,b::Body,vb::PluckerMotion)
    for i in 1:numpts(b)
        vp = velocity_in_body_coordinates_2d(b.x̃[i],b.ỹ[i],vb)
        u[i] = vp[2]
        v[i] = vp[3]
    end
end

function velocity_in_inertial_coordinates_2d!(u::AbstractVector,v::AbstractVector,b::Body,vb::PluckerMotion,Xb_to_0::MotionTransform)
    Rb_to_0 = rotation_transform(Xb_to_0)
    for i in 1:numpts(b)
        vp = Rb_to_0*velocity_in_body_coordinates_2d(b.x̃[i],b.ỹ[i],vb)
        u[i] = vp[2]
        v[i] = vp[3]
    end
end



#=
"""
    RigidBodyMotion(Up::Union{ComplexF64,Tuple},Ω[;pivot::Union{ComplexF64,Tuple}=(0.0,0.0)])

Create an instance of constant rigid-body motion with translational velocity Up
and angular velocity Ω with respect to a specified point P.
This point is established by its initial position `pivot`. By default, this
initial position is (0,0). `Up` and `pivot` can be specified by either Tuple or by complex value.
"""
RigidBodyMotion(Up::Union{Number,Tuple}, Ω::Number;kwargs...) = (kin = Constant(Up, Ω; kwargs...); RigidBodyMotion(kin(0), kin))

"""
    RigidBodyMotion(kin::AbstractKinematics)

Create an instance of rigid-body motion with kinematics `kin`.
"""
RigidBodyMotion(kin::AbstractKinematics) = RigidBodyMotion(kin(0), kin)
(m::RigidBodyMotion)(t) = m.kin(t)


function (m::RigidBodyMotion)(t,x̃::Union{Number,Tuple})
  # This expects coordinates in a comoving coordinate system
  # based at c.
  z̃ = complex(x̃...)
  #m.c, m.ċ, m.c̈, m.α, m.α̇, m.α̈ = m.kin(t)
  k = m.kin(t)
  c = complex_translational_position(k)
  ċ = complex_translational_velocity(k)
  c̈ = complex_translational_acceleration(k)
  α = angular_position(k)
  α̇ = angular_velocity(k)
  α̈ = angular_acceleration(k)
  z = exp(im*α)*z̃
  return KinematicData(t, c + z, ċ + im*α̇*z, c̈ + (im*α̈-α̇^2)*z, α, α̇, α̈)
end

"""
    motion_velocity(b::Body,m::RigidBodyMotion,t::Real)

Return the velocity components (as a vector) of a `RigidBodyMotion`
at the given time `t`.
"""
function motion_velocity(b::Body,motion::RigidBodyMotion,t::Real)
   k = motion(t)
  #_,ċ,_,_,α̇,_ = motion(t)
  return [translational_velocity(k)...,angular_velocity(k)]
end

"""
    motion_state(b::Body,m::RigidBodyMotion)

Return the current state vector of body `b` associated with
rigid body motion `m`. It returns the current coordinates
of the body centroid and the angle of the body.
"""
function motion_state(b::Body,m::RigidBodyMotion)
    return [b.cent...,b.α]
end

"""
    update_body!(b::Body,x::AbstractVector,m::RigidBodyMotion)

Update body `b` with the rigid-body motion state vector `x`. The information
in `m` is used for checking lengths only.
"""
function update_body!(b::Body,x::AbstractVector,m::RigidBodyMotion)
    length(x) == length(motion_state(b,m)) || error("wrong length for motion state vector")
    T = RigidTransform(x)
    T(b)
    return b
end


"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 body::Body,motion::AbstractMotion,t::Real[;inertial=true])

Assign the components of body velocity `u` and `v` (in inertial coordinate system)
at surface positions described by inertial coordinates in body `body` at time `t`,
based on supplied motions in the `motion` for the body. If, instead,
`inertial=false`, then it is assumed that `u`, `v` and the surface positions in `body` are all
expressed in comoving coordinates (but velocities are relative to inertial
  frame).
"""
surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,m::RigidBodyMotion,t::Real;kwargs...) =
                 surface_velocity!(u,v,b.x,b.y,b.cent...,b.α,m,t;kwargs...)

"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     x::AbstractVector{Float64},y::AbstractVector{Float64},
                     xc::Real,yc::Real,α::Real,
                     motion::RigidBodyMotion,t::Real[;inertial=true])

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at surface positions described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body. The current position and orientation
of the rigid-body coordinate system are supplied as `xc`, `yc`, and `α`. If, instead,
`inertial=false`, then it is assumed that `u`, `v`, `x`, `y`, `xc`, `yc` are all
expressed in comoving coordinates (but velocities are relative to inertial
  frame).
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                          x::AbstractVector{Float64},y::AbstractVector{Float64},
                          xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real;inertial=true)

  length(u) == length(v) == length(x) == length(y) || error("Inconsistent lengths of vectors")

  k = m(t)
  ċ = complex_translational_velocity(k;inertial=inertial)
  α̇ = angular_velocity(k)
  uc = ċ .+ im*α̇*((x .- xc) .+ im*(y .- yc))
  u .= real.(uc)
  v .= imag.(uc)

  return u, v
end

"""
    surface_velocity(x::AbstractVector{Float64},y::AbstractVector{Float64},
                    xc::Real,yc::Real,α::Real,motion::RigidBodyMotion,t::Real[;inertial=true])

Return the components of rigid body velocities (in inertial components) at surface positions
described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body. If, instead,
`inertial=false`, then it is assumed that `u`, `v`, `x`, `y`, `xc`, `yc` are all
expressed in comoving coordinates (but velocities are relative to inertial
  frame).
"""
surface_velocity(x::AbstractVector{Float64},y::AbstractVector{Float64},
                xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real;kwargs...) =
                surface_velocity!(similar(x),similar(y),x,y,xc,yc,α,m,t;kwargs...)


function show(io::IO, m::RigidBodyMotion)
    println(io, "Rigid Body Motion:")
    print(io, "  $(m.kin)")
end
=#
