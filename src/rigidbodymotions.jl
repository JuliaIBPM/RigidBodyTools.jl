# Linked systems #

function _default_exogenous_function!(a_edof,u,p,t) end

function _default_unconstrained_function!(a_udof,u,p,t) end

"""
    RigidBodyMotion

Type containing the connectivities and motions of rigid bodies, linked to the
inertial system and possibly to each other, via joints. The basic constructor is
`RigidBodyMotion(joints::Vector{Joint},nbody::Int)`, in which `joints` contains
a vector of joints of [`Joint`](@ref) type, each specifying the connection between
a parent and a child body. (The parent may be the inertial coordinate system.)
"""
struct RigidBodyMotion{ND} <: AbstractMotion

    "Number of distinct linked systems"
    nls :: Int

    "Number of bodies"
    nbody :: Int

    "List of linked bodies, starting with the one connected to inertial system"
    lslists :: Vector{Vector{Int}}

    "Vector of joint specifications and motions"
    joints :: Vector{Joint}

    "Vector of prescribed deformations of bodies"
    deformations :: Vector{AbstractDeformationMotion}

    "Vector specifying the parent body of a body (0 = inertial system)"
    parent_body :: Vector{Int}

    "Vector specifying the parent joint of a body"
    parent_joint :: Vector{Int}

    "Vector of vectors, specifying the list of child bodies of a body (empty if no children)"
    child_bodies :: Vector{Vector{Int}}

    "Vector of vectors, specifying the list of child joints of a body (empty if no children)"
    child_joints :: Vector{Vector{Int}}

    "Set of indices providing the joint position coordinates in the global state vector"
    position_indices :: Vector{Int}

    "Set of indices providing the joint velocity coordinates in the global state vector"
    vel_indices :: Vector{Int}

    "Sets of indices providing the body point coordinates in the global state vector"
    deformation_indices :: Vector{Vector{Int}}

    "Function to specify exogenous behavior, must be mutating and have signature (a,u,p,t)"
    exogenous_function! :: Function

    "Buffers for exogenous and unconstrained behavior"
    a_edof_buffer :: Vector{Float64}

    "Buffers for exogenous and unconstrained behavior"
    a_udof_buffer :: Vector{Float64}

end


function RigidBodyMotion(joints::Vector{<:Joint},bodies::BodyList,deformations::Vector{<:AbstractDeformationMotion};
                              exogenous::Function = _default_exogenous_function!)
    nbody = length(bodies)
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

    position_indices = mapreduce(jid -> _getrange(joints,position_and_vel_dimension,jid)[1:position_dimension(joints[jid])],
                              vcat,eachindex(joints))
    vel_indices = mapreduce(jid -> _getrange(joints,position_and_vel_dimension,jid)[position_dimension(joints[jid])+1:position_and_vel_dimension(joints[jid])],
                              vcat,eachindex(joints))

    index_i = length(position_indices) + length(vel_indices)
    deformation_indices = Vector{Int}[]
    for (b,m) in zip(bodies,deformations)
        def_length = length(motion_state(b,m))
        push!(deformation_indices,collect(index_i+1:index_i+def_length))
        index_i += def_length
    end

    a_edof_buffer = _zero_joint(joints,exogenous_dimension)
    a_udof_buffer = _zero_joint(joints,unconstrained_dimension)


    RigidBodyMotion{physical_dimension(joints)}(lscnt,nbody,lslists,joints,deformations,parent_body,parent_joint,
                                              child_bodies,child_joints,position_indices,vel_indices,deformation_indices,
                                              exogenous,a_edof_buffer,a_udof_buffer)
end

RigidBodyMotion(joints::Vector{<:Joint},bodies::BodyList;kwargs...) = RigidBodyMotion(joints,bodies,[NullDeformationMotion() for bi in 1:length(bodies)];kwargs...)

RigidBodyMotion(joint::Joint,body::Body,def::AbstractDeformationMotion;kwargs...) = RigidBodyMotion([joint],BodyList([body]),[def];kwargs...)
RigidBodyMotion(joint::Joint,body::Body;kwargs...) = RigidBodyMotion([joint],BodyList([body]);kwargs...)


function Base.show(io::IO, ls::RigidBodyMotion)
    println(io, "$(ls.nls) linked system(s) of bodies")
    println(io, "   $(ls.nbody) bodies")
    println(io, "   $(length(ls.joints)) joints")
end

"""
    ismoving(m::RigidBodyMotion)

Checks if any joint in `m` is in motion
"""
function ismoving(m::RigidBodyMotion)
    @unpack joints, deformations = m
    any(map(joint -> ismoving(joint),joints)) || any(map(def -> ismoving(def),deformations))
end

"""
    is_system_in_relative_motion(lsid::Int,ls::RigidBodyMotion) -> Bool

Returns `true` if the linked system with ID `lsid` is in relative motion,
i.e., if one or more of the joints (not connected to the inertial system)
returns `true` for the function `ismoving(joint)`.
"""
function is_system_in_relative_motion(lsid::Int,m::RigidBodyMotion)
    @unpack lslists = m
    bodyid = first_body(lsid,m)
    ismove = false
    return _ismoving_tree(ismove,bodyid,m)
end

function _ismoving_tree(ismove::Bool,bodyid::Int,m::RigidBodyMotion)
    @unpack joints, deformations = m
    ismove |= ismoving(deformations[bodyid])
    for jid in child_joints_of_body(bodyid,m)
        ismove |= ismoving(joints[jid])
        bodyid = child_body_of_joint(jid,m)
        ismove = _ismoving_tree(ismove,bodyid,m)
    end
    return ismove
end

function _list_of_linked_bodies!(lslist,bodyid,child_bodies)
    push!(lslist,bodyid)
    for childid in child_bodies[bodyid]
        _list_of_linked_bodies!(lslist,childid,child_bodies)
    end
    nothing
end

function _check_for_only_one_body(ls::RigidBodyMotion)
    @assert ls.nbody == 1 "Linked system must only contain one body"
end

"""
    parent_joint_of_body(bodyid,ls::RigidBodyMotion)

Return the index of the parent joint of body `bodyid` in system `ls`.
"""
parent_joint_of_body(bodyid::Int,ls::RigidBodyMotion) = _parent_joint_of_body(bodyid,ls.joints)

"""
    child_joints_of_body(bodyid,ls::RigidBodyMotion)

Return the indices (if any) of the child joints of body `bodyid` in system `ls`.
This function also accepts `bodyid = 0`, and returns the indices of joints
that connect bodies to the inertial coordinate system.
"""
child_joints_of_body(bodyid::Int,ls::RigidBodyMotion) = _child_joints_of_body(bodyid,ls.joints)

"""
    parent_body_of_joint(jid,ls::RigidBodyMotion)

Return the index of the parent body of joint `jid` in system `ls`.
"""
parent_body_of_joint(jid::Int,ls::RigidBodyMotion) = ls.joints[jid].parent_id

"""
    child_body_of_joint(jid,ls::RigidBodyMotion)

Return the index of the child body of joint `jid` in system `ls`.
"""
child_body_of_joint(jid::Int,ls::RigidBodyMotion) = ls.joints[jid].child_id



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




for f in [:number_of_dofs,:position_dimension,:position_and_vel_dimension,:constrained_dimension,
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
    view(q::AbstractVector,ls::RigidBodyMotion,jid::Int[;dimfcn=position_dimension]) -> SubArray

Provide a view of the range of values in vector `q` corresponding to the position
of the joint with index `jid` in a RigidBodyMotion `ls`. The optional argument `dimfcn`
can be set to `position_dimension`, `constrained_dimension`, `unconstrained_dimension`,
or `exogenous_dimension`.
"""
function Base.view(q::AbstractVector,ls::RigidBodyMotion,jid::Int;dimfcn::Function=position_dimension)
    return view(q,getrange(ls,dimfcn,jid))
end

"""
    body_transforms(x::AbstractVector,ls::RigidBodyMotion) -> MotionTransformList

Parse the overall state vector `x` into the individual joints and construct
the inertial system-to-body transforms for every body. Return these
transforms in a `MotionTransformList`.
"""
function body_transforms(x::AbstractVector,ls::RigidBodyMotion{ND}) where {ND}
    q = position_vector(x,ls)
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
    body_velocities(x::AbstractVector,t::Real,ls::RigidBodyMotion) -> PluckerMotionList

Compute the velocities of bodies for the rigid body system `ls`, expressing each in its own coordinate system.
To carry this out, the function evaluates velocities of dofs with prescribed kinematics at time `t` and
and obtains the remaining free dofs (exogenous and unconstrained) from the state vector `x`.
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

function _child_velocity_from_parent!(vl,vp::AbstractPluckerMotionVector,jid::Int,x::AbstractVector,t::Real,ls::RigidBodyMotion)
    @unpack joints, child_joints = ls

    joint = joints[jid]
    xJ = view(x,ls,jid;dimfcn=position_and_vel_dimension)
    q = position_vector(x,ls)
    vJ = joint_velocity(xJ,t,joint)
    qJ = view(q,ls,jid)

    Xp_to_ch = parent_to_child_transform(qJ,joint)
    #xJ_to_ch = Xp_to_ch*inv(joint.Xp_to_j)
    xJ_to_ch = inv(joint.Xch_to_j)


    vch = _child_velocity_from_parent(Xp_to_ch,vp,xJ_to_ch,vJ)
    bid = joint.child_id
    vl[bid] = deepcopy(vch)

    for jcid in child_joints[bid]
        _child_velocity_from_parent!(vl,vch,jcid,x,t,ls)
    end
    nothing

end


function _child_velocity_from_parent(Xp_to_ch::MotionTransform,vp::AbstractPluckerMotionVector,XJ_to_ch::MotionTransform,vJ::AbstractPluckerMotionVector)
    vch = Xp_to_ch*vp + XJ_to_ch*vJ
end


"""
    zero_joint(ls::RigidBodyMotion[;dimfcn=position_and_vel_dimension])

Create a vector of zeros for some aspect of the state of the linked system(s) `ls`, based
on the argument `dimfcn`. By default, it uses `position_and_vel_dimension` and creates a zero vector sized
according to the state of the joint. Alternatively, one
can use `position_dimension`, `constrained_dimension`, `unconstrained_dimension`,
or `exogenous_dimension`.
"""
zero_joint(ls::RigidBodyMotion;dimfcn=position_and_vel_dimension) = _zero_joint(ls.joints,dimfcn)


_zero_joint(joints::Vector{<:Joint},dimfcn) = mapreduce(joint -> zero_joint(joint;dimfcn=dimfcn),vcat,joints)

"""
    zero_motion_state(bl::BodyList,ls::RigidBodyMotion)

Create a vector of zeros for the state of the linked system(s) `ls`.
"""
function zero_motion_state(bl::BodyList,ls::RigidBodyMotion)
    @unpack deformations = ls
    x = zero_joint(ls)
    for (b,m) in zip(bl,deformations)
      xi = zero_motion_state(b,m)
      append!(x,xi)
    end
    return x
end

function zero_motion_state(b::Body,ls::RigidBodyMotion)
    _check_for_only_one_body(ls)
    zero_motion_state(BodyList([b]),ls)
end

"""
    init_motion_state(bl::BodyList,ls::RigidBodyMotion[;tinit = 0.0])

Initialize the global linked system state vector, using
the prescribed motions for constrained degrees of freedom to initialize
the position components (evaluated at `tinit`, which by default is 0).
The other degrees of freedom are initialized to zero, and can be
set manually.
"""
function init_motion_state(bl::BodyList,ls::RigidBodyMotion;kwargs...)
    @unpack deformations = ls
    x = mapreduce(joint -> init_joint(joint;kwargs...),vcat,ls.joints)
    for (b,m) in zip(bl,deformations)
      xi = motion_state(b,m)
      append!(x,xi)
    end
    return x
end

function init_motion_state(b::Body,ls::RigidBodyMotion;kwargs...)
    _check_for_only_one_body(ls)
    init_motion_state(BodyList([b]),ls;kwargs...)
end

"""
    position_vector(x::AbstractVector,ls::RigidBodyMotion)

Returns a view of the global state vector for a linked system containing only the position.
"""
function position_vector(x::AbstractVector,ls::RigidBodyMotion)
   @unpack position_indices = ls
   return view(x,position_indices)
end

"""
    position_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)

Returns a view of the global state vector for a linked system containing only the position
of joint `jid`.
"""
function position_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)
   @unpack position_indices = ls
   q = view(x,position_indices)
   return view(q,ls,jid)
end

"""
    exogenous_position_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)

Returns a view of the exogenous dof position(s), if any, of joint `jid` in
the global state vector `x`.
"""
function exogenous_position_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)
  @unpack joints = ls
  qj = position_vector(x,ls,jid)
  return view(qj,joints[jid].edofs)
end

"""
    unconstrained_position_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)

Returns a view of the unconstrained dof position(s), if any, of joint `jid` in
the global state vector `x`.
"""
function unconstrained_position_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)
  @unpack joints = ls
  qj = position_vector(x,ls,jid)
  return view(qj,joints[jid].udofs)
end


"""
    velocity_vector(x::AbstractVector,ls::RigidBodyMotion)

Returns a view of the global state vector for a linked system containing only the velocity.
"""
function velocity_vector(x::AbstractVector,ls::RigidBodyMotion)
   @unpack vel_indices = ls
   return view(x,vel_indices)
end

"""
    velocity_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)

Returns a view of the global state vector for a linked system containing only the velocity
of joint `jid`.
"""
function velocity_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)
   @unpack joints = ls
   ir = _getrange(joints,position_and_vel_dimension,jid)[position_dimension(joints[jid])+1:position_and_vel_dimension(joints[jid])]
   return view(x,ir)
end

"""
    exogenous_velocity_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)

Returns a view of the exogenous dof position(s), if any, of joint `jid` in
the global state vector `x`.
"""
function exogenous_velocity_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)
  @unpack joints = ls
  qdotj = velocity_vector(x,ls,jid)
  return view(qdotj,1:exogenous_dimension(joints[jid]))
end

"""
    unconstrained_velocity_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)

Returns a view of the unconstrained dof position(s), if any, of joint `jid` in
the global state vector `x`.
"""
function unconstrained_velocity_vector(x::AbstractVector,ls::RigidBodyMotion,jid::Int)
  @unpack joints = ls
  qdotj = velocity_vector(x,ls,jid)
  j = joints[jid]
  return view(qdotj,exogenous_dimension(j)+1:exogenous_dimension(j)+unconstrained_dimension(j))
end

"""
    deformation_vector(x::AbstractVector,ls::RigidBodyMotion,bid::Int)

Returns a view of the global state vector for a linked system containing only
the body surface positions of the body with id `bid`.
"""
function deformation_vector(x::AbstractVector,ls::RigidBodyMotion,bid::Int)
  @unpack deformation_indices = ls
  return view(x,deformation_indices[bid])
end

### motion_rhs! APIs ###

"""
    motion_rhs!(dxdt::AbstractVector,x::AbstractVector,p::Tuple{RigidBodyMotion,BodyList},t::Real)

Sets the right-hand side vector `dxdt` (mutating) for linked system `ls` of bodies `bl`, using the current state vector `x`,
the current time `t`.
"""
function motion_rhs!(dxdt::AbstractVector,x::AbstractVector,p::Tuple{RigidBodyMotion,Union{Body,BodyList}},t::Real)
    ls, bl = p
    @unpack a_edof_buffer, a_udof_buffer, exogenous_function! = ls
    fill!(a_udof_buffer,0.0)
    exogenous_function!(a_edof_buffer,x,ls,t)
    motion_rhs!(dxdt,x,t,a_edof_buffer,a_udof_buffer,ls,bl)
end


function motion_rhs!(dxdt::AbstractVector,x::AbstractVector,t::Real,a_edof::AbstractVector,a_udof::AbstractVector,ls::RigidBodyMotion,bl::BodyList)
    @unpack joints, deformations = ls
    for (jid,joint) in enumerate(joints)
        dxdt_j = view(dxdt,ls,jid;dimfcn=position_and_vel_dimension)
        x_j = view(x,ls,jid;dimfcn=position_and_vel_dimension)
        a_edof_j = view(a_edof,ls,jid;dimfcn=exogenous_dimension)
        a_udof_j = view(a_udof,ls,jid;dimfcn=unconstrained_dimension)
        joint_rhs!(dxdt_j,x_j,t,a_edof_j,a_udof_j,joint)
    end

    for (bid,b) in enumerate(bl)
        dxdt_b = deformation_vector(dxdt,ls,bid)
        dxdt_b .= deformation_velocity(b,deformations[bid],t)
    end
    nothing
end

function motion_rhs!(dxdt::AbstractVector,x::AbstractVector,t::Real,a_edof::AbstractVector,a_udof::AbstractVector,ls::RigidBodyMotion,b::Body)
  _check_for_only_one_body(ls)
  motion_rhs!(dxdt,x,t,a_edof,a_udof,ls,BodyList([b]))
end

### updating exogenous variables ###

"""
    zero_exogenous(ls::RigidBodyMotion)

Generate a zero vector of exogenous acclerations for system `ls`
"""
zero_exogenous(ls::RigidBodyMotion) = zero(ls.a_edof_buffer)


"""
    update_exogenous!(ls::RigidBodyMotion,a_edof::AbstractVector)

Mutates the exogenous buffer of `ls` with the supplied vector of exogenous
accelerations `a_edof`. This function is useful for supplying time-varying exogenous
values to the integrator from an outer loop.
"""
function update_exogenous!(ls::RigidBodyMotion,a_edof::AbstractVector)
  @assert length(ls.a_edof_buffer) == length(a_edof) "Incorrect length of exogenous vector"
  ls.a_edof_buffer .= a_edof
  nothing
end


"""
    rebase_from_inertial_to_reference(Xi_to_A::MotionTransform,x::AbstractVector,m::RigidBodyMotion,reference_body::Int)

If `Xi_to_A` represents a motion transform from the inertial coordinates to some
other coordinates A, such as that of a body, then this function returns a motion transform from
the reference body coordinates to coordinates A. It uses the current joint
state vector `x` to construct the transform list of the system `m`.
"""
function rebase_from_inertial_to_reference(Xi_to_A::MotionTransform{ND},x::AbstractVector,m::RigidBodyMotion,reference_body::Int) where {ND}
    Xl = body_transforms(x,m)
    Xi_to_ref = _reference_transform(Xl,ND,Val(reference_body))
    Xref_to_A = Xi_to_A*inv(Xi_to_ref)
    return Xref_to_A
end

### Surface velocity API ###

"""
    surface_velocity!(u::AbstractVector,v::AbstractVector,bl::BodyList,x::AbstractVector,m::RigidBodyMotion,t::Real[;axes=0,frame=axes,motion_part=:full])

Calculate the surface velocity components `u` and `v` for the points on bodies `bl`. The function evaluates
prescribed kinematics at time `t` and extracts non-prescribed (exogenous and unconstrained) velocities from
state vector `x`. There are three keyword arguments that can change the behavior of this function:
`axes` determines whether the components returned are expressed in the inertial coordinate
system (0, the default) or another body's coordinates (the body ID); `frame`
specifies that the motion should be measured relative to another reference body's velocity (by default, 0); and `motion_part` determines
whether the entire reference body velocity is used (`:full`, the default), only the angular part (`:angular`),
or only the linear part (`:linear`).
"""
function surface_velocity!(u::AbstractVector,v::AbstractVector,bl::BodyList,x::AbstractVector,m::RigidBodyMotion{ND},t::Real;
                           axes=0,frame=0,motion_part=:full) where {ND}
    @unpack deformations = m
    Xl = body_transforms(x,m)
    vl = body_velocities(x,t,m)
    Xref = _reference_transform(Xl,ND,Val(axes))
    Xframe = _reference_transform(Xl,ND,Val(frame))
    vref_0 = inv(Xframe)*_reference_velocity(vl,ND,motion_part,Val(frame))
    for (bid,b) in enumerate(bl)
        ub = view(u,bl,bid)
        vb = view(v,bl,bid)
        R = _rotation_transform(Xl[bid],Xref,Val(bid==axes))
        vrel = vl[bid] - Xl[bid]*vref_0    # relative velocity
        _surface_velocity!(ub,vb,b,vrel,deformations[bid],R,t)
    end
end

function surface_velocity!(u::AbstractVector,v::AbstractVector,b::Body,x::AbstractVector,m::RigidBodyMotion,t::Real;kwargs...)
    _check_for_only_one_body(m)
    surface_velocity!(u,v,BodyList([b]),x,m,t;kwargs...)
end

_reference_transform(ml,ND,::Val{0}) = MotionTransform{ND}()
_reference_transform(ml,ND,::Val{N}) where N = ml[N]

_reference_velocity(vl,ND,motion_part,::Val{0}) = PluckerMotion{ND}()
_reference_velocity(vl,ND,motion_part,::Val{N}) where N = _body_velocity(vl[N],Val(motion_part))

_rotation_transform(X::MotionTransform,Xref::MotionTransform,::Val{false}) = rotation_transform(Xref*inv(X))
_rotation_transform(X::MotionTransform{ND},Xref::MotionTransform,::Val{true}) where {ND} = MotionTransform{ND}()

_body_velocity(v::PluckerMotion,::Val{:full}) = v
_body_velocity(v::PluckerMotion,::Val{:angular}) = angular_only(v)
_body_velocity(v::PluckerMotion,::Val{:linear}) = linear_only(v)

"""
    velocity_in_body_coordinates_2d(x,y,v::AbstractPluckerMotionVector)

Evaluate the rigid-body velocity `v` at position ``(x,y)``, in the same
coordinate system as `v` itself.
"""
function velocity_in_body_coordinates_2d(x̃,ỹ,vb::AbstractPluckerMotionVector{2})
    Xb_to_p = MotionTransform(x̃,ỹ,0.0)
    Xb_to_p*vb
end

function velocity_in_inertial_coordinates_2d(x̃,ỹ,vb::AbstractPluckerMotionVector{2},Xb_to_0::MotionTransform{2})
    rotation_transform(Xb_to_0)*velocity_in_body_coordinates_2d(x̃,ỹ,vb)
end

function velocity_in_body_coordinates_2d!(u::AbstractVector,v::AbstractVector,b::Body,vb::AbstractPluckerMotionVector{2})
    for i in 1:numpts(b)
        vp = velocity_in_body_coordinates_2d(b.x̃[i],b.ỹ[i],vb)
        u[i], v[i] = vp.linear
    end
end

function _surface_velocity!(u::AbstractVector,v::AbstractVector,b::Body,vb::AbstractPluckerMotionVector{2},deformation::AbstractDeformationMotion,Rb_to_0::MotionTransform{2},t)
    u .= 0.0
    v .= 0.0
    _surface_velocity!(u,v,b,deformation,t)
    for i in 1:numpts(b)
        vp = velocity_in_body_coordinates_2d(b.x̃[i],b.ỹ[i],vb)
        vp += PluckerMotion{2}(linear=[u[i],v[i]]) # add the deformation to the rigid-body plucker vector
        vp = Rb_to_0*vp # rotate to the inertial system
        u[i], v[i] = vp.linear
    end
end



### update_body! API ###

"""
    update_body!(bl::Union{Body,BodyList},x::AbstractVector,m::RigidBodyMotion)

Update body/bodies in `bl` with the rigid-body motion `m` and state vector `x`.
"""
function update_body!(bl::BodyList,x::AbstractVector,m::RigidBodyMotion)
    @unpack deformations = m
    for (bid,b) in enumerate(bl)
        _update_body!(b,deformation_vector(x,m,bid),deformations[bid])
    end
    ml = body_transforms(x,m)
    update_body!(bl,ml)
    return bl
end

function update_body!(b::Body,x::AbstractVector,m::RigidBodyMotion)
    _check_for_only_one_body(m)
    @unpack deformations = m
    _update_body!(b,deformation_vector(x,m,1),deformations[1])
    T = body_transforms(x,m)[1]
    update_body!(b,T)
    return b
end

"""
    update_body!(bl::Union{Body,BodyList},x::AbstractVector,m::RigidBodyMotion)

Transform body/bodies in `bl` with the rigid-body motion `m` and state vector `x`.
"""
function transform_body!(bl::BodyList,x::AbstractVector,m::RigidBodyMotion)
    @unpack deformations = m
    for (bid,b) in enumerate(bl)
        _update_body!(b,deformation_vector(x,m,bid),deformations[bid])
    end
    ml = body_transforms(x,m)
    transform_body!(bl,ml)
    return bl
end

function transform_body!(b::Body,x::AbstractVector,m::RigidBodyMotion)
    _check_for_only_one_body(m)
    @unpack deformations = m
    _update_body!(b,deformation_vector(x,m,1),deformations[1])
    T = body_transforms(x,m)[1]
    transform_body!(b,T)
    return b
end

"""
    maxvelocity(b::Union{Body,BodyList},x::AbstractVector,m::RigidBodyMotion[,tmax=10,dt=0.05])

Search through the given motion state vector `x` and motion `m` applied to body `b` and return `(umax,i,t,bid)`,
the maximum velocity magnitude, the global index of the body point where it
occurs, the time at which it occurs, and the body on which it occurs.
"""
function maxvelocity(bl::BodyList,x::AbstractVector,m::RigidBodyMotion;tf=10.0,dt=0.05)
    u, v = zeros(numpts(bl)), zeros(numpts(bl))
    i = 1
    bid = 1
    umax = 0.0
    tmax = 0.0
    for t in 0.0:dt:tf
        surface_velocity!(u,v,bl,x,m,t)
        umag = sqrt.(u.^2+v.^2)
        umax_t, i_t = findmax(umag)
        bid_t, iloc_t = global_to_local_index(i_t,bl)
        if umax_t > umax
            umax, i, tmax, bid  = umax_t, i_t, t, bid_t
        end
    end
    return umax, i, tmax, bid
end

maxvelocity(b::Body,x::AbstractVector,m::RigidBodyMotion;kwargs...) = maxvelocity(BodyList([b]),x,m;kwargs...)
