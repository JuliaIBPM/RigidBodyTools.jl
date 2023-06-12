# Linked systems #

struct LinkedSystem{ND}
    nls :: Int
    nbody :: Int

    "List of linked bodies, starting with the one connected to inertial system"
    lslists :: Vector{Vector{Int}}
    joints :: Vector{Joint}
    parent_body :: Vector{Int}
    parent_joint :: Vector{Int}
    child_bodies :: Vector{Vector{Int}}
    child_joints :: Vector{Vector{Int}}
end


function LinkedSystem(joints::Vector{<:Joint},nbody::Int)
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
    @assert total_members == nbody "Not all bodies are members of linked lists"

    LinkedSystem{physical_dimension(joints)}(lscnt,nbody,lslists,joints,parent_body,parent_joint,child_bodies,child_joints)
end

function Base.show(io::IO, ls::LinkedSystem)
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

# For iteration purposes
number_of_linked_systems(ls::LinkedSystem) = ls.nls

"""
    first_body(lsid::Int,ls::LinkedSystem)

Return the index of the first body in linked system `lsid` in
the overall set of linked systems `ls`.  This body's parent is the
inertial coordinate system.
"""
function first_body(lsid::Int,ls::LinkedSystem)
  @unpack lslists, parent_joint = ls
  first(lslists[lsid])
end

"""
    first_joint(lsid::Int,ls::LinkedSystem)

Return the index of the first joint in linked system `lsid` in
the overall set of linked systems `ls`. This joint's parent is the
inertial coordinate system.
"""
function first_joint(lsid::Int,ls::LinkedSystem)
  @unpack lslists, parent_joint = ls
  parent_joint[first_body(lsid,ls)]
end


number_of_dofs(ls::LinkedSystem) = mapreduce(x -> number_of_dofs(x),+,ls.joints)
state_dimension(ls::LinkedSystem) = mapreduce(x -> state_dimension(x),+,ls.joints)


"""
    getrange(ls::LinkedSystem,jid::Int) -> Range

Return the subrange of indices in the global state vector
for the state corresponding to joint `jid` in linked system `ls`.
"""
function getrange(ls::LinkedSystem,jid::Int)
    @unpack joints = ls
    0 < jid <= length(joints) || error("Unavailable joint")
    first = 1
    j = 1
    while j < jid
        first += state_dimension(joints[j])
        j += 1
    end
    last = first+state_dimension(joints[jid])-1
    return first:last
end

"""
    view(q::AbstractVector,ls::LinkedSystem,jid::Int) -> SubArray

Provide a view of the range of values in vector `q` corresponding to the state
of the joint with index `jid` in a LinkedSystem `ls`.
"""
function Base.view(q::AbstractVector,ls::LinkedSystem,jid::Int)
    length(q) == state_dimension(ls) || error("Inconsistent size of data for viewing")
    return view(q,getrange(ls,jid))
end


"""
    linked_system_transform(q::AbstractVector,ls::LinkedSystem) -> MotionTransformList

Return a MotionTransformList by parsing the overall state vector `q`
into the individual joints.
"""
function linked_system_transform(q::AbstractVector,ls::LinkedSystem{ND}) where {ND}
    X = MotionTransform{ND}()
    ml = [deepcopy(X) for jb in 1:ls.nbody]
    for lsid in 1:number_of_linked_systems(ls)
      jid = first_joint(lsid,ls)
      _joint_descendants_transform!(ml,X,jid,q,ls)
    end
    return MotionTransformList(ml)
end

function _joint_descendants_transform!(ml,Xp::MotionTransform,jid::Int,q::AbstractVector,ls::LinkedSystem)
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
