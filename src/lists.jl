import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect,view,findall

export BodyList, MotionList, RigidTransformList, JointList, getrange

abstract type SetOfBodies end

const LISTS = [:BodyList, :Body],
              [:MotionList, :AbstractMotion],
              [:RigidTransformList, :RigidTransform]

"""
    BodyList([b1,b2,...])

Create a list of bodies
""" BodyList

"""
    MotionList([m1,m2,...])

Create a list of motions
""" MotionList

"""
    RigidTransformList([t1,t2,...])

Create a list of rigid transforms
""" RigidTransformList


for (listtype,listelement) in LISTS

  @eval struct $listtype <: SetOfBodies
      list :: Vector{$listelement}
  end

  @eval Base.eltype(f::$listtype) = Base.eltype(f.list)

  @eval $listtype() = $listtype($listelement[])

  @eval @propagate_inbounds getindex(A::$listtype, i::Int) = A.list[i]
  @eval @propagate_inbounds setindex!(A::$listtype, v::$listelement, i::Int) = A.list[i] = v
  @eval @propagate_inbounds getindex(A::$listtype, I...) = A.list[I...]
  @eval @propagate_inbounds setindex!(A::$listtype, v, I...) = A.list[I...] = v

  @eval iterate(A::$listtype) = iterate(A.list)
  @eval iterate(A::$listtype,I) = iterate(A.list,I)
  @eval size(A::$listtype) = size(A.list)
  @eval length(A::$listtype) = length(A.list)

  @eval findall(f::Function,A::$listtype) = findall(f,A.list)

  @eval push!(bl::$listtype,b::$listelement) = push!(bl.list,b)


end


numpts(bl::BodyList) = mapreduce(numpts,+,bl;init=0)

"""
    collect(bl::bodylist[,endpoints=false][,ref=false]) -> Vector{Float64}, Vector{Float64}

Collect the inertial-space coordinates of all of the Lagrange points comprising
the bodies in body list `bl` and return each assembled set of coordinates as a vector.
By default, `endpoints=false` and `ref=false`, which means this collects
the midpoints of segments in the inertial coordinates. If `endpoints=true`
it collects segment endpoints instead. If `ref=true` it collects
the coordinates in the body coordinate system.
"""
collect(bl::BodyList;endpoints=false,ref=false) = _collect(bl,Val(endpoints),Val(ref))

collect(body::Body;kwargs...) = collect(BodyList([body]);kwargs...)


function _collect(bl,::Val{false},::Val{false})
    xtmp = Float64[]
    ytmp = Float64[]
    for b in bl
        append!(xtmp,b.x)
        append!(ytmp,b.y)
    end
    return xtmp,ytmp
end

function _collect(bl::BodyList,::Val{true},::Val{false})
    xtmp = Float64[]
    ytmp = Float64[]
    for b in bl
        append!(xtmp,b.xend)
        append!(ytmp,b.yend)
    end
    return xtmp,ytmp
end

function _collect(bl,::Val{false},::Val{true})
    xtmp = Float64[]
    ytmp = Float64[]
    for b in bl
        append!(xtmp,b.x̃)
        append!(ytmp,b.ỹ)
    end
    return xtmp,ytmp
end

function _collect(bl::BodyList,::Val{true},::Val{true})
    xtmp = Float64[]
    ytmp = Float64[]
    for b in bl
        append!(xtmp,b.x̃end)
        append!(ytmp,b.ỹend)
    end
    return xtmp,ytmp
end



"""
    getrange(bl::BodyList,i::Int) -> Range

Return the subrange of indices in the global set of surface point data
corresponding to body `i` in a BodyList `bl`.
"""
function getrange(bl::BodyList,i::Int)
    i <= length(bl) || error("Unavailable body")
    first = 1
    j = 1
    while j < i
        first += numpts(bl[j])
        j += 1
    end
    last = first+numpts(bl[i])-1
    return first:last
end

"""
    getrange(bl::BodyList,ml::MotionList,i::Int) -> Range

Return the subrange of indices in the global vector of motion state data
corresponding to body `i` in a BodyList `bl` with corresponding motion list `ml`.
"""
function getrange(bl::BodyList,ml::MotionList,i::Int)
    i <= length(bl) || error("Unavailable body")
    first = 1
    j = 1
    while j < i
        first += length(motion_state(bl[j],ml[j]))
        j += 1
    end
    last = first+length(motion_state(bl[i],ml[i]))-1
    return first:last
end

"""
    view(f::AbstractVector,bl::BodyList,i::Int) -> SubArray

Provide a view of the range of values in vector `f` corresponding to the Lagrange
points of the body with index `i` in a BodyList `bl`.
"""
function Base.view(f::AbstractVector,bl::BodyList,i::Int)
    length(f) == numpts(bl) || error("Inconsistent size of data for viewing")
    return view(f,getrange(bl,i))
end

"""
    sum(f::AbstractVector,bl::BodyList,i::Int) -> Real

Compute a sum of the elements of vector `f` corresponding to body `i` in body
list `bl`.
"""
Base.sum(f::AbstractVector,bl::BodyList,i::Int) = sum(view(f,bl,i))

"""
    (tl::RigidTransformList)(bl::BodyList) -> BodyList

Carry out in-place transformations of each body in `bl` with the
corresponding transformation in `tl`.
"""
@inline function (tl::RigidTransformList)(bl::BodyList)
  length(tl) == length(bl) || error("Inconsistent lengths of lists")
  map((T,b) -> T(b),tl,bl)
end

"""
    vec(tl::RigidTransformList) -> Vector{Float64}

Returns a concatenation of length-3 vectors of the form [x,y,α] corresponding to the translation
and rotation specified by the given by the list of transforms `tl`.
"""
vec(tl::RigidTransformList) = mapreduce(vec,vcat,tl)

"""
    RigidTransformList(x::Vector)

Parses vector `x`, containing rigid-body configuration data ordered as
x, y, α for each body, into a `RigidTransformList`.
"""
function RigidTransformList(x::Vector{T}) where T <: Real
    nb, md = _length_and_mod(x)
    md == 0 || error("Invalid length of vector")
    first = 0
    tl = RigidTransformList()
    for i in 1:nb
        push!(tl,RigidTransform(view(x,first+1:first+CHUNK)))
        first += CHUNK
    end
    return tl
end

_length_and_mod(x::Vector{T}) where T <: Real = (n = length(x); return n ÷ CHUNK, n % CHUNK)

"""
    motion_velocity(bl::BodyList,ml::MotionList,t::Real) -> Vector

Return the aggregated velocity components (as a vector) of a `MotionList`
at the given time `t`.
"""
function motion_velocity(bl::BodyList,ml::MotionList,t::Real)
    u = Float64[]
    for (b,m) in zip(bl,ml)
      ui = motion_velocity(b,m,t)
      append!(u,ui)
    end
    return u
end

function motion_velocity(bl::BodyList,motion::AbstractMotion,t::Real)
    u = Float64[]
    for b in bl
      ui = motion_velocity(b,motion,t)
      append!(u,ui)
    end
    return u
end

"""
    motion_state(bl::BodyList,ml::MotionList)

Return the current state vector of body list `bl` associated with
motion list `ml`. It returns the aggregated state vectors
of each body.
"""
function motion_state(bl::BodyList,ml::MotionList)
    x = Float64[]
    length(bl) == length(ml) || error("body and motion lists are not the same length")
    for (b,m) in zip(bl,ml)
      xi = motion_state(b,m)
      append!(x,xi)
    end
    return x
end

function motion_state(bl::BodyList,motion::AbstractMotion)
    x = Float64[]
    for b in bl
      xi = motion_state(b,motion)
      append!(x,xi)
    end
    return x
end

"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     bl::BodyList,ml::Union{AbstractMotion,MotionList},t::Real[;inertial=true])

Assign the components of velocity `u` and `v` (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the MotionList `ml` for each body. If only one motion is specified in `ml`, then it is assumed that
  this is to be applied to all bodies.
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 bl::BodyList,ml::Union{AbstractMotion,MotionList},t::Real;kwargs...)

   for i in 1:length(bl)
      mli = isa(ml,AbstractMotion) ? ml : ml[i]
      surface_velocity!(view(u,bl,i),view(v,bl,i),bl[i],mli,t;kwargs...)
   end
   return u, v
end

"""
    surface_velocity(bl::BodyList,ml::Union{AbstractMotion,MotionList},t::Real[;inertial=true])

Return the components of rigid body velocity (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the MotionList `ml` for each body. If only one motion is specified in `ml`, then it is assumed that
  this is to be applied to all bodies.
"""
surface_velocity(bl::BodyList,ml::Union{AbstractMotion,MotionList},t::Real;kwargs...) =
    surface_velocity!(zeros(Float64,numpts(bl)),zeros(Float64,numpts(bl)),bl,ml,t;kwargs...)


"""
    update_body!(bl::BodyList,x::AbstractVector,ml::MotionList)

Update the bodies in list `bl` with the given motion state vector `x`.
The argument `ml` simply provides the information needed to parse
the vector into each body.
"""
function update_body!(bl::BodyList,x::AbstractVector,ml::MotionList)
    for i in 1:length(bl)
        update_body!(bl[i],x[getrange(bl,ml,i)],ml[i])
    end
    return bl
end

"""
    maxlistvelocity(bl::BodyList,ml::Union{AbstractMotion,List}[,tmax=100,dt=0.01])

Search through the given motions `ml` applied to bodies `bl` and return `(umax,i,t,bodyindex)`,
the maximum velocity magnitude, the index of the body points where it
occurs, the time at which it occurs, and the body index it occurs on.
"""
function maxlistvelocity(bl::BodyList,ml::Union{AbstractMotion,MotionList};kwargs...)
    i = 1
    umax = 0.0
    tmax = 0.0
    bmax = 1
    for j in 1:length(bl)
        mlj = isa(ml,AbstractMotion) ? ml : ml[j]
        umax_j,i_j,t_j = maxvelocity(bl[j],mlj,kwargs...)
        if umax_j > umax
            umax, i, tmax, bmax = umax_j, i_j, t_j, j
        end
    end
    return umax, i, tmax, bmax
end
