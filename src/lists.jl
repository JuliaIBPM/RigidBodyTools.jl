

export BodyList, RigidTransformList, MotionTransformList, PluckerMotionList, getrange

abstract type SetOfBodies end

const LISTS = [:BodyList, :Body],
              [:RigidTransformList, :RigidTransform],
              [:MotionTransformList, :MotionTransform],
              [:PluckerMotionList, :PluckerMotion]

"""
    BodyList([b1,b2,...])

Create a list of bodies
""" BodyList


"""
    RigidTransformList([t1,t2,...])

Create a list of rigid transforms
""" RigidTransformList

"""
    MotionTransformList([t1,t2,...])

Create a list of motion transforms
""" MotionTransformList

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

zero_body(b::Union{Body,BodyList}) = zeros(Float64,numpts(b))


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
    getrange(bl::BodyList,bid::Int) -> Range

Return the subrange of indices in the global set of surface point data
corresponding to body `bid` in a BodyList `bl`.
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
    global_to_local_index(i::Int,bl::BodyList) -> (Int, Int)

Return the ID `bid` of the body in body list `bl` on which global index `i` sits,
as well as the local index of `iloc` on that body and return as `(bid,iloc)`
"""
function global_to_local_index(i::Int,bl::BodyList)
    bid = 1
    ir = getrange(bl,bid)
    while !(in(i,ir)) && bid <= length(bl)
        bid += 1
        ir = getrange(bl,bid)
    end

    return bid, i-first(ir)+1
end

#=
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
=#

"""
    view(f::AbstractVector,bl::BodyList,bid::Int) -> SubArray

Provide a view of the range of values in vector `f` corresponding to the Lagrange
points of the body with index `bid` in a BodyList `bl`.
"""
function Base.view(f::AbstractVector,bl::BodyList,bid::Int)
    length(f) == numpts(bl) || error("Inconsistent size of data for viewing")
    return view(f,getrange(bl,bid))
end

"""
    sum(f::AbstractVector,bl::BodyList,i::Int) -> Real

Compute a sum of the elements of vector `f` corresponding to body `i` in body
list `bl`.
"""
Base.sum(f::AbstractVector,bl::BodyList,i::Int) = sum(view(f,bl,i))

"""
    (tl::MotionTransformList)(bl::BodyList) -> BodyList

Carry out transformations of each body in `bl` with the
corresponding transformation in `tl`, creating a new body list.
"""
@inline function (tl::MotionTransformList)(bl::BodyList)
  bl_transform = deepcopy(bl)
  update_body!(bl_transform,tl)
end

"""
    update_body!(bl::BodyList,tl::MotionTransformList) -> BodyList

Carry out in-place transformations of each body in `bl` with the
corresponding transformation in `tl`.
"""
@inline function update_body!(bl::BodyList,tl::MotionTransformList)
  length(tl) == length(bl) || error("Inconsistent lengths of lists")
  map((T,b) -> update_body!(b,T),tl,bl)
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
