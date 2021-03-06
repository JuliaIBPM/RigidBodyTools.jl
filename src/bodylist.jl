import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect,view

export BodyList, RigidMotionList, RigidTransformList, DirectlySpecifiedMotionList,
          getrange, numpts

abstract type SetOfBodies end

const LISTS = [:BodyList, :Body],
              [:RigidMotionList, :RigidBodyMotion],
              [:RigidTransformList, :RigidTransform],
              [:DirectlySpecifiedMotionList,:DirectlySpecifiedMotion]



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

  @eval push!(bl::$listtype,b::$listelement) = push!(bl.list,b)


end


numpts(bl::BodyList) = mapreduce(numpts,+,bl)

"""
    collect(bl::bodylist) -> Vector{Float64}, Vector{Float64}

Collect the inertial-space coordinates of all of the Lagrange points comprising
the bodies in body list `bl` and return each assembled set of coordinates as a vector.
"""
function collect(bl::BodyList)
    xtmp = Float64[]
    ytmp = Float64[]
    for b in bl
        append!(xtmp,b.x)
        append!(ytmp,b.y)
    end
    return xtmp,ytmp
end
collect(body::Body) = collect(BodyList([body]))

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
        first += length(bl[j])
        j += 1
    end
    last = first+length(bl[i])-1
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
vec(tl::RigidTransformList) where {N} = mapreduce(vec,vcat,tl)

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
    rigidbodyvelocity(ml::RigidMotionList,t::Real) -> Vector

Return the velocity components (as a vector) of a `RigidMotionList`
at the given time `t`.
"""
function rigidbodyvelocity(ml::Union{RigidMotionList,DirectlySpecifiedMotionList},t::Real)
    u = Float64[]
    for m in ml
      ui = rigidbodyvelocity(m,t)
      append!(u,ui)
    end
    return u
end
