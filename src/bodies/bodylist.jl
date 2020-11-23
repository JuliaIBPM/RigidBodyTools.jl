import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect,view

export BodyList, RigidMotionList, RigidTransformList, getrange, numpts

abstract type SetOfBodies end

const LISTS = [:BodyList, :Body],
              [:RigidMotionList, :RigidBodyMotion],
              [:RigidTransformList, :RigidTransform]



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



numpts(A::Body{N}) where {N} = N
numpts(A::BodyList) = mapreduce(b -> length(b.x),+,A)

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
