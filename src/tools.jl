### Tools ###

import Base:diff,length
export midpoints,dlength,normal,dlengthmid,centraldiff,normalmid


# Evaluate some geometric details of a body
"""
    length(body::Body)

Return the number of points on the body perimeter
"""
length(::Body{N}) where {N} = N

for f in [:diff,:midpoints,:centraldiff]
    _f = Symbol("_"*string(f))
    @eval $f(b::Body;ref=false) = $_f(b,Val(ref))

    @eval $_f(b::Body{N,C},::Val{false}) where {N,C<:BodyClosureType} = $_f(b.xend,b.yend,C)
    @eval $_f(b::Body{N,C},::Val{true}) where {N,C<:BodyClosureType}  = $_f(b.x̃end,b.ỹend,C)

    @eval function $f(bl::BodyList;kwargs...)
        xl = Float64[]
        yl = Float64[]
        for b in bl
            xb, yb = $f(b;kwargs...)
            append!(xl,xb)
            append!(yl,yb)
        end
        return xl, yl
    end
end

"""
    diff(body::Body[,ref=false]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y differences of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body (if `ref=false`),
or at the reference `x̃` and `ỹ` coordinates (body-fixed space) if `ref=true`. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
""" diff(::Body)

"""
    diff(bl::BodyList[,ref=false]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the `diff` on each constituent body in `bl`.
""" diff(::BodyList)


function _diff(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})
  N = length(x)
  @assert N == length(y)

  dxtmp = circshift(x,-1) .- x
  dytmp = circshift(y,-1) .- y

  return dxtmp, dytmp

end

function _diff_ccw(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})
  N = length(x)
  @assert N == length(y)

  dxtmp = x .- circshift(x,1)
  dytmp = y .- circshift(y,1)

  return dxtmp, dytmp

end

function _diff(x::Vector{Float64},y::Vector{Float64},::Type{OpenBody})
  @assert length(x) == length(y)
  return diff(x), diff(y)

end

"""
    midpoints(body::Body[,ref=false]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y midpoints of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body (if `ref=false`),
or at the reference `x̃` and `ỹ` coordinates (body-fixed space) if `ref=true`. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
""" midpoints(::Body)

"""
    midpoints(bl::BodyList[,ref=false]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the `midpoints` on each constituent body in `bl`.
""" midpoints(::BodyList)

function _midpoints(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})

  N = length(x)
  @assert N == length(y)

  xc = 0.5*(x .+ circshift(x,-1))
  yc = 0.5*(y .+ circshift(y,-1))

  return xc, yc

end

function _midpoints_ccw(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})

  N = length(x)
  @assert N == length(y)

  xc = 0.5*(x .+ circshift(x,1))
  yc = 0.5*(y .+ circshift(y,1))

  return xc, yc

end

function _midpoints(x::Vector{Float64},y::Vector{Float64},::Type{OpenBody})

  @assert length(x) == length(y)

  xc = 0.5*(x[1:end-1]+x[2:end])
  yc = 0.5*(y[1:end-1]+y[2:end])

  return xc, yc

end

"""
    centraldiff(body::Body) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the circular central differences of coordinates on body `body` (or
on each body in list `body`).
""" centraldiff(::Body)

"""
    centraldiff(bl::BodyList[,ref=false]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the `centraldiff` on each constituent body in `bl`.
""" centraldiff(::BodyList)

_centraldiff(x,y,C) = _diff(x,y,C)
#=
function _centraldiff(xend::Vector{Float64},yend::Vector{Float64},::Type{ClosedBody})

  xc, yc = _midpoints(x,y,ClosedBody)
  xc .= circshift(xc,1)
  yc .= circshift(yc,1)

  return _diff(xend,yend,ClosedBody)
end

function _centraldiff(xend::Vector{Float64},yend::Vector{Float64},::Type{OpenBody})

  xc, yc = _midpoints(x,y,OpenBody)
  dxc, dyc = _diff(xc,yc,OpenBody)

  return [xc[1]-x[1];dxc;x[end]-xc[end]], [yc[1]-y[1];dyc;y[end]-yc[end]]
end
=#

"""
    dlength(body::Body/BodyList) -> Vector{Float64}

Compute the lengths of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body. Face 1
corresponds to the face between points 1 and 2, for example. For an `OpenBody`,
this provides a vector that is one element shorter than the number of points, to
ensure that `sum(dlength(body))` is equal to the arclength of the body.
"""
function dlength(b::Union{Body,BodyList})
  dx, dy = diff(b)
  return sqrt.(dx.^2+dy.^2)
end


"""
    dlengthmid(body::Body/BodyList) -> Vector{Float64}

Compute the lengths of the faces formed between the face midpoints on the
perimeter of body `body`. The indexing of these midpoint faces is consistent
with that of the regular vertex points adjacent to both midpoints.
Midpoint face 2 corresponds to the face between midpoints 1 and 2, for example.
For an `OpenBody`, the lengths for the first and last points are calculated
to the adjoining midpoints, to
ensure that `sum(dlength(body))` is equal to the arclength of the body.
"""
dlengthmid(b::Union{Body,BodyList}) = dlength(b)
#=
function dlengthmid(b::Union{Body,BodyList})
  dx, dy = centraldiff(b)
  return sqrt.(dx.^2+dy.^2)
end
=#

"""
    normal(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals in inertial components (if `ref=false`) or body-
  fixed components (if `ref=true`) of the faces on the perimeter
of body `body`, whose ends are at the current `xend` and `yend` coordinates (in inertial space)
of the body. Face 1 corresponds to the face between points 1 and 2, for example. For an `OpenBody`,
this provides a vector that is one element shorter than the number of points.
"""
function normal(b::Union{Body,BodyList};kwargs...)
  dx, dy = diff(b;kwargs...)
  ds = dlength(b)
  return dy./ds, -dx./ds
end

"""
    normalmid(body::Body/BodyList[,ref=false]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals in inertial components (if `ref=false`) or body-
  fixed components (if `ref=true`) of the faces formed between
endpoints on the perimeter of body `body` (or each body in list `body`).
"""
normalmid(b::Union{Body,BodyList};kwargs...) = normal(b;kwargs...)

#=
function normalmid(b::Union{Body,BodyList};kwargs...)
  dx, dy = centraldiff(b;kwargs...)
  ds = dlengthmid(b)
  return dy./ds, -dx./ds
end
=#
