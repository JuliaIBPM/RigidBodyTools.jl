### Tools ###

import Base:diff, length
export midpoints,dlength,normal,dlengthmid,centraldiff,normalmid,numpts,
        arccoord,arccoordmid,arclength,tangent,tangentmid,curvature

#=
A few notes on indexing of points on bodies.
- the face i has midpoint x[i], y[i], and endpoints with
  coordinates xend[i],yend[i],xend[i+1],yend[i+1]
- for a closed body, the dimensions of x,y and xend, yend are equal
- for an open body, the dimensions of x,y are one smaller than xend, yend
=#


# Evaluate some geometric details of a body
"""
    length(body::Body)

Return the number of points on the body perimeter
"""
length(::Body{N}) where {N} = N

@deprecate length numpts

for f in [:diff,:midpoints,:centraldiff]
    _f = Symbol("_"*string(f))
    @eval $f(b::Body;axes=:inertial) = $_f(b,Val(axes))

    @eval $_f(b::Body,::Val{:inertial}) = $_f(b,b.xend,b.yend)
    @eval $_f(b::Body,::Val{:body}) = $_f(b,b.x̃end,b.ỹend)

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


for f in [:diff,:midpoints]
  _f = Symbol("_"*string(f))
  _f_ccw = Symbol("_"*string(f)*"_ccw")

  @eval $_f(b::Body{N,C},x::AbstractVector{T},y::AbstractVector{T}) where {N,C<:BodyClosureType,T} = $_f(x,y,N,C) 
  @eval $_f_ccw(b::Body{N,C},x::AbstractVector{T},y::AbstractVector{T}) where {N,C<:BodyClosureType,T} = $_f_ccw(x,y,N,C) 

  @eval function $_f(bl::BodyList,x::AbstractVector{T},y::AbstractVector{T}) where {T}
    xl = Float64[]
    yl = Float64[]
    for (i,b) in enumerate(bl)
        xb, yb = $_f(b,view(x,bl,i),view(y,bl,i))
        append!(xl,xb)
        append!(yl,yb)
    end
    return xl, yl
  end

  @eval function $_f_ccw(bl::BodyList,x::AbstractVector{T},y::AbstractVector{T}) where {T}
    xl = Float64[]
    yl = Float64[]
    for (i,b) in enumerate(bl)
        xb, yb = $_f_ccw(b,view(x,bl,i),view(y,bl,i))
        append!(xl,xb)
        append!(yl,yb)
    end
    return xl, yl
  end

end


"""
    diff(body::Body[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y differences of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body (if `axes=:inertial`),
or at the reference `x̃` and `ỹ` coordinates (body-fixed space) if `axes=:body`. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
""" diff(::Body)

"""
    diff(bl::BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the `diff` on each constituent body in `bl`.
""" diff(::BodyList)


# This function is appropriate for differencing of face end points to
# compute midpoint-centered values.
function _diff(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{ClosedBody})
  N = length(x)
  @assert N == length(y)

  dxtmp = circshift(x,-1) .- x
  dytmp = circshift(y,-1) .- y

  return dxtmp, dytmp

end

# This function is appropriate for differencing of face midpoints points to
# compute endpoint-centered values.
function _diff_ccw(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{ClosedBody})
  N = length(x)
  @assert N == length(y)

  dxtmp = x .- circshift(x,1)
  dytmp = y .- circshift(y,1)

  return dxtmp, dytmp

end

function _diff(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{OpenBody})
  Nx = length(x)
  @assert Nx == length(y)
  _diff_openbody(x,y,Val(Nface==Nx))
end

# differencing of end point data to center point data -> length N+1 to N
function _diff_openbody(x,y,::Val{false})
  return diff(x), diff(y)
end

# differencing of center point data to end point data -> length N to N+1
function _diff_openbody(x,y,::Val{true})
  return [0.0;diff(x);0.0], [0.0;diff(y);0.0]
end

_diff_ccw(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{OpenBody}) = 
    _diff(x,y,Nface,OpenBody)


"""
    midpoints(body::Body[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y midpoints of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body (if `axes=:inertial`),
or at the reference `x̃` and `ỹ` coordinates (body-fixed space) if `axes=:body`. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
""" midpoints(::Body)

"""
    midpoints(bl::BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the `midpoints` on each constituent body in `bl`.
""" midpoints(::BodyList)

function _midpoints(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{ClosedBody})

  N = length(x)
  @assert N == length(y)

  xc = 0.5*(x .+ circshift(x,-1))
  yc = 0.5*(y .+ circshift(y,-1))

  return xc, yc

end

function _midpoints_ccw(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{ClosedBody})

  N = length(x)
  @assert N == length(y)

  xc = 0.5*(x .+ circshift(x,1))
  yc = 0.5*(y .+ circshift(y,1))

  return xc, yc

end

function _midpoints(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{OpenBody})
  Nx = length(x)
  @assert Nx == length(y)

  _midpoints_openbody(x,y,Val(Nx==Nface))

end

# averaging of end point data to center point data -> length Nface+1 to Nface
function _midpoints_openbody(x,y,::Val{false})
  xc = 0.5*(x[1:end-1]+x[2:end])
  yc = 0.5*(y[1:end-1]+y[2:end])
  return xc, yc
end

# averaging of center point data to end point data -> length Nface to Nface+1
function _midpoints_openbody(x,y,::Val{true})
  xc = 0.5*(x[1:end-1]+x[2:end])
  yc = 0.5*(y[1:end-1]+y[2:end])
  return [x[1];xc;x[end]], [y[1];yc;y[end]]
end

_midpoints_ccw(x::AbstractVector{Float64},y::AbstractVector{Float64},Nface,::Type{OpenBody}) = 
    _midpoints(x,y,Nface,OpenBody)

"""
    centraldiff(body::Body[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the circular central differences of coordinates on body `body` (or
on each body in list `body`). If `axes=:body`, uses the reference coordinates
in body-fixed space.
""" centraldiff(::Body)

"""
    centraldiff(bl::BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the `centraldiff` on each constituent body in `bl`.  If `axes=:body`, uses the reference coordinates
in body-fixed space.
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
    dlength(body::Body/BodyList[;axes=:inertial]) -> Vector{Float64}

Compute the lengths of the faces on the perimeter of body `body`, whose ends
are at the current `xend` and `yend` coordinates (in inertial space) of the body. Face 1
corresponds to the face between endpoints 1 and 2, for example. If `axes=:body`,
uses the reference coordinates in body-fixed space.
"""
function dlength(b::Union{Body,BodyList};kwargs...)
  dx, dy = diff(b;kwargs...)
  return sqrt.(dx.^2+dy.^2)
end


"""
    dlengthmid(body::Body/BodyList[;axes=:inertial]) -> Vector{Float64}

Same as [`dlength`](@ref).
"""
dlengthmid(b::Union{Body,BodyList};kwargs...) = dlength(b;kwargs...)
#=
function dlengthmid(b::Union{Body,BodyList})
  dx, dy = centraldiff(b)
  return sqrt.(dx.^2+dy.^2)
end
=#

"""
    tangent(body::Body/BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current tangent in inertial components (if `axes=:inertial`) or body-
  fixed components (if `axes=:body`) of the faces on the perimeter
of body `body`, whose ends are at the current `xend` and `yend` coordinates (in inertial space)
of the body. Face 1 corresponds to the face between points 1 and 2, for example. For an `OpenBody`,
this provides a vector that is one element shorter than the number of points.
"""
function tangent(b::Union{Body,BodyList};kwargs...)
  dx, dy = diff(b;kwargs...)
  ds = dlength(b)
  return dx./ds, dy./ds
end

"""
    tangentmid(body::Body/BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current tangents in inertial components (if `axes=:inertial`) or body-
  fixed components (if `axes=:body`) of the faces formed between
endpoints on the perimeter of body `body` (or each body in list `body`).
"""
tangentmid(b::Union{Body,BodyList};kwargs...) = tangent(b;kwargs...)


"""
    normal(body::Body/BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals in inertial components (if `axes=:inertial`) or body-
  fixed components (if `axes=:body`) of the faces on the perimeter
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
    normalmid(body::Body/BodyList[;axes=:inertial]) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals in inertial components (if `axes=:inertial`) or body-
  fixed components (if `axes=:body`) of the faces formed between
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

"""
    curvature(body::Body/BodyList[;axes=:inertial]) -> Vector{Float64}

Compute the current curvature of the faces on the perimeter
of body `body`. Face 1 corresponds to the face between points 1 and 2, for example. For an `OpenBody`,
this provides a vector that is one element shorter than the number of points.
"""
function curvature(b::Body;kwargs...)

  tx, ty = tangent(b;kwargs...) # tau, based on face centers
  dtx, dty = _diff(b,tx,ty) # dtau, based on face endpoints
  # dtx and dty "live" at face endpoints. For open bodies,
  dtxc, dtyc = _midpoints_ccw(b,dtx,dty) # average of adjacent dtau to face center

  ds = dlength(b) # lengths of faces

  return sqrt.(dtxc.^2 .+ dtyc.^2)./ds
end

function curvature(bl::BodyList;kwargs...)
  κl = Float64[]
  for b in bl
    κ = curvature(b;kwargs...)
    append!(κl,κ)
  end
  return κl
end

"""
    arccoordmid(body::Body/BodyList[;axes=:inertial]) -> Vector{Float64}

Returns a vector containing the arclength coordinate along the surface of `body`, evaluated at the
midpoints between the ends of faces. So, e.g., the first coordinate would be half
of the length of face 1, the second would be half of face 2 plus all of face 1,
etc. Use inertial components (if `axes=:inertial`) or body-
  fixed components (if `axes=:body`). If this is a body list, restart
  the origin of the coordinates on each body in the list.
"""
function arccoordmid(b::Body;kwargs...)
    ds = dlength(b;kwargs...)
    s = 0.5*ds
    pop!(pushfirst!(ds,0.0))
    s .+= accumulate(+,ds)
    return s
end

"""
    arccoord(body::Body/BodyList[;axes=:inertial]) -> Vector{Float64}

Returns a vector containing the arclength coordinate along the surface of `body`, evaluated at the
second endpoint of each face. So, e.g., the first coordinate would be the length
of face 1, the second the length of face 2, and the last would
be total length of all of the faces. Use inertial components (if `axes=:inertial`) or body-fixed components (if `axes=:body`).
If this is a body list, restart
the origin of the coordinates on each body in the list.
"""
function arccoord(b::Body;kwargs...)
    ds = dlength(b;kwargs...)
    s = accumulate(+,ds)
    return s
end

for f in [:arccoord,:arccoordmid]
    @eval function $f(bl::BodyList;kwargs...)
        sl = Float64[]
        for b in bl
            sb = $f(b;kwargs...)
            append!(sl,sb)
        end
        return sl
    end
end

"""
    arclength(body::Body[;axes=:inertial])

Compute the total arclength of `body`, from the sum of the lengths of the
faces. If `axes=:body`, use the body-fixed coordinates.
"""
arclength(b::Body;kwargs...) = sum(dlength(b;kwargs...))

"""
    arclength(bl::BodyList[;axes=:inertial]) -> Vector{Float64}

Compute the total arclength of each body in `bl` and assemble the
results into a vector. If `axes=:body`, use the body-fixed coordinates.
"""
function arclength(bl::BodyList;kwargs...)
    [arclength(b;kwargs...) for b in bl]
end
