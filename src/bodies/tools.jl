### Tools ###

import Base:diff,length
export midpoints,dlength,normal,dlengthmid,centraldiff,normalmid


# Evaluate some geometric details of a body
"""
    length(body::Body)

Return the number of points on the body perimeter
"""
length(::Body{N}) where {N} = N

"""
    diff(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y differences of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
"""
diff(b::Body{N,C}) where {N,C<:BodyClosureType} = _diff(b.x,b.y,C)

function diff(bl::BodyList)
    dx = Float64[]
    dy = Float64[]
    for b in bl
        dxb, dyb = diff(b)
        append!(dx,dxb)
        append!(dy,dyb)
    end
    return dx, dy
end

function _diff(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})
  N = length(x)
  @assert N == length(y)

  ip1(i) = 1+mod(i,N)
  dxtmp = [x[ip1(i)] - x[i] for i = 1:N]
  dytmp = [y[ip1(i)] - y[i] for i = 1:N]

  return dxtmp, dytmp

end

function _diff(x::Vector{Float64},y::Vector{Float64},::Type{OpenBody})
  @assert length(x) == length(y)
  return diff(x), diff(y)

end

"""
    midpoints(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the x and y midpoints of the faces on the perimeter of body `body`, whose ends
are at the current `x` and `y` coordinates (in inertial space) of the body. Face 1
corresponds to the face between points 1 and 2, for example.

If `body` is a `BodyList`, then it computes the differences separately on each
constituent body.
"""
midpoints(b::Body{N,C}) where {N,C<:BodyClosureType} = _midpoints(b.x,b.y,C)

function midpoints(bl::BodyList)
    xc = Float64[]
    yc = Float64[]
    for b in bl
        xcb, ycb = midpoints(b)
        append!(xc,xcb)
        append!(yc,ycb)
    end
    return xc, yc
end

function _midpoints(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})

  N = length(x)
  @assert N == length(y)

  ip1(i) = 1+mod(i,N)
  xc = 0.5*[x[ip1(i)] + x[i] for i = 1:N]
  yc = 0.5*[y[ip1(i)] + y[i] for i = 1:N]

  return xc, yc

end

function _midpoints(x::Vector{Float64},y::Vector{Float64},::Type{OpenBody})

  @assert length(x) == length(y)

  xc = 0.5*(x[1:end-1]+x[2:end])
  yc = 0.5*(y[1:end-1]+y[2:end])

  return xc, yc

end

"""
    centraldiff(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the circular central differences of coordinates on body `body` (or
on each body in list `body`).
"""
centraldiff(b::Body{N,C}) where {N,C<:BodyClosureType} = _centraldiff(b.x,b.y,C)

function centraldiff(bl::BodyList)
    dx = Float64[]
    dy = Float64[]
    for b in bl
        dxb, dyb = centraldiff(b)
        append!(dx,dxb)
        append!(dy,dyb)
    end
    return dx, dy
end

function _centraldiff(x::Vector{Float64},y::Vector{Float64},::Type{ClosedBody})

  xc, yc = _midpoints(x,y,ClosedBody)
  xc .= circshift(xc,1)
  yc .= circshift(yc,1)

  return _diff(xc,yc,ClosedBody)
end

function _centraldiff(x::Vector{Float64},y::Vector{Float64},::Type{OpenBody})

  xc, yc = _midpoints(x,y,OpenBody)
  dxc, dyc = _diff(xc,yc,OpenBody)

  return [xc[1]-x[1];dxc;x[end]-xc[end]], [yc[1]-y[1];dyc;y[end]-yc[end]]
end

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
function dlengthmid(b::Union{Body,BodyList})
  dx, dy = centraldiff(b)
  return sqrt.(dx.^2+dy.^2)
end

"""
    normal(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals (in inertial components) of the faces on the perimeter
of body `body`, whose ends are at the current `x` and `y` coordinates (in inertial space)
of the body. Face 1 corresponds to the face between points 1 and 2, for example. For an `OpenBody`,
this provides a vector that is one element shorter than the number of points.
"""
function normal(b::Union{Body,BodyList})
  dx, dy = diff(b)
  ds = dlength(b)
  return -dy./ds, dx./ds
end

"""
    normalmid(body::Body/BodyList) -> Tuple{Vector{Float64},Vector{Float64}}

Compute the current normals (in inertial components) of the faces formed between
midpoints on the perimeter of body `body` (or each body in list `body`). For an `OpenBody`,
the normals for the first and last points are calculated for the face adjoining
with the adjacent midpoints. 
"""
function normalmid(b::Union{Body,BodyList})
  dx, dy = centraldiff(b)
  ds = dlengthmid(b)
  return -dy./ds, dx./ds
end
