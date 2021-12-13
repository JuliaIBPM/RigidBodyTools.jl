using Dierckx
using Statistics: mean

using Elliptic
using Roots
using LinearAlgebra

#=
Approach:
- every body shape must contain xend, yend, x̃end, ỹend, x, y, x̃, ỹ
- end points of segments are moved (dx̃end/dt) transformed (x̃end -> xend).
- For polygons, the endpoints coincide with vertices.
- the Lagrange points (x, y) are the midpoints of segments, and are to be
  computed from (xend, yend) using midpoint after every transform.
- may not need the x̃, ỹ coordinates.
=#

export BasicBody,Ellipse,Circle,Rectangle,Square,Plate,ThickPlate,SplinedBody,NACA4,Polygon

"""
    BasicBody(x,y[,closuretype=ClosedBody]) <: Body

Construct a body by simply passing in the `x` and `y` coordinate vectors. The last
point will be automatically connected to the first point. The coordinate vectors
are assumed to be expressed in the body-fixed coordinate system. The optional
`closuretype` specifies whether the body is closed (`ClosedBody`) or open (`OpenBody`).
If closed, then the first and last points are assumed joined in operations that
require neighbor points.
"""
mutable struct BasicBody{N,C<:BodyClosureType} <: Body{N,C}
  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

  x̃end :: Vector{Float64}
  ỹend :: Vector{Float64}

  xend :: Vector{Float64}
  yend :: Vector{Float64}

end

function BasicBody(xend::Vector{T},yend::Vector{T};closuretype::Type{<:BodyClosureType}=ClosedBody) where {T <: Real}
    @assert length(xend) == length(yend)
    x, y = _midpoints(xend,yend,closuretype)
    BasicBody{length(x),closuretype}((0.0,0.0),0.0,x,y,x,y,xend,yend,xend,yend)
end

function Base.show(io::IO, body::BasicBody{N,C}) where {N,C}
    println(io, "Basic pointwise-specified body with $N points")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end



#### Ellipses and circles ####

"""
    Ellipse(a,b,n) <: Body

Construct an elliptical body with semi-major axis `a` and semi-minor axis `b`,
with `n` points distributed on the body perimeter.
"""
mutable struct Ellipse{N} <: Body{N,ClosedBody}
  a :: Float64
  b :: Float64
  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

  x̃end :: Vector{Float64}
  ỹend :: Vector{Float64}

  xend :: Vector{Float64}
  yend :: Vector{Float64}


end

#=
endpointson=false ensures that the midpoints lie on the surface
of the ellipse (and the endpoints are outside)
=#
function Ellipse(a::Real,b::Real,N::Int;endpointson=false,shifted=false)
    #N = 4*round(Int,N_given/4) # ensure that N is divisible by 4

    shiftedtype = shifted ? Shifted : Unshifted
    Nqrtr = div(N,4) + 1
    maj, min = a > b ? (a, b) : (b, a)
    m = 1 - min^2/maj^2
    smax = maj*Elliptic.E(m)
    Δs = smax/(Nqrtr-1)
    ϴ, _ = _get_ellipse_thetas(Δs,smax,maj,m,shiftedtype)
    _x, _y = min*cos.(ϴ), maj*sin.(ϴ)
    x, y = a > b ? (reverse(_y),reverse(_x)) : (_x,_y)
    adj = _last_segment_decrement(shiftedtype)
    xt = vcat(x,-reverse(x[1:end-adj]),-x[1+adj:end], reverse(x[1+adj:end-adj]))
    yt = vcat(y, reverse(y[1:end-adj]),-y[1+adj:end],-reverse(y[1+adj:end-adj]))

    # if this is endpointson = false, then xt, yt are the desired midpoints
    # otherwise, these are the endpoints
    if endpointson
      xend, yend = copy(xt), copy(yt)
      x, y = _midpoints(xt,yt,ClosedBody)
    else
       x, y = copy(xt), copy(yt)
       midinv = midpoint_inverse(length(x))
       xend = midinv*x
       yend = midinv*y
    end

    Ellipse{length(x)}(a,b,(0.0,0.0),0.0,x,y,x,y,xend,yend,xend,yend)
end

function _get_ellipse_thetas(Δs,smax,maj,m,shiftedtype)

    slast = _start_s(Δs,shiftedtype)
    θlast = 0.0
    s = maj*Elliptic.E(θlast,m)
    θ = _init_Θ(shiftedtype)
    ds = []
    while s < _max_s(smax,Δs,shiftedtype)
        θi = find_zero(x -> maj*Elliptic.E(x,m) - slast - Δs,θlast,atol=1e-12)
        s = maj*Elliptic.E(θi,m)
        push!(ds,s-slast)
        push!(θ,θi)
        θlast = θi
        slast = s
    end
    θ, ds
end

_start_s(Δs,::Type{Shifted}) = -0.5*Δs
_start_s(Δs,::Type{Unshifted}) = 0.0

_init_Θ(::Type{Shifted}) = Float64[]
_init_Θ(::Type{Unshifted}) = Float64[0.0]

_max_s(smax,Δs,::Type{Shifted}) = smax-Δs
_max_s(smax,Δs,::Type{Unshifted}) = smax

_last_segment_decrement(::Type{Shifted}) = 0
_last_segment_decrement(::Type{Unshifted}) = 1

midpoint_inverse(N) = _midpoint_inverse(N,Val(mod(N,2)==0))

function _midpoint_inverse(n::Int,::Val{true}) # even N
  #mod(n,2) == 0 || error("midpoint inverse only verified for even n")
  d = 1/n
  u = Matrix(undef,1,n)
  num = 1.0-d
  sgn = 1.0
  u[1] = u[n] = sgn*num
  for i in 2:n÷2
      num -= 2d
      sgn *= -1.0
      u[i] = u[n-i+1] = sgn*num
  end
  A = Matrix(undef,n,n)
  A[1,:] = u
  for i in 2:n
    u .= circshift(u,(0,1))
    A[i,:] = u
  end
  return A
end

function _midpoint_inverse(n::Int,::Val{false}) # odd N
  #mod(n,2) == 0 || error("midpoint inverse only verified for even n")
  u = Matrix(undef,1,n)
  sgn = 1.0
  u[1] = sgn
  for i in 2:n
      sgn *= -1.0
      u[i] = sgn
  end
  A = copy(u)
  for i in 2:n
    u .= circshift(u,(0,1))
    A = vcat(A,u)
  end
  return A
end

"""
    Circle(a,n) <: Body

Construct a circular body with radius `a`
and with `n` points distributed on the body perimeter.
""" Circle(::Real,::Int)

"""
    Circle(a,targetsize::Float64) <: Body

Construct a circular body with radius `a` with spacing between points set
approximately to `targetsize`.
""" Circle(::Real,::Real)

Circle(a::Real,arg;kwargs...) = Ellipse(a,a,arg;kwargs...)


function Base.show(io::IO, body::Ellipse{N}) where {N}
    if body.a == body.b
      println(io, "Circular body with $N points and radius $(body.a)")
    else
      println(io, "Elliptical body with $N points and semi-axes ($(body.a),$(body.b))")
    end
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end

#### Rectangles and squares ####

"""
    Rectangle(a,b,n) <: Body

Construct a rectangular body with x̃ side half-length `a` and ỹ side half-length `b`,
with approximately `n` points distributed along the perimeter. The centroid
of the rectangle is placed at the origin (so that the lower left corner is at (-a,-b)).

Points are not placed at the corners, but rather, are shifted
by half a segment. This ensures that all normals are perpendicular to the sides.
"""
Rectangle(a::Real,b::Real,arg;kwargs...) = Polygon([-a,a,a,-a],[-b,-b,b,b],arg;kwargs...)

"""
    Rectangle(a,b,ds) <: Body

Construct a rectangular body with x̃ side half-length `a` and ỹ side half-length `b`,
with approximate spacing `ds` between points.
""" Rectangle(::Real,::Real,::Real)

"""
    Square(a,n) <: Body

Construct a square body with side half-length `a`
and with approximately `n` points distributed along the perimeter.
"""
Square(a::Real,arg;kwargs...) = Rectangle(a,a,arg;kwargs...)

"""
    Square(a,ds) <: Body

Construct a square body with side half-length `a`,
with approximate spacing `ds` between points.
""" Square(::Real,::Real)


#### Plates ####

"""
    Plate(a,n) <: Body

Construct a flat plate of zero thickness with length `a`,
divided into `n` equal segments.
"""
Plate(len::Real,arg) = Polygon([-0.5*len,0.5*len],[0.0,0.0],arg,closuretype=OpenBody)

"""
    Plate(a,ds) <: Body

Construct a flat plate of zero thickness with length `a`,
with approximate spacing `ds` between points.
""" Plate(::Real,::Real)



#### Polygons ####

"""
    Polygon(x::Vector,y::Vector,n[,closuretype=ClosedBody])

Create a polygon shape with vertices `x` and `y`, with approximately `n` points distributed along
the perimeter.
"""
mutable struct Polygon{N,NV,C<:BodyClosureType} <: Body{N,C}
  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

  x̃end :: Vector{Float64}
  ỹend :: Vector{Float64}

  xend :: Vector{Float64}
  yend :: Vector{Float64}

end

"""
    Polygon(x::Vector,y::Vector,ds::Float64[,closuretype=ClosedBody])

Create a polygon shape with vertices `x` and `y`, with approximate spacing `ds` between points.
""" Polygon(::AbstractVector,::AbstractVector,::Real)


function Polygon(xv::AbstractVector{T},yv::AbstractVector{T},a::Float64;closuretype=ClosedBody) where {T<:Real}

  xend, yend, x, y  = _polygon(xv,yv,a,closuretype)
  Polygon{length(x),length(xv),closuretype}((0.0,0.0),0.0,x,y,x,y,xend,yend,xend,yend)
end

@inline Polygon(xv,yv,n::Int;closuretype=ClosedBody) =
            Polygon(xv,yv,polygon_perimeter(xv,yv,closuretype)/n,closuretype=closuretype)

function polygon_perimeter(xv,yv,closuretype)
  xvcirc = _extend_array(xv,closuretype)
  yvcirc = _extend_array(yv,closuretype)
  len = 0.0
  for i in 1:length(xvcirc)-1
     len += sqrt((xvcirc[i+1]-xvcirc[i])^2+(yvcirc[i+1]-yvcirc[i])^2)
  end
  return len
end

function Base.show(io::IO, body::Polygon{N,NV,CS}) where {N,NV,CS}
    cb = CS == ClosedBody ? "Closed " : "Open "
    println(io, cb*"polygon with $NV vertices and $N points")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end

function _polygon(xv::AbstractVector{T},yv::AbstractVector{T},ds::Float64,closuretype) where {T<:Real}
    xvcirc = _extend_array(xv,closuretype)
    yvcirc = _extend_array(yv,closuretype)
    xend, yend = Float64[], Float64[]
    for i in 1:length(xvcirc)-2
        xi, yi = _line_points(xvcirc[i],yvcirc[i],xvcirc[i+1],yvcirc[i+1],ds)
        append!(xend,xi[1:end-1])
        append!(yend,yi[1:end-1])
    end
    xi, yi = _line_points(xvcirc[end-1],yvcirc[end-1],xvcirc[end],yvcirc[end],ds)
    append!(xend,xi[1:end-_last_segment_decrement(closuretype)])
    append!(yend,yi[1:end-_last_segment_decrement(closuretype)])

    x, y = _midpoints(xend,yend,closuretype)

    return xend, yend, x, y
end

_extend_array(x,::Type{ClosedBody}) = [x;x[1]]
_extend_array(x,::Type{OpenBody}) = x

_last_segment_decrement(::Type{ClosedBody}) = 1
_last_segment_decrement(::Type{OpenBody}) = 0

function _line_points(x1,y1,x2,y2,n::Int)

  dx, dy = x2-x1, y2-y1
  len = sqrt(dx^2+dy^2)

  Δs = len/(n-1)

  x = zeros(n)
  y = zeros(n)

  range = (0:n-1)/(n-1)

  @. x = x1 + dx*range
  @. y = y1 + dy*range

  return x, y

end

function _line_points(x1,y1,x2,y2,ds::Float64)
    dx, dy = x2-x1, y2-y1
    len = sqrt(dx^2+dy^2)
    return _line_points(x1,y1,x2,y2,round(Int,len/ds)+1)
end



#### Splined body ####


"""
    SplinedBody(X,Δx[,closuretype=ClosedBody]) -> BasicBody

Using control points in `X` (assumed to be N x 2, where N is the number of points), create a set of points
that are uniformly spaced (with spacing `Δx`) on a curve that passes through the control points. A cubic
parametric spline algorithm is used. If the optional parameter `closuretype` is set
to `OpenBody`, then the end points are not joined together.
"""
function SplinedBody(Xpts_raw::Array{Float64,2},Δx::Float64;closuretype::Type{<:BodyClosureType}=ClosedBody)
    # Assume Xpts are in the form N x 2
    Xpts = copy(Xpts_raw)
    if Xpts[1,:] != Xpts[end,:]
        Xpts = vcat(Xpts,Xpts[1,:]')
    end

    spl = (closuretype == ClosedBody) ? ParametricSpline(Xpts',periodic=true) : ParametricSpline(Xpts')
    tfine = range(0,1,length=1001)
    dX = derivative(spl,tfine)

    np = ceil(Int,sqrt(sum(dX.^2)*(tfine[2]-tfine[1]))/Δx)

    tsamp = range(0,1,length=np)
    x = [X[1] for X in spl.(tsamp[1:end-1])]
    y = [X[2] for X in spl.(tsamp[1:end-1])]

    return BasicBody(x,y,closuretype=closuretype)
end

"""
    SplinedBody(x,y,Δx[,closuretype=ClosedBody]) -> BasicBody

Using control points in `x` and `y`, create a set of points
that are uniformly spaced (with spacing `Δx`) on a curve that passes through the control points. A cubic
parametric spline algorithm is used. If the optional parameter `closuretype` is set
to `OpenBody`, then the end points are not joined together.
"""
SplinedBody(x::AbstractVector{Float64},y::AbstractVector{Float64},Δx::Float64;kwargs...) =
      SplinedBody(hcat(x,y),Δx;kwargs...)


"""
    ThickPlate(length,thick,n,[λ=1.0]) <: Body

Construct a flat plate with length `length` and thickness `thick`,
with `n` points distributed along one side.

The optional parameter `λ` distributes the points differently. Values between `0.0`
and `1.0` are accepted.

Alternatively, the shape can be specified with a target point spacing in place of `n`.
"""
mutable struct ThickPlate{N} <: Body{N,ClosedBody}
  len :: Float64
  thick :: Float64
  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

  x̃end :: Vector{Float64}
  ỹend :: Vector{Float64}

  xend :: Vector{Float64}
  yend :: Vector{Float64}

end

function ThickPlate(len::Real,thick::Real,N::Int;λ::Float64=1.0)
    # input N is the number of panels on one side only

    # set up points on flat sides
    Δϕ = π/N
    Jϕa = [sqrt(sin(ϕ)^2+λ^2*cos(ϕ)^2) for ϕ in range(π-Δϕ/2,stop=Δϕ/2,length=N)]
    Jϕ = len*Jϕa/Δϕ/sum(Jϕa)
    xtopface = -0.5*len .+ Δϕ*cumsum([0.0; Jϕ])
    xtop = 0.5*(xtopface[1:N] + xtopface[2:N+1])

    Δsₑ = Δϕ*Jϕ[1]
    Nₑ = 2*floor(Int,0.25*π*thick/Δsₑ)
    xedgeface = [0.5*len + 0.5*thick*cos(ϕ) for ϕ in range(π/2,stop=-π/2,length=Nₑ+1)]
    yedgeface = [          0.5*thick*sin(ϕ) for ϕ in range(π/2,stop=-π/2,length=Nₑ+1)]
    xedge = 0.5*(xedgeface[1:Nₑ]+xedgeface[2:Nₑ+1])
    yedge = 0.5*(yedgeface[1:Nₑ]+yedgeface[2:Nₑ+1])

    x̃end = Float64[]
    ỹend = Float64[]
    for xi in xtop
      push!(x̃end,xi)
      push!(ỹend,0.5*thick)
    end
    for i = 1:Nₑ
      push!(x̃end,xedge[i])
      push!(ỹend,yedge[i])
    end
    for xi in reverse(xtop,dims=1)
      push!(x̃end,xi)
      push!(ỹend,-0.5*thick)
    end
    for i = Nₑ:-1:1
      push!(x̃end,-xedge[i])
      push!(ỹend,yedge[i])
    end
    x̃, ỹ = _midpoints(x̃end,ỹend,ClosedBody)

    ThickPlate{length(x̃)}(len,thick,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ,x̃end,ỹend,x̃end,ỹend)

end

function Base.show(io::IO, body::ThickPlate{N}) where {N}
    println(io, "Thick plate with $N points and length $(body.len) and thickness $(body.thick)")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end

#### NACA 4-digit airfoil ####

"""
    NACA4(cam,pos,thick,np,[len=1.0]) <: Body{N}

Generates points in the shape of a NACA 4-digit airfoil of chord length 1. The
relative camber is specified by `cam`, the position of
maximum camber (as fraction of chord) by `pos`, and the relative thickness
by `thick`. The parameter `np` specifies the number of points on the upper
or lower surface. The optional parameter `len` specifies the chord length,
which defaults to 1.0.

# Example

```jldoctest
julia> b = NACA4(0.0,0.0,0.12);
```
"""
mutable struct NACA4{N} <: Body{N,ClosedBody}
  len :: Float64
  camber :: Float64
  pos :: Float64
  thick :: Float64

  cent :: Tuple{Float64,Float64}
  α :: Float64

  x̃ :: Vector{Float64}
  ỹ :: Vector{Float64}

  x :: Vector{Float64}
  y :: Vector{Float64}

  x̃end :: Vector{Float64}
  ỹend :: Vector{Float64}

  xend :: Vector{Float64}
  yend :: Vector{Float64}

end


function NACA4(cam::Real,pos::Real,t::Real,np::Int;len=1.0)

# Here, cam is the fractional camber, pos is the fractional chordwise position
# of max camber, and t is the fractional thickness.

npan = 2*np-2

# Trailing edge bunching
an = 1.5
anp = an+1
x = zeros(np)

θ = zeros(size(x))
yc = zeros(size(x))

for j = 1:np
    frac = Float64((j-1)/(np-1))
    x[j] = 1 - anp*frac*(1-frac)^an-(1-frac)^anp;
    if x[j] < pos
        yc[j] = cam/pos^2*(2*pos*x[j]-x[j]^2)
        if pos > 0
            θ[j] = atan(2*cam/pos*(1-x[j]/pos))
        end
    else
        yc[j] = cam/(1-pos)^2*((1-2*pos)+2*pos*x[j]-x[j]^2)
        if pos > 0
            θ[j] = atan(2*cam*pos/(1-pos)^2*(1-x[j]/pos))
        end
    end
end

xu = zeros(size(x))
yu = xu
xl = xu
yl = yu

yt = t/0.20*(0.29690*sqrt.(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4)

xu = x-yt.*sin.(θ)
yu = yc+yt.*cos.(θ)

xl = x+yt.*sin.(θ)
yl = yc-yt.*cos.(θ)

rt = 1.1019*t^2;
θ0 = 0
if pos > 0
    θ0 = atan(2*cam/pos)
end
# Center of leading edge radius
xrc = rt*cos(θ0)
yrc = rt*sin(θ0)
θle = collect(0:π/50:2π)
xlec = xrc .+ rt*cos.(θle)
ylec = yrc .+ rt*sin.(θle)

# Assemble data
coords = [xu yu xl yl x yc]
cole = [xlec ylec]

# Close the trailing edge
xpanold = [0.5*(xl[np]+xu[np]); reverse(xl[2:np-1],dims=1); xu[1:np-1]]
ypanold = [0.5*(yl[np]+yu[np]); reverse(yl[2:np-1],dims=1); yu[1:np-1]]

xpan = zeros(npan)
ypan = zeros(npan)
for ipan = 1:npan
    if ipan < npan
        xpan1 = xpanold[ipan]
        ypan1 = ypanold[ipan]
        xpan2 = xpanold[ipan+1]
        ypan2 = ypanold[ipan+1]
    else
        xpan1 = xpanold[npan]
        ypan1 = ypanold[npan]
        xpan2 = xpanold[1]
        ypan2 = ypanold[1]
    end
    xpan[ipan] = 0.5*(xpan1+xpan2)
    ypan[ipan] = 0.5*(ypan1+ypan2)
end
w = ComplexF64[1;reverse(xpan,dims=1)+im*reverse(ypan,dims=1)]*len
w .-= mean(w)

x̃end = real.(w)
ỹend = imag.(w)

x̃, ỹ = _midpoints(x̃end,ỹend,ClosedBody)


NACA4{length(x̃)}(len,cam,pos,t,(0.0,0.0),0.0,x̃,ỹ,x̃,ỹ,x̃end,ỹend,x̃end,ỹend)

end


function Base.show(io::IO, body::NACA4{N}) where {N}
    println(io, "NACA 4-digit airfoil with $N points and length $(body.len) and thickness $(body.thick)")
    println(io, "   Current position: ($(body.cent[1]),$(body.cent[2]))")
    println(io, "   Current angle (rad): $(body.α)")
end


####


function _adjustnumber(targetsize::Real,shapefcn::Type{T},params...;kwargs...) where {T <: Body}
    ntrial = 501
    return floor(Int,ntrial*mean(dlength(shapefcn(params...,ntrial;kwargs...)))/targetsize)
end


for shape in (:Ellipse,:ThickPlate,:NACA4)
    @eval RigidBodyTools.$shape(params...;kwargs...) =
        $shape(Base.front(params)...,
          _adjustnumber(Base.last(params),$shape,Base.front(params)...;kwargs...);kwargs...)
end
