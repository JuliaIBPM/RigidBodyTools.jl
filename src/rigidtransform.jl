# Rigid-body transformation routines

export RigidTransform

"""
    RigidTransform(x::Tuple{Float64,Float64},α::Float64)

Construct a rigid-body transform operator, with rotation by angle `α` and
translation specified by `x`. The translation coordinates are specified in the
target coordinate system.

The resulting transform can be used as an operator on pairs of coordinate vectors,
`x` and `y`, or on bodies. For transformation of bodies, it only overwrites the
`x` and `y` fields of the body, but leaves the `x̃` and `ỹ` (body coordinates) intact.

The translation can be provided as either a tuple `(x,y)` or as a complex number.

# Constructors

- `RigidTransform((x,y),α)`
- `RigidTransform(u::Vector{Float64})`
- `RigidTransform(u::NTuple{3,Float64})`
- `RigidTransform.(u)` where `u` is a collection of vectors or tuples.

# Example

```jldoctest
julia> body = Bodies.Ellipse(0.5,0.1,100)
Elliptical body with 100 points and semi-axes (0.5,0.1)
   Current position: (0.0,0.0)
   Current angle (rad): 0.0

julia> T = RigidTransform((1.0,1.0),π/4)
Rigid-body transform
  Translation: (1.0,1.0)
  Rotation angle (rad): 0.7853981633974483

julia> T(body)
Elliptical body with 100 points and semi-axes (0.5,0.1)
   Current position: (1.0,1.0)
   Current angle (rad): 0.7853981633974483
```
"""
struct RigidTransform
   α   :: Float64
   rot :: Matrix{Float64}
   trans  :: Tuple{Float64,Float64}
end

function RigidTransform(x::Tuple{Float64,Float64},α::Float64)
    rot = [cos(α) -sin(α)
           sin(α) cos(α)]
    RigidTransform(α,rot,x)
end

RigidTransform(x::Union{AbstractVector{T},NTuple{3,T}}) where {T<:Real} = RigidTransform((x[1],x[2]),x[3])

RigidTransform(c::ComplexF64,α::Float64) = RigidTransform((real(c),imag(c)),α)


function Base.show(io::IO, T::RigidTransform)
    name = "Rigid-body transform"
    println(io, name)
    println(io, "  Translation: ($(T.trans[1]),$(T.trans[2]))")
    println(io, "  Rotation angle (rad): $(T.α)")
end

"""
    vec(T::RigidTransform) -> Vector{Float64}

Returns a length-3 vector of the form [x,y,α] corresponding to the translation
and rotation specified by the given transform `T`.
"""
vec(T::RigidTransform) = [T.trans[1],T.trans[2],T.α]



function (T::RigidTransform)(x̃::Float64,ỹ::Float64)
    Xr = T.rot*[x̃,ỹ]
    return T.trans .+ (Xr[1],Xr[2])
end
function (T::RigidTransform)(x̃::AbstractVector{Float64},ỹ::AbstractVector{Float64})
    x = deepcopy(x̃)
    y = deepcopy(ỹ)
    for i = 1:length(x̃)
        x[i],y[i] = T(x̃[i],ỹ[i])
    end
    return x, y
end

function (T::RigidTransform)(b::Body{N}) where {N}
  b.x, b.y = T(b.x̃,b.ỹ)
  b.α = T.α
  b.cent = T.trans
  return b
end


# Things to do here:
#  - add an inverse operation
#  - should be composite operations
#  - probably should distinguish in-place and non-in-place versions
