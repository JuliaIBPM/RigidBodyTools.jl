# Rigid-body transformation routines

export RigidTransform, motion_transform_matrix, motion_transform_matrix_2d,
          force_transform_matrix, force_transform_matrix_2d,
          rotation_about_x, rotation_about_y, rotation_about_z, rotation_from_quaternion,
          quaternion, rotation_about_axis

"""
    RigidTransform(x::Tuple{Real,Real},α::Real)

Construct a rigid-body transform operator, with rotation by angle `α` and
translation specified by `x`. The translation coordinates are specified in the
target coordinate system.

The resulting transform can be used as an operator on pairs of coordinate vectors,
`x` and `y`, or on bodies. For transformation of bodies, it only overwrites the
`x` and `y` fields of the body, but leaves the `x̃` and `ỹ` (body coordinates) intact.

The translation can be provided as either a tuple `(x,y)` or as a complex number.

# Constructors

- `RigidTransform((x,y),α)`
- `RigidTransform(u::Vector{Real})`
- `RigidTransform(u::NTuple{3,Real})`
- `RigidTransform.(u)` where `u` is a collection of vectors or tuples.

# Example

```jldoctest
julia> body = RigidTransform.Ellipse(0.5,0.1,100)
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

function RigidTransform(x::Tuple{Real,Real},α::Real)
    rot = [cos(α) -sin(α)
           sin(α) cos(α)]
    RigidTransform(α,rot,x)
end

RigidTransform(x::Union{AbstractVector{T},NTuple{3,T}}) where {T<:Real} = RigidTransform((x[1],x[2]),x[3])

RigidTransform(c::Complex{T},α::Real) where T<:Real = RigidTransform((real(c),imag(c)),α)


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



function (T::RigidTransform)(x̃::Real,ỹ::Real)
    Xr = T.rot*[x̃,ỹ]
    return T.trans .+ (Xr[1],Xr[2])
end
function (T::RigidTransform)(x̃::AbstractVector{S},ỹ::AbstractVector{S}) where {S<:Real}
    x = deepcopy(x̃)
    y = deepcopy(ỹ)
    for i = 1:length(x̃)
        x[i],y[i] = T(x̃[i],ỹ[i])
    end
    return x, y
end

function (T::RigidTransform)(b::Body{N,C}) where {N,C}
  b.xend, b.yend = T(b.x̃end,b.ỹend)
  b.x, b.y = _midpoints(b.xend,b.yend,C)

  b.α = T.α
  b.cent = T.trans
  return b
end


### Plucker transform matrices ###

"""
    motion_transform_matrix(xA_to_B::SVector,RA_to_B::SMatrix)

Computes the 6 x 6 Plucker transform matrix for motion vectors, transforming
from system A to system B. The input `xA_to_B` is the Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `RA_to_B` is the rotation
matrix transforming coordinates in system A to those in system B. The resulting matrix has the form

``{}^B T^{(m)}_A = \\begin{bmatrix} R & 0 \\\\ 0 & R \\end{bmatrix} \\begin{bmatrix} 1 & 0 \\\\ -x^\\times & 1 \\end{bmatrix}``
"""
function motion_transform_matrix(xA_to_B::SVector{3},RA_to_B::SMatrix{3,3})
    xM = cross_matrix(xA_to_B)
    id = SMatrix{3,3}(I)
    O = SMatrix{3,3}(zeros(9))
    Mtrans = SMatrix{6,6}([id O; -xM id])
    Mrot = SMatrix{6,6}([RA_to_B O; O RA_to_B])
    return Mrot*Mtrans
end

"""
    force_transform_matrix(xA_to_B::SVector,RA_to_B::SMatrix)

Computes the 6 x 6 Plucker transform matrix for force vectors, transforming
from system A to system B. The input `xA_to_B` is the Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `RA_to_B` is the rotation
matrix transforming coordinates in system A to those in system B. The resulting matrix has the form

``{}^B T^{(f)}_A = \\begin{bmatrix} R & 0 \\\\ 0 & R \\end{bmatrix} \\begin{bmatrix} 1 & -x^\\times \\\\ 0 & 1 \\end{bmatrix}``
"""
function force_transform_matrix(xA_to_B::SVector{3},RA_to_B::SMatrix{3,3})
    xM = cross_matrix(xA_to_B)
    id = SMatrix{3,3}(I)
    O = SMatrix{3,3}(zeros(9))
    Mtrans = SMatrix{6,6}([id -xM; O id])
    Mrot = SMatrix{6,6}([RA_to_B O; O RA_to_B])
    return Mrot*Mtrans
end

for f in [:motion, :force]

    fname = Symbol(f,"_transform_matrix")
    fname2d = Symbol(fname,"_2d")

    @eval $fname(x::AbstractVector,R::AbstractMatrix) = $fname(SVector{3}(x),SMatrix{3,3}(R))

    @eval function $fname2d(xA_to_B::SVector{2},θ::Real)
        x_3d = SVector{3}([xA_to_B... 0.0])

        # We need to rotate in -Θ direction to get the A-to-B rotation matrix
        R_3d = rotation_about_z(-θ)
        T3d = $fname(x_3d,R_3d)
        return SMatrix{3,3}(T3d[3:5,3:5])
    end

    @eval $fname2d(x::Vector,θ::Real) = $fname2d(SVector{2}(x),θ)

    @eval $fname2d(x::Tuple,θ::Real) = $fname2d(SVector{2}(x...),θ)

    @eval $fname2d(x::Real,y::Real,θ::Real) = $fname2d(SVector{2}([x,y]),θ)

    @eval $fname2d(v::Vector) = $fname2d(v...)

    @eval $fname2d(T::RigidTransform) = $fname2d(vec(T))


end
"""
    motion_transform_matrix_2d(xA_to_B,θ::Real) -> SMatrix{3,3}

Computes the 3 x 3 2D Plucker transform matrix for motion vectors, transforming
from system A to system B. The input `xA_to_B` is the 2-d Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `θ` is the angle of system B relative
to system A. `xA_to_B` can be in the form
of a static vector, a vector, or a tuple.
""" motion_transform_matrix_2d

"""
    motion_transform_matrix_2d(T::RigidTransform) -> SMatrix{3,3}

Computes the 3 x 3 2D Plucker transform matrix for motion vectors, transforming
from system A to system B, from the rigid transform `T`.
""" motion_transform_matrix_2d(::RigidTransform)


"""
    force_transform_matrix_2d(xA_to_B,θ::Real) -> SMatrix{3,3}

Computes the 3 x 3 2D Plucker transform matrix for force vectors, transforming
from system A to system B. The input `xA_to_B` is the 2-d Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `θ` is the angle of system B relative
to system A. `xA_to_B` can be in the form
of a static vector, a vector, or a tuple.
""" force_transform_matrix_2d

"""
    force_transform_matrix_2d(T::RigidTransform) -> SMatrix{3,3}

Computes the 3 x 3 2D Plucker transform matrix for force vectors, transforming
from system A to system B, from the rigid transform `T`.
""" force_transform_matrix_2d(::RigidTransform)



### Rotation matrices ###

"""
    rotation_about_z(θ::Real) -> SMatrix{3,3}

Constructs the rotation matrix corresponding to rotation about the z axis
by angle θ.  The matrix transforms the coordinates of a vector
in system A to coordinates in system B, where Θ is the angle of A with respect to B.
(Flip the sign of the input Θ if you want it to correspond to B's angle relative to A.)
"""
function rotation_about_z(θ::Real)
    cth, sth = cos(θ), sin(θ)
    R = @SMatrix[cth -sth 0.0; sth cth 0.0; 0.0 0.0 1.0]
    return R
end

"""
    rotation_about_y(θ::Real) -> SMatrix{3,3}

Constructs the rotation matrix corresponding to rotation about the y axis
by angle θ.  The matrix transforms the coordinates of a vector
in system A to coordinates in system B, where Θ is the angle of A with respect to B.
(Flip the sign of the input Θ if you want it to correspond to B's angle relative to A.)
"""
function rotation_about_y(θ::Real)
    cth, sth = cos(θ), sin(θ)
    R = @SMatrix[cth 0.0 sth; 0.0 1.0 0.0; -sth 0.0 cth]
    return R
end

"""
    rotation_about_x(θ::Real) -> SMatrix{3,3}

Constructs the rotation matrix corresponding to rotation about the x axis
by angle θ. The matrix transforms the coordinates of a vector
in system A to coordinates in system B, where Θ is the angle of A with respect to B.
(Flip the sign of the input Θ if you want it to correspond to B's angle relative to A.)
"""
function rotation_about_x(θ::Real)
    cth, sth = cos(θ), sin(θ)
    R = @SMatrix[1.0 0.0 0.0; 0.0 cth -sth; 0.0 sth cth]
    return R
end

"""
    rotation_about_axis(θ::Real,v::Vector) -> SMatrix{3,3}

Constructs the rotation matrix corresponding to rotation about the axis `v`
by angle θ. The matrix transforms the coordinates of a vector
in system A to coordinates in system B, where Θ is the angle of A with respect to B.
(Flip the sign of the input Θ if you want it to correspond to B's angle relative to A.)
"""
rotation_about_axis(ϴ::Real,v::SVector{3}) = rotation_from_quaternion(quaternion(ϴ,v))

rotation_about_axis(Θ::Real,v::Vector) = rotation_about_axis(Θ,SVector{3}(v))

"""
    rotation_from_quaternion(q::SVector{4}) -> SMatrix{3,3}

Constructs a rotation matrix from a given quaternion `q` (a 4-dimensional vector).
The matrix transforms the coordinates of a vector in system A to coordinates in system B,
where A is rotated with respect to B
by θ.
"""
function rotation_from_quaternion(q::SVector{4})
    w, x, y, z = q
    R = @SMatrix[1-2y^2-2z^2 2*x*y-2*z*w     2*x*z+2*y*w;
                 2*x*y+2*z*w     1-2x^2-2z^2 2*y*z-2*x*w;
                 2*x*z-2*y*w     2*y*z+2*x*w     1-2x^2-2y^2]
    return R
end

rotation_from_quaternion(q::Vector) = rotation_from_quaternion(SVector{4}(q))

"""
    quaternion(θ::Real,v::SVector{3}) -> SVector{4}

Form the quaternion ``q = w + x\\hat{i} + y\\hat{j} + z\\hat{k}``
from the given rotation angle `θ` and rotation axis `v`. Note that `v` need not be a unit vector.
"""
function quaternion(θ::Real,v::SVector{3})
    vnorm = v/norm(v)
    q = SVector{4}([cos(θ/2),sin(θ/2)*vnorm...])
end

quaternion(θ::Real,v::Vector) = quaternion(θ,SVector{3}(v))



"""
    cross_matrix(v::SVector) -> SMatrix

Takes a vector `v` and forms the correspsonding cross-product matrix `M`, so that
``v \\times w`` is equivalent to ``M\\cdot w``.
"""
function cross_matrix(v::SVector{3})
    M = @SMatrix[0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    return M
end

cross_matrix(v::Vector) = cross_matrix(SVector{3}(v))



# Things to do here:
#  - add an inverse operation
#  - should be composite operations
#  - probably should distinguish in-place and non-in-place versions
