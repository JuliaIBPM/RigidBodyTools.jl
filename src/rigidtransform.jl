# Rigid-body transformation routines

# Things to do here:
#  - Combine RigidTransform with MotionTransform
#  - probably should distinguish in-place and non-in-place versions


const O3VECTOR = SVector{3}(zeros(Float64,3))
const I3 = SMatrix{3,3}(I)
const O3 = SMatrix{3,3}(zeros(Float64,9))

plucker_dimension(::Val{2}) = 3
plucker_dimension(::Val{3}) = 6

### Plucker vectors ###
abstract type AbstractPluckerVector{ND} end

"""
    PluckerMotion(data::AbstractVector) -> SVector

Creates an instance of a Plucker motion vector,

``v = \\begin{bmatrix} \\Omega \\\\ U \\end{bmatrix}``

using the vector in `data`. If `data` is of length 6, then it creates a 3d motion vector,
and the first 3 entries are assumed to comprise the rotational component `\\Omega` and the last 3 entries the
translational component `U`. If `data` is of length 3, then it creates
a 2d motion vector, assuming that the first entry in `data` represents
the rotational component and the second and third entries the
x and y translational components.
"""
struct PluckerMotion{ND} <: AbstractPluckerVector{ND}
  data :: SVector
end

function show(io::IO, p::PluckerMotion{3})
  print(io, "3d Plucker motion vector, Ω = $(p.data[1:3]), v = $(p.data[4:6])")
end

function show(io::IO, p::PluckerMotion{2})
  print(io, "2d Plucker motion vector, Ω = $(p.data[1]), U = $(p.data[2:3])")
end

"""
    PluckerForce(data::AbstractVector) -> SVector

Creates an instance of a Plucker force vector,

``f = \\begin{bmatrix} M \\\\ F \\end{bmatrix}``

using the vector in `data`. If `data` is of length 6, then it creates a 3d force vector,
and the first 3 entries are assumed to comprise the moment component `M` and the last 3 entries the
force component `F`. If `data` is of length 3, then it creates
a 2d force vector, assuming that the first entry in `data` represents
the moment component and the second and third entries the
x and y force components.
"""
struct PluckerForce{ND} <: AbstractPluckerVector{ND}
  data :: SVector
end

function show(io::IO, p::PluckerForce{3})
  print(io, "3d Plucker force vector, M = $(p.data[1:3]), F = $(p.data[4:6])")
end

function show(io::IO, p::PluckerForce{2})
  print(io, "2d Plucker force vector, M = $(p.data[1]), F = $(p.data[2:3])")
end

for typename in [:PluckerMotion,:PluckerForce]

  fname_underscore = Symbol("_",lowercase(string(typename)))

  @eval function $typename{ND}() where {ND}
    pd = plucker_dimension(Val(ND))
    $typename{ND}(SVector{pd}(zeros(Float64,pd)))
  end

  @eval $typename(v::SVector{3}) = $typename{2}(v)

  @eval $typename(Ω::Real,U::AbstractVector) = $typename([Ω,U...])

  @eval $typename(Ω::AbstractVector,U::AbstractVector) = $typename([Ω...,U...])


  @eval $typename(v::SVector{6}) = $typename{3}(v)

  @eval $typename(v::AbstractVector) = $fname_underscore(v,Val(length(v)))

  @eval $fname_underscore(v,::Val{3}) = $typename(SVector{3}(v))

  @eval $fname_underscore(v,::Val{6}) = $typename(SVector{6}(v))

  @eval (+)(a::$typename{ND},b::$typename{ND}) where {ND} = $typename{ND}(a.data+b.data)

  @eval (-)(a::$typename{ND},b::$typename{ND}) where {ND} = $typename{ND}(a.data-b.data)

  @eval (-)(a::$typename{ND}) where {ND} = $typename{ND}(-a.data)

  @eval @propagate_inbounds getindex(A::$typename, i::Int) = A.data[i]
  @eval @propagate_inbounds getindex(A::$typename, I...) = A.data[I...]

  @eval iterate(A::$typename) = iterate(A.data)
  @eval iterate(A::$typename,I) = iterate(A.data,I)
  @eval size(A::$typename) = size(A.data)
  @eval length(A::$typename) = length(A.data)


end

"""
    dot(f::PluckerForce,v::PluckerMotion) -> Real

Calculate the scalar product between force `f` and motion `v`. The
commutation of this is also possible, `dot(v,f)`.
"""
dot(f::PluckerForce{ND},v::PluckerMotion{ND}) where {ND} = dot(f.data,v.data)

dot(v::PluckerMotion{ND},f::PluckerForce{ND}) where {ND} = dot(f,v)




### Plucker transform matrices ###
abstract type AbstractTransformOperator{ND} end

struct MotionTransform{ND} <: AbstractTransformOperator{ND}
   x :: SVector
   R :: SMatrix
   matrix :: SMatrix
end

struct ForceTransform{ND} <: AbstractTransformOperator{ND}
   x :: SVector
   R :: SMatrix
   matrix :: SMatrix
end

const RigidTransform = MotionTransform{ND} where {ND}

function show(io::IO, p::MotionTransform{3})
  print(io, "3d motion transform, x = $(p.x), R = $(p.R)")
end

function show(io::IO, p::ForceTransform{3})
  print(io, "3d force transform, x = $(p.x), R = $(p.R)")
end

function show(io::IO, p::MotionTransform{2})
  print(io, "2d motion transform, x = $(p.x[1:2]), R = $(p.R[1:2,1:2])")
end

function show(io::IO, p::ForceTransform{2})
  print(io, "2d force transform, x = $(p.x[1:2]), R = $(p.R[1:2,1:2])")
end

translation(T::AbstractTransformOperator{3}) = T.x
rotation(T::AbstractTransformOperator{3}) = T.R

translation(T::AbstractTransformOperator{2}) = SVector{2}(T.x[1:2])
rotation(T::AbstractTransformOperator{2}) = SMatrix{2,2}(T.R[1:2,1:2])


(*)(T::AbstractTransformOperator,v) = T.matrix*v
(*)(v,T::AbstractTransformOperator) = v*T.matrix

(*)(T::MotionTransform,v::PluckerMotion) = PluckerMotion(T*v.data)
(*)(T::ForceTransform,v::PluckerForce) = PluckerForce(T*v.data)


function transpose(T::MotionTransform{ND}) where {ND}
    x = cross_vector(T.R*cross_matrix(-T.x)*T.R')
    ForceTransform{ND}(x,transpose(T.R),transpose(T.matrix))
end
function transpose(T::ForceTransform{ND}) where {ND}
    x = cross_vector(T.R*cross_matrix(-T.x)*T.R')
    MotionTransform{ND}(x,transpose(T.R),transpose(T.matrix))
end

"""
    inv(X::AbstractTransformOperator) -> AbstractTransformOperator

Return the inverse of the motion or force transform `X`.
"""
inv(T::MotionTransform{ND}) where {ND} = transpose(ForceTransform{ND}(T.x,T.R))
inv(T::ForceTransform{ND}) where {ND} = transpose(MotionTransform{ND}(T.x,T.R))

rotation_transform(T::MotionTransform{ND}) where {ND} = MotionTransform{ND}(O3VECTOR,T.R)
rotation_transform(T::ForceTransform{ND}) where {ND} = ForceTransform{ND}(O3VECTOR,T.R)


vec(T::AbstractTransformOperator{2}) = [T.x[1],T.x[2],_get_angle_of_2d_transform(T)]


function _get_angle_of_2d_transform(T::AbstractTransformOperator{2})
  eith = T.R[1,1] + im*T.R[1,2]
  return real(-im*log(eith))
end


### MotionTransform acting on position coordinate data ###

function (T::MotionTransform{2})(x̃::Real,ỹ::Real)
    Xr = T.R'*[x̃,ỹ,0.0]
    return T.x[1] + Xr[1], T.x[2] + Xr[2]
end
function (T::MotionTransform{2})(x̃::AbstractVector{S},ỹ::AbstractVector{S}) where {S<:Real}
    x = deepcopy(x̃)
    y = deepcopy(ỹ)
    for i = 1:length(x̃)
        x[i],y[i] = T(x̃[i],ỹ[i])
    end
    return x, y
end

"""
    (T::MotionTransform)(b::Body) -> Body

Transforms a body `b` using the given `MotionTransform`, creating a copy of this body
with the new configuration. In using this
transform `T` (which defines a transform from system A to system B), A is interpreted as an inertial coordinate
system and B as the body system. Thus, the position vector in `T` is interpreted
as the relative position of the body in inertial coordinates and the inverse of the rotation
operator is applied to transform body-fixed coordinates to the inertial frame.
"""
function (T::MotionTransform{2})(b::Body)
    b_transform = deepcopy(b)
    update_body!(b_transform,T)
end

"""
    update_body!(b::Body,t::MotionTransform)

Transforms a body (in-place) using the given `MotionTransform`. In using this
transform `T` (which defines a transform from system A to system B), A is interpreted as an inertial coordinate
system and B as the body system. Thus, the position vector in `T` is interpreted
as the relative position of the body in inertial coordinates and the inverse of the rotation
operator is applied to transform body-fixed coordinates to the inertial frame.
"""
function update_body!(b::Body{N,C},T::MotionTransform{2}) where {N,C}
  b.xend, b.yend = T(b.x̃end,b.ỹend)
  b.x, b.y = _midpoints(b.xend,b.yend,C)

  b.α = _get_angle_of_2d_transform(T)
  b.cent = (T.x[1], T.x[2])
  return b
end


###


function _motion_transform_matrix(xA_to_B::SVector{3},RA_to_B::SMatrix{3,3})
    xM = cross_matrix(xA_to_B)
    Mtrans = SMatrix{6,6}([I3 O3;
                          -xM I3])
    Mrot = SMatrix{6,6}([RA_to_B O3;
                         O3 RA_to_B])
    return Mrot*Mtrans
end

function _force_transform_matrix(xA_to_B::SVector{3},RA_to_B::SMatrix{3,3})
    xM = cross_matrix(xA_to_B)
    Mtrans = SMatrix{6,6}([I3 -xM;
                           O3 I3])
    Mrot = SMatrix{6,6}([RA_to_B O3;
                          O3   RA_to_B])
    return Mrot*Mtrans
end

for f in [:motion, :force]

    fname_underscore = Symbol("_",f,"_transform_matrix")
    fname2d_underscore = Symbol(fname_underscore,"_2d")
    typename = Symbol(uppercasefirst(string(f)),"Transform")


    @eval $typename{3}(x_3d::SVector{3},R_3d::SMatrix{3,3}) = $typename(x_3d,R_3d)

    @eval $typename{2}(x_3d::SVector{3},R_3d::SMatrix{3,3}) = $typename(SVector{2}(x_3d[1:2]),R_3d)

    @eval $typename{ND}() where {ND} = $typename{ND}(O3VECTOR,rotation_identity())


    @eval function $typename(x_3d::SVector{3},R_3d::SMatrix{3,3})
        M = $fname_underscore(x_3d,R_3d)
        return $typename{3}(x_3d,R_3d,M)
    end

    @eval function $typename(x_2d::SVector{2},R_3d::SMatrix{3,3})
        x_3d = SVector{3}([x_2d... 0.0])
        M = $fname2d_underscore(x_3d,R_3d)
        return $typename{2}(x_3d,R_3d,M)
    end

    @eval function $typename(x_2d::SVector{2},R_2d::SMatrix{2,2})
      R_3d = SMatrix{3,3}([R_2d zeros(Float64,2,1); zeros(Float64,1,2) 1.0])
      $typename(x_2d,R_3d)
    end

    @eval function $typename(x::AbstractVector,R::AbstractMatrix)
        lenx = length(x)
        nx, ny = size(R)
        @assert (lenx == 2 || lenx == 3) "x has inconsistent length"
        @assert (nx == ny == 3 || nx == ny == 2) "Rotation matrix has inconsistent dimensions"

        $typename(SVector{lenx}(x),SMatrix{nx,ny}(R))
    end

    @eval $typename(x::SVector{2},Θ::Real) = $typename(x,rotation_about_z(-Θ))

    @eval $typename(x::Vector,θ::Real) = $typename(SVector{2}(x),θ)

    @eval $typename(x::Tuple,θ::Real) = $typename(SVector{2}(x...),θ)

    @eval $typename(x::Real,y::Real,θ::Real) = $typename(SVector{2}([x,y]),θ)

    @eval $typename(c::Complex{T},θ::Real) where T<:Real = $typename((real(c),imag(c)),θ)

    @eval $typename(v::AbstractVector) = $typename(v...)

    @eval $typename(T::RigidTransform) = $typename(vec(T))

    @eval function $fname2d_underscore(x_3d::SVector{3},R_3d::SMatrix{3,3})
        M_3d = $fname_underscore(x_3d,R_3d)
        return SMatrix{3,3}(M_3d[3:5,3:5])
    end

    @eval function (*)(T1::$typename{ND},T2::$typename{ND}) where {ND}
      x12 = cross_vector(T2.R'*cross_matrix(T1.x)*T2.R + cross_matrix(T2.x))
      R12 = T1.R*T2.R
      $typename{ND}(x12,R12,T1.matrix*T2.matrix)
    end

end


"""
    MotionTransform(xA_to_B::SVector,RA_to_B::SMatrix) -> MotionTransform

Computes the Plucker transform matrix for motion vectors, transforming
from system A to system B. The input `xA_to_B` is the Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `RA_to_B` is the rotation
matrix transforming coordinates in system A to those in system B. The resulting matrix has the form

``{}^B T^{(m)}_A = \\begin{bmatrix} R & 0 \\\\ 0 & R \\end{bmatrix} \\begin{bmatrix} 1 & 0 \\\\ -x^\\times & 1 \\end{bmatrix}``

One can also provide `xA_to_B` as a standard vector and `RA_to_B` as a standard
3 x 3 matrix.

If `xA_to_B` has length 3, then a three-dimensional transform (a 6 x 6 Plucker transform)
is created. If `xA_to_B` has length 2, then a two-dimensional transform (3 x 3 Plucker transform)
is returned.
""" MotionTransform(::SVector{3},::SMatrix)


"""
    MotionTransform(xA_to_B,θ::Real) -> MotionTransform

Computes the 3 x 3 2D Plucker transform matrix for motion vectors, transforming
from system A to system B. The input `xA_to_B` is the 2-d Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `θ` is the angle of system B relative
to system A. `xA_to_B` can be in the form
of a static vector, a vector, or a tuple.
""" MotionTransform(::AbstractVector,::Real)

"""
    MotionTransform(T::RigidTransform) -> MotionTransform

Computes the 3 x 3 2D Plucker transform matrix for motion vectors, transforming
from system A to system B, from the rigid transform `T`.
""" MotionTransform(::RigidTransform)


"""
    ForceTransform(xA_to_B::SVector,RA_to_B::SMatrix) -> ForceTransform

Computes the 6 x 6 Plucker transform matrix for force vectors, transforming
from system A to system B. The input `xA_to_B` is the Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `RA_to_B` is the rotation
matrix transforming coordinates in system A to those in system B. The resulting matrix has the form

``{}^B T^{(f)}_A = \\begin{bmatrix} R & 0 \\\\ 0 & R \\end{bmatrix} \\begin{bmatrix} 1 & -x^\\times \\\\ 0 & 1 \\end{bmatrix}``
""" ForceTransform(::SVector{3},::SMatrix)

"""
    ForceTransform(xA_to_B,θ::Real) -> ForceTransform

Computes the 3 x 3 2D Plucker transform matrix for force vectors, transforming
from system A to system B. The input `xA_to_B` is the 2-d Euclidean vector from the origin of A to
the origin of B, expressed in A coordinates, and `θ` is the angle of system B relative
to system A. `xA_to_B` can be in the form
of a static vector, a vector, or a tuple.
""" ForceTransform(::AbstractVector,::Real)

"""
    ForceTransform(T::RigidTransform) -> ForceTransform

Computes the 3 x 3 2D Plucker transform matrix for force vectors, transforming
from system A to system B, from the rigid transform `T`.
""" ForceTransform(::RigidTransform)



### Rotation matrices ###

rotation_identity() = I3

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

quaternion(q::AbstractVector) = quaternion(q[1],SVector{3}(q[2:4]))



"""
    cross_matrix(v::SVector) -> SMatrix

Takes a vector `v` and forms the corresponding cross-product matrix `M`, so that
``v \\times w`` is equivalent to ``M\\cdot w``.
"""
function cross_matrix(v::SVector{3})
    M = @SMatrix [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    return M
end

cross_matrix(v::Vector) = cross_matrix(SVector{3}(v))

"""
    cross_vector(M::SMatrix) -> SVector

Takes a matrix `M` and forms the corresponding axial vector `v` from the matrix's
anti-symmetric part (``M_a = (M-M^{T})/2``), so that ``v \\times w`` is equivalent
to ``M_a\\cdot w``.
"""
function cross_vector(M::SMatrix{3,3})
    Ma = 0.5*(M - transpose(M))
    return @SVector [Ma[3,2],Ma[1,3],Ma[2,1]]
end

cross_vector(M::Matrix) = cross_vector(SMatrix{3,3}(M))
