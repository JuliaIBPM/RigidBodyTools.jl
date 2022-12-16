#=
Rigid body motions
=#


"""
    RigidBodyMotion <: AbstractMotion

A type to store the body's current kinematics

# Fields

- `data`: current kinematic data
- `kin`: a [`Kinematics`](@ref) structure

The first six fields are meant as a cache of the current kinematics
while the `kin` field can be used to find the plate kinematics at any time.
"""
mutable struct RigidBodyMotion <: AbstractMotion
    data::KinematicData
    kin::Kinematics
end

"""
    RigidBodyMotion(Up::Union{ComplexF64,Tuple},Ω[;pivot::Union{ComplexF64,Tuple}=(0.0,0.0)])

Create an instance of constant rigid-body motion with translational velocity Up
and angular velocity Ω with respect to a specified point P.
This point is established by its initial position `pivot`. By default, this
initial position is (0,0). `Up` and `pivot` can be specified by either Tuple or by complex value.
"""
RigidBodyMotion(Up::Union{Number,Tuple}, Ω::Number;kwargs...) = (kin = Constant(Up, Ω; kwargs...); RigidBodyMotion(kin(0), kin))

"""
    RigidBodyMotion(kin::Kinematics)

Create an instance of rigid-body motion with kinematics `kin`.
"""
RigidBodyMotion(kin::Kinematics) = RigidBodyMotion(kin(0), kin)
(m::RigidBodyMotion)(t) = m.kin(t)


function (m::RigidBodyMotion)(t,x̃::Union{Number,Tuple})
  # This expects coordinates in a comoving coordinate system
  # based at c.
  z̃ = complex(x̃...)
  #m.c, m.ċ, m.c̈, m.α, m.α̇, m.α̈ = m.kin(t)
  k = m.kin(t)
  c = complex_translational_position(k)
  ċ = complex_translational_velocity(k)
  c̈ = complex_translational_acceleration(k)
  α = angular_position(k)
  α̇ = angular_velocity(k)
  α̈ = angular_acceleration(k)
  z = exp(im*α)*z̃
  return KinematicData(t, c + z, ċ + im*α̇*z, c̈ + (im*α̈-α̇^2)*z, α, α̇, α̈)
end


"""
    motion_velocity(b::Body,m::RigidBodyMotion,t::Real)

Return the velocity components (as a vector) of a `RigidBodyMotion`
at the given time `t`.
"""
function motion_velocity(b::Body,motion::RigidBodyMotion,t::Real)
   k = motion(t)
  #_,ċ,_,_,α̇,_ = motion(t)
  return [translational_velocity(k)...,angular_velocity(k)]
end

"""
    motion_state(b::Body,m::RigidBodyMotion)

Return the current state vector of body `b` associated with
rigid body motion `m`. It returns the current coordinates
of the body centroid and the angle of the body.
"""
function motion_state(b::Body,m::RigidBodyMotion)
    return [b.cent...,b.α]
end

"""
    update_body!(b::Body,x::AbstractVector,m::RigidBodyMotion)

Update body `b` with the rigid-body motion state vector `x`. The information
in `m` is used for checking lengths only.
"""
function update_body!(b::Body,x::AbstractVector,m::RigidBodyMotion)
    length(x) == length(motion_state(b,m)) || error("wrong length for motion state vector")
    T = RigidTransform(x)
    T(b)
    return b
end


"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 body::Body,motion::AbstractMotion,t::Real[;inertial=true])

Assign the components of body velocity `u` and `v` (in inertial coordinate system)
at surface positions described by inertial coordinates in body `body` at time `t`,
based on supplied motions in the `motion` for the body. If, instead,
`inertial=false`, then it is assumed that `u`, `v` and the surface positions in `body` are all
expressed in comoving coordinates (but velocities are relative to inertial
  frame).
"""
surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,m::RigidBodyMotion,t::Real;kwargs...) =
                 surface_velocity!(u,v,b.x,b.y,b.cent...,b.α,m,t;kwargs...)

"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     x::AbstractVector{Float64},y::AbstractVector{Float64},
                     xc::Real,yc::Real,α::Real,
                     motion::RigidBodyMotion,t::Real[;inertial=true])

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at surface positions described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body. The current position and orientation
of the rigid-body coordinate system are supplied as `xc`, `yc`, and `α`. If, instead,
`inertial=false`, then it is assumed that `u`, `v`, `x`, `y`, `xc`, `yc` are all
expressed in comoving coordinates (but velocities are relative to inertial
  frame).
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                          x::AbstractVector{Float64},y::AbstractVector{Float64},
                          xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real;inertial=true)

  length(u) == length(v) == length(x) == length(y) || error("Inconsistent lengths of vectors")

  k = m(t)
  ċ = complex_translational_velocity(k;inertial=inertial)
  α̇ = angular_velocity(k)
  uc = ċ .+ im*α̇*((x .- xc) .+ im*(y .- yc))
  u .= real.(uc)
  v .= imag.(uc)

  return u, v
end

"""
    surface_velocity(x::AbstractVector{Float64},y::AbstractVector{Float64},
                    xc::Real,yc::Real,α::Real,motion::RigidBodyMotion,t::Real[;inertial=true])

Return the components of rigid body velocities (in inertial components) at surface positions
described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body. If, instead,
`inertial=false`, then it is assumed that `u`, `v`, `x`, `y`, `xc`, `yc` are all
expressed in comoving coordinates (but velocities are relative to inertial
  frame).
"""
surface_velocity(x::AbstractVector{Float64},y::AbstractVector{Float64},
                xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real;kwargs...) =
                surface_velocity!(similar(x),similar(y),x,y,xc,yc,α,m,t;kwargs...)


function show(io::IO, m::RigidBodyMotion)
    println(io, "Rigid Body Motion:")
    print(io, "  $(m.kin)")
end
