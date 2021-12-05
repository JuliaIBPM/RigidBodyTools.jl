




"""
    RigidBodyMotion <: AbstractMotion

A type to store the body's current kinematics

# Fields

- `c`: current centroid position (relative to initial position)
- `ċ`: current centroid velocity
- `c̈`: current centroid acceleration
- `α`: current angle (relative to initial angle)
- `α̇`: current angular velocity
- `α̈`: current angular acceleration
- `kin`: a [`Kinematics`](@ref) structure

The first six fields are meant as a cache of the current kinematics
while the `kin` field can be used to find the plate kinematics at any time.
"""
mutable struct RigidBodyMotion <: AbstractMotion
    c::ComplexF64
    ċ::ComplexF64
    c̈::ComplexF64
    α::Float64
    α̇::Float64
    α̈::Float64

    kin::Kinematics
end

RigidBodyMotion(ċ, α̇) = RigidBodyMotion(0.0im, complex(ċ), 0.0im, 0.0, float(α̇),
                                          0.0, Constant(ċ, α̇))
RigidBodyMotion(kin::Kinematics) = RigidBodyMotion(kin(0)..., kin)
(m::RigidBodyMotion)(t) = m.kin(t)


function (m::RigidBodyMotion)(t,x̃::Tuple{Real,Real})
  # This expects coordinates in body's own coordinate system
  #
  z̃ = ComplexF64(x̃[1],x̃[2])
  m.c, m.ċ, m.c̈, m.α, m.α̇, m.α̈ = m.kin(t)
  z = exp(im*m.α)*z̃
  return m.c + z, m.ċ + im*m.α̇*z, m.c̈ + (im*m.α̈-m.α̇^2)*z
end

"""
    motion_velocity(m::RigidBodyMotion,t::Real)

Return the velocity components (as a vector) of a `RigidBodyMotion`
at the given time `t`.
"""
function motion_velocity(motion::RigidBodyMotion,t::Real)
  _,ċ,_,_,α̇,_ = motion(t)
  return [real(ċ),imag(ċ),α̇]
end

"""
    assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     x::AbstractVector{Float64},y::AbstractVector{Float64},
                     xc::Real,yc::Real,α::Real,
                     motion::RigidBodyMotion,t::Real)

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at positions described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
function assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                          x::AbstractVector{Float64},y::AbstractVector{Float64},
                          xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real)

  length(u) == length(v) == length(x) == length(y) || error("Inconsistent lengths of vectors")

  _,ċ,_,_,α̇,_ = m(t)
  uc = ċ .+ im*α̇*((x .- xc) .+ im*(y .- yc))
  u .= real.(uc)
  v .= imag.(uc)
  #=
  for i = 1:length(x)
      Δz = (x[i]-xc)+im*(y[i]-yc)
      ċi = ċ + im*α̇*Δz
      u[i] = real(ċi)
      v[i] = imag(ċi)
  end
  =#
  return u, v
end

"""
    assign_velocity(x::AbstractVector{Float64},y::AbstractVector{Float64},
                    xc::Real,yc::Real,α::Real,motion::RigidBodyMotion,t::Real)

Return the components of rigid body velocities (in inertial components) at positions
described by coordinates `x`, `y` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
assign_velocity(x::AbstractVector{Float64},y::AbstractVector{Float64},
                xc::Real,yc::Real,α::Real,m::RigidBodyMotion,t::Real) =
                assign_velocity!(similar(x),similar(y),x,y,xc,yc,α,m,t)


function show(io::IO, m::RigidBodyMotion)
    println(io, "Rigid Body Motion:")
    println(io, "  ċ = $(round(m.ċ, digits=2))")
    println(io, "  c̈ = $(round(m.c̈, digits=2))")
    println(io, "  α̇ = $(round(m.α̇, digits=2))")
    println(io, "  α̈ = $(round(m.α̈, digits=2))")
    print(io, "  $(m.kin)")
end
