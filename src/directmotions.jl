export DirectlySpecifiedMotion, BasicDirectMotion

abstract type DirectlySpecifiedMotion <: AbstractMotion end

#=
To create a subtype of DirectlySpecifiedMotion, one must
extend `motion_velocity(m,t)` to convert the instantaneous
velocity components into a vector (e.g., to be used for
advancing the associated points in a time-marching scheme),
`motion_state(b,m)` to convert the instantaneous coordinates
of body `b` into a state vector,
and `surface_velocity!(u,v,body,m,t)`, to supply the
surface values of the motion in-place in vectors `u` and `v`,
which must be of the same length as `b.x` and `b.y`.
=#

struct BasicDirectMotion{VT} <: DirectlySpecifiedMotion
    u :: VT
    v :: VT
end


"""
    motion_velocity(m::BasicDirectMotion,t::Real)

Return the velocity components (as a vector) of a `BasicDirectMotion`
at the given time `t`.
"""
motion_velocity(m::BasicDirectMotion,t::Real) = vcat(m.u,m.v)


"""
    motion_state(b::Body,m::BasicDirectMotion)

Return the current state vector of body `b` associated with
direct motion `m`. It returns the concatenated coordinates
of the body surface (in the inertial coordinate system).
"""
function motion_state(b::Body,m::BasicDirectMotion)
    return vcat(b.x,b.y)
end


"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     b::Body,motion::BasicDirectMotion,t::Real)

Assign the components of velocity `u` and `v` (in inertial coordinate system)
at surface positions described by points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                           b::Body,m::DirectlySpecifiedMotion,t::Real)
     u .= m.u
     v .= m.v
     return u, v
end


"""
    surface_velocity(b::Body,motion::BasicDirectMotion,t::Real)

Return the components of velocities (in inertial components) at surface positions
described by points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
surface_velocity(b::Body,m::BasicDirectMotion,t::Real) =
                surface_velocity!(similar(b.x),similar(b.y),b,m,t)
