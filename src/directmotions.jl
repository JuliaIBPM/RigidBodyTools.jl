export DirectlySpecifiedMotion, BasicDirectMotion

abstract type DirectlySpecifiedMotion <: AbstractMotion end

#=
To create a subtype of DirectlySpecifiedMotion, one must
extend `surface_velocity!(u,v,body,m,t)`, to supply the
surface values of the motion in-place in vectors `u` and `v`,
which must be of the same length as `b.x` and `b.y`.
=#

"""
    BasicDirectMotion(u::Vector{Float64},v::Vector{Float64})

Create an instance of basic directly-specified (constant)
velocity, to be associated with a body whose length
is the same as `u` and `v`.
"""
struct BasicDirectMotion{VT} <: DirectlySpecifiedMotion
    u :: VT
    v :: VT
end

"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     b::Body,motion::BasicDirectMotion,t::Real)

Assign the components of velocity `u` and `v` (in inertial coordinate system)
at surface positions described by points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                           b::Body,m::BasicDirectMotion,t::Real)
     u .= m.u
     v .= m.v
     return u, v
end


"""
    motion_velocity(b::Body,m::DirectlySpecifiedMotion,t::Real)

Return the velocity components (as a vector) of a `DirectlySpecifiedMotion`
at the given time `t`.
"""
function motion_velocity(b::Body,m::DirectlySpecifiedMotion,t::Real)
    u, v = zero(b.x), zero(b.y)
    surface_velocity!(u,v,b,m,t)
    return vcat(u,v)
end


"""
    motion_state(b::Body,m::DirectlySpecifiedMotion)

Return the current state vector of body `b` associated with
direct motion `m`. It returns the concatenated coordinates
of the body surface (in the inertial coordinate system).
"""
function motion_state(b::Body,m::DirectlySpecifiedMotion)
    return vcat(b.x,b.y)
end




#=
"""
    surface_velocity(b::Body,motion::DirectlySpecifiedMotion,t::Real)

Return the components of velocities (in inertial components) at surface positions
described by points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
surface_velocity(b::Body,m::DirectlySpecifiedMotion,t::Real) =
                surface_velocity!(similar(b.x),similar(b.y),b,m,t)

=#
