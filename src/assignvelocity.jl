"""
    surface_velocity(body::Body,motion::AbstractMotion,t::Real[;inertial=true])

Return the components of rigid body velocity (in inertial coordinate system)
at surface positions described by inertial coordinates in body `body` at time `t`,
based on supplied motions in `motion` for the body. If `inertial=false`,
then velocities are computed and body positions are assumed to be in comoving coordinates.
"""
surface_velocity(b::Body,a...;kwargs...) = surface_velocity!(zero(b.x),zero(b.y),b,a...;kwargs...)
