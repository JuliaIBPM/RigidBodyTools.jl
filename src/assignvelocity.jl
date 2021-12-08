"""
    surface_velocity(body::Body,motion::AbstractMotion,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in `motion` for the body.

As a shorthand, you can also apply this as `motion(t,body)`.
"""
surface_velocity(b::Body,a...) = surface_velocity!(zero(b.x),zero(b.y),b,a...)


#(m::RigidBodyMotion)(t::Real,b::Body) = surface_velocity(b,m,t)

#(m::DirectlySpecifiedMotion)(t::Real,b::Body) = surface_velocity(b,m,t)
