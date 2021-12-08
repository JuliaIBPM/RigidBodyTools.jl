"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 body::Body,motion::AbstractMotion,t::Real)

Assign the components of body velocity `u` and `v` (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in the `motion` for the body.
"""
surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,m::RigidBodyMotion,t::Real) =
                 surface_velocity!(u,v,b.x,b.y,b.cent...,b.α,m,t)




"""
    surface_velocity(body::Body,motion::AbstractMotion,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in the RigidBodyMotion `motion` for the body.

As a shorthand, you can also apply this as `motion(t,body)`.
"""
surface_velocity(b::Body,a...) = surface_velocity!(zero(b.x),zero(b.y),b,a...)


(m::RigidBodyMotion)(t::Real,b::Body) = surface_velocity(b,m,t)

(m::DirectlySpecifiedMotion)(t::Real,b::Body) = surface_velocity(b,m,t)


"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     bl::BodyList,ml::MotionList,t::Real)

Assign the components of velocity `u` and `v` (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the MotionList `ml` for each body.
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 bl::BodyList,ml::MotionList,t::Real)

   for i in 1:length(bl)
      surface_velocity!(view(u,bl,i),view(v,bl,i),bl[i],ml[i],t)
   end
   return u, v
end

"""
    surface_velocity(bl::BodyList,ml::MotionList,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at surface positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the MotionList `ml` for each body.

As a shorthand, you an also apply this as `ml(t,bl)`.
"""
surface_velocity(bl::BodyList,ml::MotionList,t::Real) =
    surface_velocity!(zeros(Float64,numpts(bl)),zeros(Float64,numpts(bl)),bl,ml,t)

(ml::MotionList)(t::Real,bl::BodyList) = surface_velocity(bl,ml,t)
