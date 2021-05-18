"""
    assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 body::Body,motion::RigidBodyMotion,t::Real)

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in the RigidBodyMotion `motion` for the body.
"""
assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,m::Union{RigidBodyMotion,DirectlySpecifiedMotion},t::Real) =
                 assign_velocity!(u,v,b.x,b.y,b.cent...,b.α,m,t)

assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,m::DirectlySpecifiedMotion,t::Real) =
                 assign_velocity!(u,v,b.x,b.y,b.cent...,b.α,m,t)


"""
    assign_velocity(body::Body,motion::RigidBodyMotion,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in the RigidBodyMotion `motion` for the body.

As a shorthand, you can also apply this as `motion(t,body)`.
"""
assign_velocity(b::Body,a...) = assign_velocity!(zero(b.x),zero(b.y),b,a...)


(m::RigidBodyMotion)(t::Real,b::Body) = assign_velocity(b,m,t)

(m::DirectlySpecifiedMotion)(t::Real,b::Body) = assign_velocity(b,m,t)


"""
    assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     bl::BodyList,ml::RigidMotionList,t::Real)

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the RigidMotionList `ml` for each body.
"""
function assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 bl::BodyList,ml::Union{RigidMotionList,DirectlySpecifiedMotionList},t::Real)

   for i in 1:length(bl)
      assign_velocity!(view(u,bl,i),view(v,bl,i),bl[i],ml[i],t)
   end
   return u, v
end

"""
    assign_velocity(bl::BodyList,ml::RigidMotionList,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the RigidMotionList `ml` for each body.

As a shorthand, you an also apply this as `ml(t,bl)`.
"""
assign_velocity(bl::BodyList,ml::Union{RigidMotionList,DirectlySpecifiedMotionList},t::Real) =
    assign_velocity!(zeros(Float64,numpts(bl)),zeros(Float64,numpts(bl)),bl,ml,t)

(ml::RigidMotionList)(t::Real,bl::BodyList) = assign_velocity(bl,ml,t)

(ml::DirectlySpecifiedMotionList)(t::Real,bl::BodyList) = assign_velocity(bl,ml,t)
