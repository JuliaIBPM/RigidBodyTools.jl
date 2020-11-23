"""
    assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 body::Body,T::RigidTransform,motion::RigidBodyMotion,t::Real)

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in the RigidBodyMotion `motion` for the body.
"""
assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,T::RigidTransform,m::RigidBodyMotion,t::Real) =
                assign_velocity!(u,v,b.x,b.y,vec(T)...,m,t)

"""
    assign_velocity(body::Body,T::RigidTransform,motion::RigidBodyMotion,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at positions described by coordinates inertial coordinates in body in `body` at time `t`,
based on supplied motions in the RigidBodyMotion `motion` for the body.
"""
assign_velocity(b::Body,a...) = assign_velocity!(zero(b.x),zero(b.y),b,a...)


"""
    assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     bl::BodyList,tl::RigidTransformList,ml::RigidMotionList,t::Real)

Assign the components of rigid body velocity `u` and `v` (in inertial coordinate system)
at positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the RigidMotionList `ml` for each body.
"""
function assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 bl::BodyList,tl::RigidTransformList,ml::RigidMotionList,t::Real)

   for i in 1:length(bl)
      assign_velocity!(view(u,bl,i),view(v,bl,i),bl[i],tl[i],ml[i],t)
   end
   return u, v
end

"""
    assign_velocity(bl::BodyList,tl::RigidTransformList,ml::RigidMotionList,t::Real)

Return the components of rigid body velocity (in inertial coordinate system)
at positions described by coordinates inertial coordinates in each body in `bl` at time `t`,
based on supplied motions in the RigidMotionList `ml` for each body.
"""
assign_velocity(bl::BodyList,tl::RigidTransformList,ml::RigidMotionList,t::Real) =
    assign_velocity!(zeros(Float64,numpts(bl)),zeros(Float64,numpts(bl)),bl,tl,ml,t)
