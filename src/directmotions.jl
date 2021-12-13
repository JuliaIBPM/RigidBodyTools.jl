abstract type AbstractDirectlySpecifiedMotion <: AbstractMotion end

#=
To create a subtype of AbstractDirectlySpecifiedMotion, one must
extend `motion_velocity(body,m,t)`, to supply the
surface values of the motion in-place in vectors `u` and `v`.
These are interpreted as an update of the body-fixed coordinates,
`b.x̃end` and `b.ỹend`, in the body's coordinate system.
=#

"""
    BasicDirectMotion(u::Vector{Float64},v::Vector{Float64})

Create an instance of basic directly-specified (constant)
velocity, to be associated with a body whose length
is the same as `u` and `v`.
"""
struct BasicDirectMotion{VT} <: AbstractDirectlySpecifiedMotion
    u :: VT
    v :: VT
end


"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     b::Body,motion::AbstractDirectlySpecifiedMotion,t::Real)

Assign the components of velocity `u` and `v` (in inertial coordinate system)
at surface positions described by `x` and `y` points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body. This function calls the
function `motion_velocity` associated with the given motion type.
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                           b::Body{N,C},m::AbstractDirectlySpecifiedMotion,t::Real) where {N,C}


     vel = motion_velocity(b,m,t)
     lenx = length(vel)

     # interpolate to the midpoints
     umid, vmid = _midpoints(vel[1:lenx÷2],vel[lenx÷2+1:lenx],C)
     u .= umid
     v .= vmid

     # Rotate to the inertial coordinate system
     T = RigidTransform(b.cent,b.α)
     for i in eachindex(u)
         Utmp = rot*[u[i],v[i]]
         u[i], v[i] = Utmp
     end

     return u, v
end


"""
    motion_velocity(b::Body,m::BasicDirectMotion,t::Real)

Return the velocity components (as a vector) of a `BasicDirectMotion`
at the given time `t`. These specify the velocities (in body cooordinates)
of the surface segments endpoints.
"""
function motion_velocity(b::Body,m::BasicDirectMotion,t::Real)
    return vcat(m.u,m.v)
end


"""
    motion_state(b::Body,m::AbstractDirectlySpecifiedMotion)

Return the current state vector of body `b` associated with
direct motion `m`. It returns the concatenated coordinates
of the body surface (in the body-fixed coordinate system).
"""
function motion_state(b::Body,m::AbstractDirectlySpecifiedMotion)
    return vcat(b.x̃end,b.ỹend)
end

"""
    update_body!(b::Body,x::AbstractVector,m::AbstractDirectlySpecifiedMotion)

Update body `b` with the motion state vector `x`, interpreted as coordinates in the body
coordinate system. The information in `m` is used for parsing only.
"""
function update_body!(b::Body{N,C},x::AbstractVector,m::AbstractDirectlySpecifiedMotion) where {N,C}
    length(x) == length(motion_state(b,m)) || error("wrong length for motion state vector")

    lenx = length(x)
    b.x̃end .= x[1:lenx÷2]
    b.ỹend .= x[lenx÷2+1:lenx]

    b.x̃, b.ỹ = _midpoints(b.x̃end,b.ỹend,C)

    # use the existing rigid transform of the body to update the
    # inertial coordinates of the surface
    T = RigidTransform(b.cent,b.α)
    T(b)

    return b

end

#=
AbstractRigidAndDirectMotion describes motions that superpose the rigid-body
motion with surface deformation. For this type of motion, the velocity
is described by the usual rigid-body components (reference point velocity,
angular velocity), plus vectors ũ and ṽ, describing the surface endpoint velocity
*in the body coordinate system*.

The motion state consists of the centroid, the angle, and the positions
x̃end and ỹend of the surface endpoints in the body coordinate system.

To create a motion of this type, we still need to supply an extension of
motion_velocity(b,m,t). However, this needs to supply only the deforming
part of the velocity, and in the body's own coordinate system.
=#

"""
    RigidAndDirectMotion(rig::RigidBodyMotion,def::AbstractDirectlySpecifiedMotion)

Create an instance of basic superposition of a rigid-body motion
and directly-specified deformation velocity in body coordinates.
"""
struct RigidAndDirectMotion{RT,DT} <: AbstractMotion
    rigidmotion :: RT
    defmotion :: DT
end

"""
    RigidAndDirectMotion(kin::Kinematics,def::AbstractDirectlySpecifiedMotion)

Create an instance of basic superposition of a rigid-body motion with kinematics `kin`,
and directly-specified deformation velocity in body coordinates.
"""
RigidAndDirectMotion(kin::Kinematics,def::AbstractDirectlySpecifiedMotion) =
                            RigidAndDirectMotion(RigidBodyMotion(kin),def)

"""
    RigidAndDirectMotion(kin::Kinematics,ũ::Vector{Float64},ṽ::Vector{Float64})

Create an instance of basic superposition of a rigid-body motion and
directly-specified (constant) deformation velocity in body coordinates, to be associated with a body whose length
is the same as `ũ` and `ṽ`.
"""
RigidAndDirectMotion(kin::Kinematics, ũ, ṽ) = RigidAndDirectMotion(RigidBodyMotion(kin),
                                                                  BasicDirectMotion(ũ,ṽ))
"""
    RigidAndDirectMotion(ċ,α̇,ũ::Vector{Float64},ṽ::Vector{Float64})

Specify constant translational `ċ` and angular `α̇` velocity and
directly-specified (constant) deformation velocity in body coordinates, to be associated with a body whose length
is the same as `ũ` and `ṽ`.
"""
RigidAndDirectMotion(ċ, α̇, ũ, ṽ) = RigidAndDirectMotion(RigidBodyMotion(ċ, α̇),
                                                        BasicDirectMotion(ũ,ṽ))

"""
    surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                     b::Body,motion::RigidAndDirectMotion,t::Real)

Assign the components of velocity `u` and `v` (in inertial coordinate system)
at surface positions described by points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body. This function calls the supplied
function for the deformation part in `motion.defmotion`.
"""
function surface_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                           b::Body,m::RigidAndDirectMotion,t::Real)

     surface_velocity!(u, v, b, m.defmotion, t)

     # Add the rigid part
     urig, vrig = similar(u), similar(v)
     surface_velocity!(urig,vrig,b,m.rigidmotion,t)

     u .+= urig
     v .+= vrig

     return u, v
end


"""
    motion_velocity(b::Body,m::RigidAndDirectMotion,t::Real)

Return the velocity components (as a vector) of a `RigidAndDirectMotion`
at the given time `t`.
"""
@inline motion_velocity(b::Body,m::RigidAndDirectMotion,t::Real) =
          vcat(motion_velocity(b,m.rigidmotion,t),motion_velocity(b,m.defmotion,t))


"""
    motion_state(b::Body,m::RigidAndDirectMotion)

Return the current state vector of body `b` associated with
rigid+direct motion `m`. It returns the concatenated coordinates
of the rigid-body mode and the body surface (in the body coordinate system).
"""
@inline motion_state(b::Body,m::RigidAndDirectMotion) =
          vcat(motion_state(b,m.rigidmotion),motion_state(b,m.defmotion))


"""
    update_body!(b::Body,x::AbstractVector,m::RigidAndDirectMotion)

Update body `b` with the motion state vector `x`. The part of the motion state
associated with surface deformation is interpreted as expressed in body coordinates.
The information in `m` is used for parsing only.
"""
function update_body!(b::Body{N,C},x::AbstractVector,m::RigidAndDirectMotion) where {N,C}
    length(x) == length(motion_state(b,m)) || error("wrong length for motion state vector")

    lenrigx = length(motion_state(b,m.rigidmotion))
    lendefx = length(x) - lenrigx

    b.x̃end .= x[lenrigx+1:lenrigx+lendefx÷2]
    b.ỹend .= x[lenrigx+lendefx÷2+1:length(x)]
    b.x̃, b.ỹ = _midpoints(b.x̃end,b.ỹend,C)

    update_body!(b,x[1:lenrigx],m.rigidmotion)

    return b

end


#=
"""
    surface_velocity(b::Body,motion::AbstractDirectlySpecifiedMotion,t::Real)

Return the components of velocities (in inertial components) at surface positions
described by points in body `b` (also in inertial coordinate system) at time `t`,
based on supplied motion `motion` for the body.
"""
surface_velocity(b::Body,m::AbstractDirectlySpecifiedMotion,t::Real) =
                surface_velocity!(similar(b.x),similar(b.y),b,m,t)

=#
