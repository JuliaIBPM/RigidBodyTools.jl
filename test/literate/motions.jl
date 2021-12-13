# # Making bodies move

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
For many applications we wish to specify the motion of a body. Motions
come in two forms: rigid-body motions, in which the body moves as a rigid
unit; and deformations, in which the surface's shape is altered. These
two types of motions can also be combined.
=#


using RigidBodyTools
using Plots

#=
Before we get started, let's define a macro that can help us
animate the motions we will be creating
=#
macro animate_motion(b,m,dt,tmax,xlim,ylim)
    return esc(quote
            bc = deepcopy($b)
            t0, x0 = 0.0, motion_state(bc,$m)
            x = copy(x0)
            @gif for t in t0:$dt:t0+$tmax
                global x += motion_velocity(bc,$m,t)*$dt
                update_body!(bc,x,$m)
                plot(bc,xlim=$xlim,ylim=$ylim)
            end every 5
        end)
end
#=
## Rigid-body motions
We will start with rigid-body motions. These are created with the `RigidBodyMotion`
constructor. For these, we have a variety
of *kinematics* classes, which allow us to constrain the motion in particular
ways. The most basic of these is just constant translational and rotational motion:
=#
m = RigidBodyMotion((1.0,2.0),π/2)

#=
Let's demonstrate this on a rectangular body
=#
b = Rectangle(1.0,0.5,0.02)
@animate_motion b m 0.01 2.0 (-2,5) (-2,5)

#=
Notice that the animation macro makes use of three key functions:
`motion_state`, `motion_velocity`, and `update_body!`. These are used
frequently:
- `motion_state(b,m)` returns the state vector that describes the body motion
- `motion_velocity(b,m,t)` returns the rate of change of the state vector at time `t`
- `update_body!(b,x,m)` updates the body with the state vector `x`.

There is one more function that is useful for some downstream application:
`surface_velocity!(u,v,b,m,t)` provides the surface velocity components `u` and `v`
*in the inertial coordinate system, at the midpoints of the surface segments*.
=#

#=
Let's try something more interesting. We will use a pitch-up maneuver on a
thick plate. This pitches up about the leading edge to 45 degrees, at a pitch
rate of 0.2 (the ratio $\dot{\alpha}c/2U_0$), where $c$ is the length of the
plate. The forward velocity $U_0$ is set to 1.
=#
kin = Pitchup(1.0,0.5,0.2,0.0,0.5,π/4)
m = RigidBodyMotion(kin)
b = ThickPlate(1.0,0.05,0.02)
@animate_motion b m 0.01 3.0 (-1,5) (-1,1)

#=
There are a number of other types of kinematics, particularly those
with oscillatory behavior: `OscillationX`, `OscillationY`, `RotationalOscillation`,
`OscillationXY`, `PitchHeave`, `Oscillation`. Each one is documented.
=#

#=
## Deforming bodies
For deforming bodies, we specify the motion of the surface directly. This
deformation velocity is expressed in the coordinate system attached to the
body, rather than our inertial coordinate system. This enables the
motion to be easily superposed with the rigid-body motion described above.

It is also important to note that the motion is applied **to the endpoints**
of the surface segments. The midpoints are then constructed from the
updated endpoints.
=#
#=
Let's see an example. We will create an oscillatory deformation of a circle.
We create the motion by creating functions for each component of velocity.
=#
ufcn(x,y,t) = 0.25*x*y*cos(t)
vfcn(x,y,t) = 0.25*(x^2-y^2)*cos(t)
m = DeformationMotion(ufcn,vfcn)
#=
Now create the body and visualize the motion
=#
b = Circle(1.0,0.02)
@animate_motion b m π/100 4π (-2,2) (-2,2)

#=
There is a more basic type of deformation motion: `ConstantDeformationMotion`,
which is not time varying, but simply a constant vector for each surface
(end) point. Let's make the circle expand radially at constant velocity
=#
u = copy(b.x̃end)
v = copy(b.ỹend)
m = ConstantDeformationMotion(u,v)
@animate_motion b m 0.01 2.0 (-5,5) (-5,5)


#=
## Rigid and deforming motions
We can easily combine rigid and deforming motions with the use of `RigidAndDeformingMotion`.
Let's make a square oscillate in rotation and undergo oscillatory deformation:
=#
b = Square(1.0,0.02)

#=
Here is the rotational oscillation, with frequency $\Omega$.
=#
Ω = 1.0
kin = RotationalOscillation(Ω,π/4,0.0)
mrig = RigidBodyMotion(kin)

#=
Here is the deformation.
=#
ufcn(x,y,t) = 0.25*(x^2+y^2)*y*cos(Ω*t)
vfcn(x,y,t) = -0.25*(x^2+y^2)*x*cos(Ω*t)
mdef = DeformationMotion(ufcn,vfcn)

#=
Now put them together with `RigidAndDeformingMotion`, and animate it:
=#
m = RigidAndDeformingMotion(mrig,mdef)
@animate_motion b m π/100 4π (-2,2) (-2,2)

#=
## Defining new motions
It is straightforward to define new types of deformation motion that don't fit into the framework
shown here. We need only do two things:
- Create a type as a subtype of `AbstractDeformationMotion`.
- Extend the function `RigidBodyTools.motion_velocity(b,m,t)` for your
  new motion type, so that it returns a concatenated vector of
  the surface segment endpoint velocity components (in body coordinate system).
=#


#md # ## Motion types
#md # ```@docs
#md # DeformationMotion
#md # ConstantDeformationMotion
#md # RigidBodyMotion
#md # RigidBodyMotion(::Kinematics)
#md # RigidBodyMotion(::Any,::Any)
#md # RigidAndDeformingMotion
#md # RigidAndDeformingMotion(::Kinematics,::AbstractDirectlySpecifiedMotion)
#md # RigidAndDeformingMotion(::Kinematics,::Any,::Any)
#md # RigidAndDeformingMotion(::Any,::Any,::Any,::Any)
#md # ```

#md # ## Surface and motion velocity functions
#md # ```@docs
#md # motion_state
#md # update_body!
#md # motion_velocity
#md # surface_velocity!
#md # surface_velocity
#md # ```



#md # ## Rigid body kinematics types
#md # ```@docs
#md # Kinematics
#md # Oscillation
#md # OscillationX
#md # OscillationY
#md # OscillationXY
#md # PitchHeave
#md # Pitchup
#md # RotationalOscillation
#md # SwitchedKinematics
#md # ```
