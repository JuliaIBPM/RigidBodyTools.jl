```@meta
EditURL = "<unknown>/literate/deformation.jl"
```

# Deforming bodies

```@meta
CurrentModule = RigidBodyTools
```

Thus far we have only shown rigid body motion. However, we can
also prescribe surface deformation as an additional component
of a body's motion.

````@example deformation
using RigidBodyTools
using Plots
````

Before we get started, let's define the same macro that we used earlier
in order to visualize our system's motion

````@example deformation
macro animate_motion(b,m,dt,tmax,xlim,ylim)
    return esc(quote
            bc = deepcopy($b)
            t0, x0 = 0.0, init_motion_state(bc,$m)
            dxdt = zero(x0)
            x = copy(x0)

            @gif for t in t0:$dt:t0+$tmax
                motion_rhs!(dxdt,x,($m,bc),t)
                global x += dxdt*$dt
                update_body!(bc,x,$m)
                plot(bc,xlims=$xlim,ylims=$ylim)
            end every 5
        end)
end
````

For deforming bodies, we specify the velocity of the surface directly. This
deformation velocity is expressed in the coordinate system attached to the
body, rather than the inertial coordinate system. This enables the
motion to be easily superposed with the rigid-body motion described earlier.

It is also important to note that the motion is applied **to the endpoints**
of the surface segments. The midpoints are then constructed from the
updated endpoints.
## Example: Basic deformation
Let's see an example. We will create an oscillatory deformation of a circle.
We create the motion by creating functions for each component of velocity.

````@example deformation
Ω = 2π
ufcn(x,y,t) = 0.25*x*y*Ω*cos(Ω*t)
vfcn(x,y,t) = 0.25*(x^2-y^2)*Ω*cos(Ω*t)
def = DeformationMotion(ufcn,vfcn)
````

We will create a simple fixed revolute joint that anchors the body's center to the
inertial system. (Note that we don't need to create a body list here, since we
are only working with one body and one joint.)

````@example deformation
Xp_to_jp = MotionTransform(0.0,0.0,0.0)
Xc_to_jc = MotionTransform(0.0,0.0,0.0)
dofs = [ConstantVelocityDOF(0.0)]
joint = Joint(RevoluteJoint,0,Xp_to_jp,1,Xc_to_jc,dofs)

body = Circle(1.0,0.02)
````

To construct the system, we supply the joint and body, as before, as well as the deformation.

````@example deformation
ls = RigidBodyMotion(joint,body,def)
````

Let's animate this motion

````@example deformation
@animate_motion body ls 0.01 4 (-2,2) (-2,2)
````

The body remains fixed, but the surface deforms!

## Example: Expanding motion
Now a circle undergoing an expansion. For this, we set constant velocity
components equal to constants, the coordinates of the surface segment endpoints

````@example deformation
body = Circle(1.0,0.02)

u = copy(body.x̃end)
v = copy(body.ỹend)
def = ConstantDeformationMotion(u,v)

ls = RigidBodyMotion(joint,body,def)
@animate_motion body ls 0.01 2 (-5,5) (-5,5)
````

## Example: Combining rigid motion and deforming motion.
Now, let's combine an oscillatory rigid-body rotation with
oscillatory deformation, this time applied to a square.

````@example deformation
Xp_to_jp = MotionTransform(0.0,0.0,0.0)
Xc_to_jc = MotionTransform(0.0,0.0,0.0)
Ω = 1.0
dofs = [OscillatoryDOF(π/4,Ω,0.0,0.0)]
joint = Joint(RevoluteJoint,0,Xp_to_jp,1,Xc_to_jc,dofs)

body = Square(1.0,0.02)

ufcn(x,y,t) = 0.25*(x^2+y^2)*y*Ω*cos(Ω*t)
vfcn(x,y,t) = -0.25*(x^2+y^2)*x*Ω*cos(Ω*t)
def = DeformationMotion(ufcn,vfcn)

ls = RigidBodyMotion(joint,body,def)

@animate_motion body ls π/100 4π (-2,2) (-2,2)
````

## Example: Defining new deformations
We can also define new types of deformation that are more specialized.
We need only define a subtype of `AbstractDeformationMotion`
and extend the function `deformation_velocity` to work with it.
The signature of this function is `deformation_velocity(body,deformation,time)`.

For example, let's define a motion on a rectangular shape that
will deform only the top side in the normal direction, but leave the rest of
the surface stationary. We will use the `side` field
of the `Polygon` shape type to access the top, and set
its vertical velocity.

````@example deformation
struct TopMotion{UT} <: AbstractDeformationMotion
    vtop :: UT
end

function RigidBodyTools.deformation_velocity(body::Polygon,def::TopMotion,t::Real)

    u, v = zero(body.x̃end), zero(body.ỹend)
    top = body.side[3]
    v[top] .= def.vtop.(body.x̃end[top],body.ỹend[top],t)
    return vcat(u,v)
end
````

Now apply it

````@example deformation
Xp_to_jp = MotionTransform(0.0,0.0,0.0)
Xc_to_jc = MotionTransform(0.0,0.0,0.0)
dofs = [ConstantVelocityDOF(0.0)]
joint = Joint(RevoluteJoint,0,Xp_to_jp,1,Xc_to_jc,dofs)

body = Rectangle(1.0,2.0,0.02)
vfcn(x,y,t) = 0.2*(1-x^2)*cos(t)
def = TopMotion(vfcn)

ls = RigidBodyMotion(joint,body,def)
````

Let's try it out

````@example deformation
@animate_motion body ls π/100 4π (-1.5,1.5) (-2.5,2.5)
````

As desired, the top surface deforms vertically, but the rest of the
surface is stationary.

## Deformation functions
```@docs
DeformationMotion
ConstantDeformationMotion
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

