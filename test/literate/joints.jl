# # Joints and body-joint systems

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
For systems consisting of one or more rigid bodies, *joints* are used to specify
the degrees of freedom permitted for the body. Even for a single rigid body,
the body is assumed to be connected by a joint with the inertial coordinate system.

In three dimensions, there are several different classes of joint:
- `RevoluteJoint`, which rotates about an axis and has only one degree of freedom
- `PrismaticJoint`, which slides along one axis and has only one degree of freedom
- `HelicalJoint`, which rotates about and slides along one axis (two degrees of freedom)
- `SphericalJoint`, which rotates freely about one point (three degrees of freedom)
- `FreeJoint`, which can move in any manner (six degrees of freedom)

However, in two spatial dimensions, there are only two
- `RevoluteJoint` (one degree of freedom: rotation)
- `FreeJoint2d` (three degrees of freedom: rotation, two translation)

=#

using RigidBodyTools
using Plots

#=
A joint serves as a connection between two bodies; one of the bodies is the *parent*
and the other is the *child*. It is also common for a body to be connected to
the inertial coordinate system, in which case the inertial system is the parent
and the body is the child. The user must specify the id of the parent body and the child body.
The id of the inertial system is 0.

Before we discuss how a joint is specified, it is important to understand that
a joint is basically just an intermediate `MotionTransform` from the parent body's system
to the child body's system.

That is, to construct the transform from body P to body C, we would compose it
as follows

$${}^C X_{P} = {}^{C}X_{J(C)} {}^{J(C)}X_{J(P)} {}^{J(P)}X_{P}$$

The transform in the middle is the joint transform and is the only one of these that
can vary in time. It maps from one
coordinate system attached (rididly) to the parent body to another coordinate
system attached (rigidly) to the child body.

The degrees of freedom of the joint constrain the behavior of this joint transform.
Below, we show how joints are constructed and used.
=#

#=
## Joint construction generalities
The user must specify where the joint is located on each body, relative to the
origin of the body's coordinate system, and what orientation the joint takes
with respect to the body's coordinate system. These are specified with
`MotionTransform`s, transforming from each body's coordinate system to the joint's
coordinate system on that body. This transform is invariant -- it never changes.

Also, the user must specify the behavior of each of the degrees of freedom
of a joint, using the tools we discussed in the previous page.

The basic signature is

`Joint(joint_type,parent_body_id,Xp_to_jp,child_body_id,Xc_to_jc,doflist)`

However, there is a specialized signature if you simply wish to place
a body in some stationary configuration `X`:

`Joint(X,body_id)`

and if there is only one body and it is stationary, then

`Joint(X)`

will do.

=#

#=
## Example
Let's see a 2d example. Suppose we have two bodies, 1 and 2. Body 1 is to
be connected to the inertial coordinate system, and prescribed with motion
that causes it to oscillate rotationally (*pitch*) about a point located
at $(-1,0)$ and in the body's coordinate system, and this pitch axis
can oscillate up and down (*heave*).

Body 2 is to be connected by a hinge to body 1, and this hinge's angle will
also oscillate. We will set the hinge to be located at $(1.02,0)$
on body 1 and $(-1.02,0)$ on body 2.
=#
#=
### Joint from inertial system to body 1
First, let's construct the joint from the inertial system to body 1. This
should be a `FreeJoint2d`, since motion is prescribed in two of the three
degrees of freedom (as well as the third one, with zero velocity). We can assume
that the joint is attached to the origin of the parent (the inertial
system). On the child (body 1), the joint is to be at `[-1,0]`.
=#
pid = 0
Xp_to_jp = MotionTransform([0,0],0)
cid = 1
Xc_to_jc = MotionTransform([-1,0],0)

#=
For the rotational and the y degrees of freedom, we need oscillatory motion.
For the rotational motion, let's set the amplitude to $\pi/4$ and the
angular frequency to $\Omega = 2\pi$, but set the phase and mean velocity both
to zero.
=#
Ar = π/4
Ω = 2π
ϕr = 0.0
vel = 0.0
kr = OscillatoryDOF(Ar,Ω,ϕr,vel)

#=
For the plunging, we will set an amplitude of 1 and a phase lag of $\pi/2$, but keep the same
frequency as the pitching.
=#
Ay = 1
ϕy = -π/2
ky = OscillatoryDOF(Ay,Ω,ϕy,vel)

#=
The x degree of freedom is simply constant velocity, set to 0, to
ensure it does not move in the $x$ direction.
=#
kx = ConstantVelocityDOF(0)
#-
#=
We put these together into a vector, to pass along to the joint constructor.
The ordering of these in the vector is important. It must be
[rotational, x, y].
=#
dofs = [kr,kx,ky];

#=
Now set the joint
=#
joint1 = Joint(FreeJoint2d,pid,Xp_to_jp,cid,Xc_to_jc,dofs)

#=
Note that this joint has three constrained degrees of freedom,
no exogenous degrees of freedom, and no unconstrained degrees of freedom.
In a later example, we will change this.
=#

#=
### Joint from body 1 to body 2
This joint is a `RevoluteJoint`. First set the joint locations on each body.
=#
pid = 1
Xp_to_jp = MotionTransform([1.02,0],0)
cid = 2
Xc_to_jc = MotionTransform([-1.02,0],0)

#=
Now set its single degree of freedom (rotation) to have oscillatory
kinematics. We will set its amplitude the same as before,
but give it a phase lag
=#
Ar = π/4
Ω = 2π
ϕr = -π/4
kr = OscillatoryDOF(Ar,Ω,ϕr,vel)
#=
Put it in a one-element vector.
=#
dofs = [kr];
#=
and construct the joint
=#
joint2 = Joint(RevoluteJoint,pid,Xp_to_jp,cid,Xc_to_jc,dofs)

#=
Group the two joints together into a vector. It doesn't matter what
order this vector is in, since the connectivity will be figured out any way
it is ordered, but it numbers the joints by the order they are provided here.
=#
joints = [joint1,joint2];

#=
## Assembling the joints and bodies
The joints and bodies comprise an overall `RigidBodyMotion` system.
When this system is constructed, all of the connectivities are
determined (or missing connectivities are revealed). The
construction requires that we have set up the bodies themselves, so
let's do that first. We will make both body 1 and body 2 an ellipse
of aspect ratio 5. Note that ordering of bodies matters here, because
the first in the list is interpreted as body 1, etc.
=#
b1 = Ellipse(1.0,0.2,0.02)
b2 = Ellipse(1.0,0.2,0.02)
bodies = BodyList([b1,b2])

#=
Now we can construct the system
=#
ls = RigidBodyMotion(joints,bodies)

#=
Before we proceed, it is useful to demonstrate some of the tools
we have to probe the connectivity of the system. For example,
to find the parent joint of body 1, we use
=#
parent_joint_of_body(1,ls)

#=
This returns 1, since we have connected body 1 to joint 1.
How about the child joint of body 1? We expect it to be 2. Since there
might be more than one child, this returns a vector:
=#
child_joints_of_body(1,ls)

#=
We can also check the body connectivity of joints. This can be very useful
for more complicated systems in which the joint numbering is less clear.
The parent body of joint 1
=#
parent_body_of_joint(1,ls)
#=
This returns 0 since we have connected joint 1 to the inertial system. The child body of joint 2:
=#
child_body_of_joint(2,ls)

#=
## The system state vector
A key concept in advancing, plotting, and doing further analysis of this
system of joints and bodies is the *state vector*, $x$. This state vector
has entries for the position of every degree of freedom of all of the joints. It also may
have further entries for other quantities that need to be advanced, but for
this example, there are no other entries.

There are two functions that are useful for constructing the state vector. The
first is `zero_motion_state`, which simply creates a vector of zeros of the
correct size.
=#
zero_motion_state(bodies,ls)
#=
The second is `init_motion_state`, which fills in initial position values for any degrees of freedom
that have been prescribed.
=#
x = init_motion_state(bodies,ls)

#=
Note that neither of these functions has any mutating effect on the arguments (`bodies` and `ls`).

Also, it is always possible for the user to modify the entries in the state
vector after this function is called. In general, it would be difficult to
determine which entry is which in this state vector, so we can use a special
function for this. For example, to get access to just
the part of the state vector for the positions of joint 1,
=#
jid = 1
x1 = position_vector(x,ls,jid)

#=
This is a view on the overall state vector. This, if you decide to change an entry of `x1`, this, in turn,
would change the correct entry in `x`.
=#

#=
We can use the system state vector to put the bodies in their proper places,
using the joint positions in `x`.
=#
update_body!(bodies,x,ls)

#=
Let's plot this just to check
=#
plot(bodies,xlims=(-4,4),ylims=(-4,4))

#=
## Advancing the state vector
Once the initial state vector is constructed, then the system can
be advanced in time. In this example, there are no exogenous or
unconstrained degrees of freedom that require extra input, so
the system is *closed* as it is. To advance the system, we need
to solve the system of equations
$$\frac{\mathrm{d}x}{\mathrm{d}t} = f(x,t)$$

The function $f(x,t)$ describing the rate of change of $x$ is given by the
function `motion_rhs!`. This function mutates its first argument,
the rate-of-change vector `dxdt`, which can then be used to update the state.
The system and bodies are passed in as a tuple, followed by time.

Using a simple forward Euler method, the state vector can be advanced
as follows
=#
t0, x0 = 0.0, init_motion_state(bodies,ls)
dxdt = zero(x0)
x = copy(x0)
dt, tmax = 0.01, 4.0
for t in t0:dt:t0+tmax
  motion_rhs!(dxdt,x,(ls,bodies),t)
  global x += dxdt*dt
end

#=
Now that we know how to advance the state vector, let's create a macro
that can be used to make a movie of the evolving system.
=#
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

#=
Let's use it here
=#
@animate_motion bodies ls 0.01 4 (-4,4) (-4,4)


#=
## Outputting the surface velocities
For use in mechanics problems, it is important to be able to output
the velocity of the points on the surface of bodies at a given system state
and time. We use the function `surface_velocity!` for this.

First, initialize vectors for the `u` and `v` components in this 2d example,
using the `zero_body` function.
=#
u, v = zero_body(bodies), zero_body(bodies)

#=
Now evaluate the velocities at time 0, with the initial state
=#
surface_velocity!(u,v,bodies,x0,ls,t0)

#=
We can plot these on each body using the `view` function for `BodyList`.
For example, the vector of u velocities on body 2 is
=#
plot(view(u,bodies,2))

#=
In this use of `surface_velocities!`, we outputted the velocities in inertial
coordinates. There is a keyword `axes` that allows us to relax this.
By default, this is set to `:inertial`. However, we can also output them in their own body coordinates with
the keyword `axes=:body`,
=#
surface_velocity!(u,v,bodies,x0,ls,t0;axes=:body)
plot(view(u,bodies,2))

#=
We can also select only to evaluate a part of the body's motion on the
surface points, using the `motion_part` keyword. This keyword defaults to `:full`,
but we can also select `:angular` or `:linear`:
=#
surface_velocity!(u,v,bodies,x0,ls,t0;axes=:body,motion_part=:angular)
plot(view(u,bodies,2))


#md # ## Joint functions
#md # ```@docs
#md # Joint
#md # zero_joint
#md # init_joint
#md # ```

#md # ## System and state functions
#md # ```@docs
#md # RigidBodyMotion
#md # zero_motion_state
#md # init_motion_state
#md # Base.view(::AbstractVector,::RigidBodyMotion,::Int)
#md # position_vector
#md # velocity_vector
#md # deformation_vector
#md # exogenous_position_vector
#md # exogenous_velocity_vector
#md # unconstrained_position_vector
#md # unconstrained_velocity_vector
#md # motion_rhs!
#md # zero_exogenous
#md # update_exogenous!
#md # surface_velocity!
#md # maxvelocity
#md # ismoving(::RigidBodyMotion)
#md # ```

#md # ## Joint types
#md # ```@docs
#md # RevoluteJoint
#md # PrismaticJoint
#md # HelicalJoint
#md # SphericalJoint
#md # FreeJoint
#md # FreeJoint2d
#md # ```
