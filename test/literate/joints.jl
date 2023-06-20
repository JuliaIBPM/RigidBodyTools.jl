# # Joints

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
=#

#=
### Example
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
Put it in a one-element vector
=#
dofs = [kr];
#=
and construct the joint
=#
joint2 = Joint(RevoluteJoint,pid,Xp_to_jp,cid,Xc_to_jc,dofs)

#=
Group the two joints together into a vector. It doesn't matter what
order this vector is.
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
form of the `view` function for this. For example, to get access to just
the part of the state vector for joint 2,
=#
jid = 2
x2 = view(x,ls,jid)

#=
Then, you might decide to change the entry of `x2`, which, in turn,
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

Using a simple forward Euler method, the state vector can be advanced
as follows
=#
t0, x0 = 0.0, init_motion_state(bodies,ls)
dxdt = zero(x0)
x = copy(x0)
dt, tmax = 0.01, 4.0
for t in t0:dt:t0+tmax
  motion_rhs!(dxdt,x,t,[],[],ls,bodies)
  global x += dxdt*dt
end

#=
Note that `motion_rhs!` takes two arguments that are given here as empty
vectors. Those represent, respectively, the accelerations of the exogenous and unconstrained
degrees of freedom. In this example, there are no such degrees of freedom. However,
it is useful to take this opportunity to point out the function `zero_joint`, which can create
a vector of zeros of length necessary for several different dimensional aspects of the
joint's state. For example, a vector of zeros for the joint's exogenous
degree-of-freedom dimensionality would be generated with
=#
zero_joint(ls,dimfcn=exogenous_dimension)

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

            a_edof = zero_joint(ls,dimfcn=exogenous_dimension)
            a_udof = zero_joint(ls,dimfcn=unconstrained_dimension)


            @gif for t in t0:$dt:t0+$tmax
                motion_rhs!(dxdt,x,t,a_edof,a_udof,$m,bc)
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
