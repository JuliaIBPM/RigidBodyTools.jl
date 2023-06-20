# # Behaviors of degrees of freedom

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
To define the motion of a joint requires that we define the behavior
of each of the joint's degrees of freedom. There are three different types of
behavior of a degree of freedom:
(1) Its motion is prescribed with a given function of time, in which case we say that the
    degree of freedom is *constrained*
(2) Its motion is specified directly, but determined by a process that lies outside
    of the body system, in which case we say the degree of freedom is *exogenous*.
(3) Its motion is determined indirectly by some sort of force model for its behavior,
    such as a spring or damper, in which case we say the degree of freedom is *unconstrained*.

Here, we will discuss how to set a degree of freedom with one of these behaviors,
and how to set the prescribed motion, if desired.
=#

using RigidBodyTools
using Plots

#=
## Constrained behavior (i.e., prescribed motion)
When a degree of freedom is constrained, then its behavior over time
is set by some time-varying function. There are a number of different
types of predefined time-varying behaviors.
=#
#=
### Constant velocity
To specify constant velocity with some velocity, e.g. $U=1$, set
=#
U = 1.0
k = ConstantVelocityDOF(U)

#=
For any prescribed motion, we can evaluate it at any specified time.
It returns data of type `DOFKinematicData`.
=#
t = 1.0
kt = k(t)

#=
The kinematic data can be parsed with
`dof_position`, `dof_velocity`, and `dof_acceleration`:
=#
dof_position(kt)
#-
dof_velocity(kt)
#-
dof_acceleration(kt)

#=
Let's plot the position over time
=#
t = range(0,3,length=301)
plot(t,dof_position.(k.(t)),xlims=(0,Inf),ylims=(0,Inf))

#=
### Oscillatory motion
We can set the position to be a sinusoidal function. For this, we
set the amplitude, the angular frequency, the phase, and the mean velocity (typically zero).
=#
A = 1.0 ## amplitude
Ω = 2π ## angular frequency
ϕ = π/2 ## phase
vel = 0  ## mean velocity
k = OscillatoryDOF(A,Ω,ϕ,vel)

#=
Plot the position, velocity, and acceleration
=#
plot(t,dof_position.(k.(t)),xlims=(0,Inf),label="x")
plot!(t,dof_velocity.(k.(t)),label="ẋ")
plot!(t,dof_acceleration.(k.(t)),label="ẍ")

#=
### Smooth ramp motion
To ramp the position from one value to another, we use the `SmoothRampDOF`.
For this, we need to specify the nominal velocity of the ramp,
the change in position, and the time at which the ramp starts.
There is an optional argument `ramp` to control the ramp's smoothness.
It defaults to `EldredgeRamp(11.0)`, an Eldredge-type ramp with smoothness
factor 11.
=#
vel = 1.0  ## nominal ramp velocity
Δx = 1.0 ## change in position
t0 = 1.0 ## time of ramp start
k = SmoothRampDOF(vel,Δx,t0)

#=
Plot the position
=#
plot(t,dof_position.(k.(t)),xlims=(0,Inf),label="x")

#=
### User-defined motion
The user can specify the time-varying position by supplying a function of time
and using `CustomDOF`. It automatically differentiates this function to
get velocity and acceleration.
For example, a quadratic behavior
=#
f(t) = 1.0*t + 2.0*t^2
k = CustomDOF(f)

#=
Plot the position
=#
plot(t,dof_position.(k.(t)),xlims=(0,Inf),label="x")

#=
and the velocity
=#
plot(t,dof_velocity.(k.(t)),xlims=(0,Inf),label="ẋ")


#=
and the acceleration
=#
plot(t,dof_acceleration.(k.(t)),xlims=(0,Inf),label="ẍ")


#=
## Exogenous and unconstrained behaviors
If the degree of freedom is to be *exogenous* or *unconstrained*,
then it can be designated as such, e.g,
=#
k = ExogenousDOF()
#=
or
=#
k = UnconstrainedDOF()


#md # ## Degree of freedom functions
#md # ```@docs
#md # ConstantVelocityDOF
#md # OscillatoryDOF
#md # SmoothRampDOF
#md # CustomDOF
#md # ExogenousDOF
#md # UnconstrainedDOF
#md # dof_position
#md # dof_velocity
#md # dof_acceleration
#md # ```
