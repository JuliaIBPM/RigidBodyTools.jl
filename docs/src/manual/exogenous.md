```@meta
EditURL = "<unknown>/literate/exogenous.jl"
```

# Exogenous degrees of freedom

```@meta
CurrentModule = RigidBodyTools
```

As we mentioned in previous pages, some degrees of freedom can be
designated as *exogenous*, meaning that their behavior is determined
by some external process that we do not model explicitly. In practice,
that means that the acceleration of such a degree of freedom must
by explicitly provided at every time step while the state vector is being
advanced.

````@example exogenous
using RigidBodyTools
using Plots
````

As usual, we will demonstrate this via an example. In the example,
a single flat plate will be commanded to pitch upward by 45 degrees
about its leading edge. It will move steadily in the +x direction.
However, its y acceleration will vary randomly via some exogenous process.

````@example exogenous
Xp_to_jp = MotionTransform(0.0,0.0,0.0)
Xc_to_jc = MotionTransform(0.5,0.0,0.0)
dofs = [SmoothRampDOF(0.4,π/4,0.5),ConstantVelocityDOF(1.0),ExogenousDOF()]
joint = Joint(FreeJoint2d,0,Xp_to_jp,1,Xc_to_jc,dofs)
joints = [joint]

b = ThickPlate(1.0,0.05,0.02)
bodies = BodyList([b])

ls = RigidBodyMotion(joints,bodies)
````

Let's initialize the state vector and its rate of change

````@example exogenous
bc = deepcopy(bodies)
dt, tmax = 0.01, 3.0
t0, x0 = 0.0, init_motion_state(bc,ls)
dxdt = zero(x0)
x = copy(x0)
````

Note that the state vector has four elements. The first two are
associated with the prescribed motions for rotation and x translation.
The third is the y position, the exogenous degree of freedom. And the
fourth is the y velocity.

Why the y velocity? Because the exogenous behavior is specified via its
acceleration. This acceleration will need to be provided at each step in
the time marching. To help with this, we create a zero vector:

````@example exogenous
a_edof = zero_joint(ls,dimfcn=exogenous_dimension)
a_udof = zero_joint(ls,dimfcn=unconstrained_dimension);
nothing #hide
````

Now, `a_edof` is not empty, as it was in the previous example, but has a single element. We will set this
element's value inside the loop, using a random value chosen from a
normal distribution. We will record the history of the state while
we advance it

````@example exogenous
xhist = []
@gif for t in t0:dt:t0+tmax
  a_edof[1] = randn()

  motion_rhs!(dxdt,x,t,a_edof,a_udof,ls,bc)
  global x += dxdt*dt
  update_body!(bc,x,ls)

  push!(xhist,copy(x))
  plot(bc,xlims=(-1,5),ylims=(-1.5,1.5))
end every 5
````

Let's plot the exogenous state and its velocity

````@example exogenous
plot(t0:dt:t0+tmax,map(x -> x[3],xhist),label="y position",xlabel="t")
plot!(t0:dt:t0+tmax,map(x -> x[4],xhist),label="y velocity")
````

The variation in velocity is quite noisy (and constitutes a random walk).
In contrast, the change in position is relatively smooth, since it represents
an integral of this velocity.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

