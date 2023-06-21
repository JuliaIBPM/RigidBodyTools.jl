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
````

Before we construct the system, we need to provide a function that will
specify the exogenous y acceleration. By default, it sets it to zero. We
will override that with a function that sets it to a random value chosen from a
normal distribution. Note that the function must be mutating and have a signature
`(a,x,p,t)`, where `a` is the exogenous acceleration vector. The arguments `x`
and `p` are a state and a parameter, which can be flexibly defined. The last
argument `t` is time. Here, we don't need any of those arguments.

````@example exogenous
function my_exogenous_function!(a,x,p,t)
    a .= randn(length(a))
end
````

We pass that along via the `exogenous` keyword argument.

````@example exogenous
ls = RigidBodyMotion(joints,bodies;exogenous=my_exogenous_function!)
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
acceleration. Let's advance the system and animate it. We include
a horizontal line along the hinge axis to show the effect of the exogenous
motion.

````@example exogenous
xhist = []
@gif for t in t0:dt:t0+tmax

  motion_rhs!(dxdt,x,(ls,bc),t)
  global x += dxdt*dt
  update_body!(bc,x,ls)

  push!(xhist,copy(x))
  plot(bc,xlims=(-1,5),ylims=(-1.5,1.5))
  hline!([0.0])
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

