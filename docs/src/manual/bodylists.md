```@meta
EditURL = "<unknown>/literate/bodylists.jl"
```

# Lists of bodies

```@meta
CurrentModule = RigidBodyTools
```

We might want to have several distinct bodies. Here, we discuss how
to combine bodies into lists, and similarly, their motions and transforms.

````@example bodylists
using RigidBodyTools
using Plots
````

Once again, we will use our animation macro:

````@example bodylists
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
````

## Body list
Suppose we have two bodies and we wish to combine them into a single list.
The advantage of doing so is that many of the operations we have presented
previously also extend to lists. We use `BodyList` to combine them.

````@example bodylists
b1 = Circle(1.0,0.02)
b2 = Rectangle(1.0,2.0,0.02)
bl = BodyList([b1,b2])
````

Another way to do this is to push each one onto the list:

````@example bodylists
bl = BodyList()
push!(bl,b1)
push!(bl,b2)
````

We can transform the list by creating a list of transforms with a `RigidTransformList`

````@example bodylists
T1 = RigidTransform((2.0,3.0),0.0)
T2 = RigidTransform((-2.0,-0.5),π/4)
tl = RigidTransformList([T1,T2])
````

The transform list can be applied to the whole body list simply with

````@example bodylists
tl(bl)
````

Let's see our effect

````@example bodylists
plot(bl)
````

It is important to note that the list points to the original bodies,
so that any change made to the list is reflected in the original bodies, e.g.

````@example bodylists
plot(b2)
````

## Utilities on lists
There are some specific utilities that are helpful for lists. For example,
to collect all of the x, y points (the segment midpoints) in the list into two
vectors, use

````@example bodylists
x, y = collect(bl)
````

In a vector comprising data on these concatenated surface points, we
can use `view` to look at just one body's part and change it:

````@example bodylists
f = zero(x)
f1 = view(f,bl,1)
f1 .= 1.0;
plot(f)
````

Also, we can sum up the values for one of the bodies:

````@example bodylists
sum(f,bl,2)
````

## Motion lists
Motions can also be assembled into lists, and most of the operations
on them extend to lists. Let's create a list of motions: one for body
1 and one for body 2.

````@example bodylists
ufcn(x,y,t) = 0.25*x*y*cos(t)
vfcn(x,y,t) = 0.25*(x^2-y^2)*cos(t)
m1 = DeformationMotion(ufcn,vfcn)
Ω = 1.0
kin = RotationalOscillation(Ω,π/4,0.0)
m2 = RigidBodyMotion(kin)
ml = MotionList([m1,m2]);
nothing #hide
````

Now let's see it in action

````@example bodylists
@animate_motion bl ml π/100 4π (-5,5) (-5,5)
````

We can also use the `surface_velocity!` function to get the velocities
of all surface points in the list. For example, to get them for
this previous list at time $t = 1.0$

````@example bodylists
x, y = collect(bl)
u, v = zero(x), zero(y)
t = 1.0
surface_velocity!(u,v,bl,ml,t)
````

We can determine the maximum velocity across the whole set of bodies:

````@example bodylists
umax, i, tmax, bmax = maxlistvelocity(bl,ml)
````

In this case, the maximum velocity occurs at t = 0 on body 2, index 301.

## Body list functions
```@docs
BodyList
getrange
Base.collect(::BodyList)
Base.sum(::AbstractVector,::BodyList,::Int)
Base.view(::AbstractVector,::BodyList,::Int)
MotionList
RigidTransformList
Base.vec(::RigidTransformList)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

