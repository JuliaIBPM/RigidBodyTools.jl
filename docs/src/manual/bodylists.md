```@meta
EditURL = "../../../test/literate/bodylists.jl"
```

# Lists of bodies and their transforms

```@meta
CurrentModule = RigidBodyTools
```

We might want to have several distinct bodies. Here, we discuss how
to combine bodies into lists, and similarly, their transforms.

````@example bodylists
using RigidBodyTools
using Plots
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

We can transform the list by creating a list of transforms with a `MotionTransformList`

````@example bodylists
X1 = MotionTransform([2.0,3.0],0.0)
X2 = MotionTransform([-2.0,-0.5],π/4)
tl = MotionTransformList([X1,X2])
````

The transform list can be applied to the whole body list simply with

````@example bodylists
tl(bl)
````

which creates a copy of the body list and transforms that, or

````@example bodylists
update_body!(bl,tl)
````

which updates each body in `bl` in place.

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

## Body and transform list functions
```@docs
BodyList
getrange
Base.collect(::BodyList)
Base.sum(::AbstractVector,::BodyList,::Int)
Base.view(::AbstractVector,::BodyList,::Int)
MotionTransformList
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

