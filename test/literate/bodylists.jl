# # Lists of bodies and their transforms

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
We might want to have several distinct bodies. Here, we discuss how
to combine bodies into lists, and similarly, their transforms.
=#

using RigidBodyTools
using Plots


#=
## Body list
Suppose we have two bodies and we wish to combine them into a single list.
The advantage of doing so is that many of the operations we have presented
previously also extend to lists. We use `BodyList` to combine them.
=#
b1 = Circle(1.0,0.02)
b2 = Rectangle(1.0,2.0,0.02)
bl = BodyList([b1,b2])

#=
Another way to do this is to push each one onto the list:
=#
bl = BodyList()
push!(bl,b1)
push!(bl,b2)

#=
We can transform the list by creating a list of transforms with a `MotionTransformList`
=#
X1 = MotionTransform([2.0,3.0],0.0)
X2 = MotionTransform([-2.0,-0.5],π/4)
tl = MotionTransformList([X1,X2])

#=
The transform list can be applied to the whole body list simply with
=#
tl(bl)

#=
which creates a copy of the body list and transforms that, or
=#
update_body!(bl,tl)

#=
which updates each body in `bl` in place.
=#

#=
Let's see our effect
=#
plot(bl)

#=
It is important to note that the list points to the original bodies,
so that any change made to the list is reflected in the original bodies, e.g.
=#
plot(b2)

#=
## Utilities on lists
There are some specific utilities that are helpful for lists. For example,
to collect all of the x, y points (the segment midpoints) in the list into two
vectors, use
=#
x, y = collect(bl)

#=
In a vector comprising data on these concatenated surface points, we
can use `view` to look at just one body's part and change it:
=#
f = zero(x)
f1 = view(f,bl,1)
f1 .= 1.0;
plot(f)

#=
Also, we can sum up the values for one of the bodies:
=#
sum(f,bl,2)


#md # ## Body and transform list functions
#md # ```@docs
#md # BodyList
#md # getrange
#md # Base.collect(::BodyList)
#md # Base.sum(::AbstractVector,::BodyList,::Int)
#md # Base.view(::AbstractVector,::BodyList,::Int)
#md # MotionTransformList
#md # ```
