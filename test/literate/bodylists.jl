# # Lists of bodies

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
We might want to have several distinct bodies. Here, we discuss how
to combine bodies into lists, and similarly, their motions and transforms.
=#

using RigidBodyTools
using Plots

#=
Once again, we will use our animation macro:
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
We can transform the list by creating a list of transforms with a `RigidTransformList`
=#
T1 = RigidTransform((2.0,3.0),0.0)
T2 = RigidTransform((-2.0,-0.5),π/4)
tl = RigidTransformList([T1,T2])

#=
The transform list can be applied to the whole body list simply with
=#
tl(bl)

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

#=
## Motion lists
Motions can also be assembled into lists, and most of the operations
on them extend to lists. Let's create a list of motions: one for body
1 and one for body 2.
=#
ufcn(x,y,t) = 0.25*x*y*cos(t)
vfcn(x,y,t) = 0.25*(x^2-y^2)*cos(t)
m1 = DeformationMotion(ufcn,vfcn)
Ω = 1.0
kin = RotationalOscillation(Ω,π/4,0.0)
m2 = RigidBodyMotion(kin)
ml = MotionList([m1,m2]);

#=
Now let's see it in action
=#
@animate_motion bl ml π/100 4π (-5,5) (-5,5)

#=
We can also use the `surface_velocity!` function to get the velocities
of all surface points in the list. For example, to get them for
this previous list at time $t = 1.0$
=#
x, y = collect(bl)
u, v = zero(x), zero(y)
t = 1.0
surface_velocity!(u,v,bl,ml,t)

#=
We can determine the maximum velocity across the whole set of bodies:
=#
umax, i, tmax, bmax = maxlistvelocity(bl,ml)

#=
In this case, the maximum velocity occurs at t = 0 on body 2, index 301.
=#

#md # ## Body list functions
#md # ```@docs
#md # BodyList
#md # getrange
#md # Base.collect(::BodyList)
#md # Base.sum(::AbstractVector,::BodyList,::Int)
#md # Base.view(::AbstractVector,::BodyList,::Int)
#md # MotionList
#md # RigidTransformList
#md # Base.vec(::RigidTransformList)
#md # ```
