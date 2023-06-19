# # Creating and transforming bodies

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
The most basic functions of this package create an object of type
`Body`. There are a variety of such functions, a few of which we will demonstrate here.
Generally speaking, we are interesting in creating the object and placing it
in a certain position and orientation. We do this in two steps: we create the
basic shape, centered at the origin with a default orientation, and then
we transform the shape to a desired location and orientation using a `RigidTransform`.

It is useful to stress that each body stores two types of points internally:
the *endpoints* of the segments that comprise the body surface, and the
*midpoints* of these segments. The midpoints are intended for use in downstream
calculations, e.g. as forcing points in the calculations on immersed layers.
The midpoints are simply the geometric averages of the endpoints, so
endpoints are the ones that are transformed first, and midpoints are updated next.
=#


using RigidBodyTools
using Plots

#=
## Creating a shape
Let's first create a shape. For any shape, we have to make a choice of the
geometric dimensions (e.g, radius of the circle, side lengths of a rectangle),
as well as the points that we use to discretely represent the surface.
For this latter choice, there are two constructor types: we can specify the
number of points (as an integer), or we can specify the nominal spacing between
points (as a floating-point number).

The second approach is usually preferable when we use these tools for constructing immersed bodies.
It is important to stress that the algorithms for placing points attempt to make the spacing
as uniform as possible.

Let's create the most basic shape, a circle of radius 1. We will discretize
it with 100 points first:
=#
b = Circle(1.0,100)

#=
Now we will create the same body with a spacing of 0.02
=#
b = Circle(1.0,0.02)

#=
This choice led to 312 points along the circumference. Quick math will tell
you that the point spacing is probably not exactly 0.02. In fact, you can
find out the actual spacing with `dlengthmid`. This function calculates
the spacing associated with each point. (It does so by first calculating
the spacing between midpoints between each point and its two adjacent points.)
=#
dlengthmid(b)

#=
It is just a bit larger than 0.02.

A few other useful functions on the shape. To simply know the number of points,
=#
length(b)

#=
To find the outward normal vectors (based on the perpendicular to the line
joining the adjacent midpoints):
=#
nx, ny = normalmid(b)

#=
We can also plot the shape
=#
plot(b)

#=
Sometimes we don't want to fill in the shape (and maybe change the line color).
In that case, we can use
=#
plot(b,fill=:false,linecolor=:black)


#=
## Other shapes
Let's see some other shapes in action, like a square and an ellipse
=#
b1, b2 = Square(1.0,0.02), Ellipse(0.6,0.1,0.02)
plot(plot(b1), plot(b2))

#=
A NACA 4412 airfoil, with chord length 1, and 0.02 spacing between points.
=#
b = NACA4(0.04,0.4,0.12,0.02)
plot(b)

#=
A flat plate with no thickness
=#
b = Plate(1.0,0.02)
plot(b)


#=
and a flat plate with a 5 percent thickness (and rounded ends)
=#
b = ThickPlate(1.0,0.05,0.01)
plot(b)

#=
There are also some generic tools for creating shapes. A `BasicBody` simply
consists of points that describe the vertices. The interface for this is very simple.
=#
x = [1.0, 1.2, 0.7, 0.6, 0.2, -0.1, 0.1, 0.4]
y = [0.1, 0.5, 0.8, 1.2, 0.8, 0.6, 0.2, 0.3]
b = BasicBody(x,y)
#-
plot(b)
scatter!(b,markersize=3,markercolor=:black)

#=
However, this function does not insert any points along the sides between vertices.
We have to do the work of specifying these points in the original call. For this reason,
there are a few functions that are more directly useful. For example, we can
create a polygon from these vertices, with a specified spacing between points
distributed along the polygon sides
=#
b = Polygon(x,y,0.02)
#-
plot(b)
scatter!(b,markersize=3,markercolor=:black)

#=
Alternatively, we can interpret
those original points as control points for splines, with a spacing between points
along the splines provided:
=#
b = SplinedBody(x,y,0.02)
#-
plot(b)
scatter!(b,markersize=3,markercolor=:black)

#md # ## Body functions

#md # ```@docs
#md # BasicBody
#md # Polygon
#md # Circle
#md # Ellipse
#md # NACA4
#md # Plate
#md # Rectangle
#md # SplinedBody
#md # Square
#md # ```


#md # ## Shape utilities
#md # ```@docs
#md # centraldiff
#md # Base.diff(::Body)
#md # Base.diff(::BodyList)
#md # dlength
#md # dlengthmid
#md # Base.length(::Body)
#md # midpoints(::Body)
#md # midpoints(::BodyList)
#md # normal
#md # normalmid
#md # arccoord
#md # arccoordmid
#md # arclength(::Body)
#md # arclength(::BodyList)
#md # ```
