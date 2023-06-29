# # Plucker vectors and coordinate transforms

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
Here we discuss the use of Plücker vectors and their transforms for
describing rigid-body motion and force. Plücker vectors succinctly describe
both the angular (rotational) and linear (translational) part of motion, and the angular (moment) and
linear (force) part of force. In three dimensions, a Plücker vector is 6-dimensional,
e.g., Plücker velocity and force vectors are

$$v = \begin{bmatrix} \Omega_x \\ \Omega_y \\ \Omega_z \\ U_x \\ U_y \\ U_z \end{bmatrix}, \qquad
f = \begin{bmatrix} M_x \\ M_y \\ M_z \\ F_x \\ F_y \\ F_z \end{bmatrix}$$

In two dimensions, there is only one angular component and two linear components, e.g.,

$$v = \begin{bmatrix} \Omega_z \\ U_x \\ U_y \end{bmatrix}, \qquad f = \begin{bmatrix} M_z \\ F_x \\ F_y \end{bmatrix}$$

We need to be able to transform these vectors from one coordinate system to another.
This requires rotating their components and shifting their center from one origin to another.
For example, a translational velocity based at system B will be different from
the translational velocity at system A because of the rotational velocity, $\Omega \times {}^Br_{A}$,
where ${}^Br_{A}$ is the vector from the origin of A to the origin of B.

Similarly, the moment about B will be different from the moment about A due to
the moment arm ${}^Br_{A} \times F$.
=#


using RigidBodyTools
using LinearAlgebra
using Plots

#=
## Plücker vectors
A Plücker vector is easily created by simply supplying a vector of its components
=#
v = PluckerMotion([1.0,2.0,3.0])

#=
This created a 2d motion vector, with angular velocity 1.0 and linear
velocity (2.0,3.0). One can also supply the angular and linear
parts separately, using keywords. If one of these keywords is
omitted, it defaults to zero for that part. Note that we also need
to write this as `PluckerMotion{2}` to specify the physical dimensionality.
For a 3d motion vector, one would write `PluckerMotion{3}` here.
=#
v2 = PluckerMotion{2}(angular=1.0,linear=[2.0,3.0])
v2 == v

#=
We can also pick off the angular and linear parts
=#
angular_only(v)
#=
and
=#
linear_only(v)

#=
Force vectors are similar
=#
f = PluckerForce([-1.0,-3.5,2.25])

#=
The vectors of the same type can be added and subtracted
=#
v3 = v + v2

#=
We can also take a scalar product of force and motion vectors
=#
dot(f,v)

#=
## Transforms
Transforms are constructed by describing the relationship between the two
coordinate systems. Consider the example in the figure below.

![CoordinateSystems.svg](./assets/CoordinateSystems.svg)

To develop the 2d transform from A to B, we supply the position $r$ and
the rotation angle $\theta$. For example, if B is shifted by [1,1]
and rotated by angle $\pi/6$ counterclockwise about A, then we construct the transform
as
=#
Xm = MotionTransform([1,1],π/6)

#=
Note that it uses the angle of rotation, $\pi/6$, to create a rotation
matrix operator.

A 2d force transform would be constructed by
=#
Xf = ForceTransform([1,1],π/6)

#=
For 3d transforms, we need to supply the rotation operator itself (as well
as the 3d translation vector). Often, this rotation is done by
rotating about a certain axis by a certain angle. We do this with the
`rotation_about_axis` function. For example, to rotate by $\pi/4$ about
an axis parallel to the vector $[1,1,1]$, then we use
=#
R = rotation_about_axis(π/4,[1,1,1])

#=
and then to translate this rotated system by $[-1,-2,-3]$,
=#
Xm = MotionTransform([-1,-2,-3],R)

#=
and similarly for a force transform.
=#

#=
We can also compute the inverses of these transforms, to transform back from
B to A
=#
inv(Xm)

#=
Transforms of the same type (motion or force) can be composed via multiplication to transform
from, e.g., A to B to C.
=#
Xm1 = MotionTransform([1.5,1.5],π/6)
Xm2 = MotionTransform([-1,1],π/3)
Xm2*Xm1

#=
## Transforming bodies
We can use motion transforms, in particular, to place bodies. We simply
apply the transform as a function, and it transforms the body's
coordinates. For example, transform `Xm1` above shifts the
body to `[1.5,1.5]` and rotates it counterclockwise by `π/6`:
=#
b = Ellipse(1.0,0.2,0.02)
plot(b,xlims=(-3,3),ylims=(-3,3),fillcolor=:gray)
plot!(Xm1(b),xlims=(-3,3),ylims=(-3,3))

#=
In the example above, we did not affect the original body by applying the
transform as a function. Rather, we created a copy of the body.

If, instead, you wish to transform the body in place, use `update_body!`
=#
update_body!(b,Xm1)

#=
One important note: a body stores a set of coordinates in its own intrinsic
coordinate system, and when a transform is applied to the body, it always
acts on these coordinates. This means that the transform's application on the body
cannot be carried out as a composite of operations, e.g. `T2(T1(b))` is not possible.
Insteady, in the application on the body, the transform is always interpreted such that system A
is the inertial coordinate system and B is the body system. Of course, the transform itself can always
be constructed from composite transforms.
=#

#=
Sometimes we need information about the normals in the body system.
For these, we can use `normalmid` with the flag `axes=:body`:
=#
nx, ny = normalmid(b,axes=:body)

#=
Finally, if you wish to transform the body's own coordinate system, rather
than use the transform to simply place the body in the inertial system, then
use `transform_body!`. This transforms the intrinsic coordinates of the body.
=#
transform_body!(b,Xm1)

#=
## Transforming Plücker vectors
Transforms can be applied to Plücker vectors to transform their components
between systems. Let's consider a 2d example in which the motion based at system A
is purely a rotation with angular velocity $\Omega = 1$, and we wish to transform this
to system B, translated by $[2,0]$ from A, but with axes aligned with B.
We expect that the velocity based at B should have the same angular velocity,
but also should have translational velocity equal to $[0,2]$ due to the angular
motion.

First we construct the motion vector at A
=#
Ω = 1.0
vA = PluckerMotion(Ω,[0,0])

#=
Now construct the transform from A to B:
=#
XA_to_B = MotionTransform([2,0],0)

#=
Now apply the transform to get the velocity at B:
=#
vB = XA_to_B*vA
#=
which gives the expected result. Now let's transform back, using the inverse,
and check that we get back to `vA`
=#
inv(XA_to_B)*vB

#md # ## Transform functions
#md # ```@docs
#md # PluckerMotion
#md # PluckerForce
#md # angular_only
#md # linear_only
#md # LinearAlgebra.dot(::PluckerForce,::PluckerMotion)
#md # MotionTransform
#md # ForceTransform
#md # Base.inv(::AbstractTransformOperator)
#md # Base.transpose(::AbstractTransformOperator)
#md # rotation_transform
#md # translation_transform
#md # update_body!
#md # transform_body!
#md # ```
