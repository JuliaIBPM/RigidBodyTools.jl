# # Deforming bodies

#md # ```@meta
#md # CurrentModule = RigidBodyTools
#md # ```

#=
Thus far we have only shown rigid body motion. However, we can
also prescribe surface deformation as an additional component
of a body's motion.
=#


using RigidBodyTools
using Plots

#=
Before we get started, let's define the same macro that we used earlier
in order to visualize our system's motion
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
For deforming bodies, we specify the velocity of the surface directly. This
deformation velocity is expressed in the coordinate system attached to the
body, rather than the inertial coordinate system. This enables the
motion to be easily superposed with the rigid-body motion described earlier.

It is also important to note that the motion is applied **to the endpoints**
of the surface segments. The midpoints are then constructed from the
updated endpoints.
=#
#=
Let's see an example. We will create an oscillatory deformation of a circle.
We create the motion by creating functions for each component of velocity.
=#
Ω = 2π
ufcn(x,y,t) = 0.25*x*y*Ω*cos(Ω*t)
vfcn(x,y,t) = 0.25*(x^2-y^2)*Ω*cos(Ω*t)
def = DeformationMotion(ufcn,vfcn)
#=
We will create a simple fixed revolute joint that anchors the body's center to the
inertial system.
=#
Xp_to_jp = MotionTransform(0.0,0.0,0.0)
Xc_to_jc = MotionTransform(0.0,0.0,0.0)
dofs = [ConstantVelocityDOF(0.0)]
joint = Joint(RevoluteJoint,0,Xp_to_jp,1,Xc_to_jc,dofs)

b = Circle(1.0,0.02)

ls = RigidBodyMotion(joint,b,def)

#=
Let's animate this motion
=#
# @animate_motion b ls 0.01 4 (-2,2) (-2,2)
