#=
Kinematics

Any set of kinematics added here must be accompanied by a function
(kin::AbstractKinematics)(t) that returns kinematic data of type `KinematicData`
holding t, c, ċ, c̈, α, α̇, α̈: the evaluation time, the reference
point position, velocity, acceleration (all complex), and angle,
angular velocity, and angular acceleration
=#
using DocStringExtensions
import ForwardDiff
import Base: +, *, -, >>, <<, show

using SpaceTimeFields
import SpaceTimeFields: Abstract1DProfile, >>, ConstantProfile, d_dt, FunctionProfile

## Individual DOF kinematics ##

abstract type AbstractDOFKinematics end

abstract type AbstractPrescribedDOFKinematics <: AbstractDOFKinematics end


struct DOFKinematicData
    t :: Float64
    x :: Float64
    ẋ :: Float64
    ẍ :: Float64
end

(k::AbstractPrescribedDOFKinematics)(t) = DOFKinematicData(t,k.x(t),k.ẋ(t),k.ẍ(t))
dof_state(kd::DOFKinematicData) = kd.x
dof_velocity(kd::DOFKinematicData) = kd.ẋ
dof_acceleration(kd::DOFKinematicData) = kd.ẍ


"""
    ConstantStateDOF(x0::Float64) <: AbstractPrescribedDOFKinematics

Set kinematics with constant state `x0`.
"""
struct ConstantStateDOF <: AbstractPrescribedDOFKinematics
    x0 :: Float64

    x :: Abstract1DProfile
    ẋ :: Abstract1DProfile
    ẍ :: Abstract1DProfile
end
function ConstantStateDOF(x0)
    x = ConstantProfile(x0)
    ẋ = d_dt(x)
    ẍ = d_dt(ẋ)
    ConstantStateDOF(x0,x,ẋ,ẍ)
end
function show(io::IO, p::ConstantStateDOF)
  print(io, "Constant state kinematics (state = $(p.x0))")
end

"""
    ConstantVelocityDOF(ẋ0::Float64) <: AbstractPrescribedDOFKinematics

Set kinematics with constant velocity `ẋ0`.
"""
struct ConstantVelocityDOF <: AbstractPrescribedDOFKinematics
    ẋ0 :: Float64

    x :: Abstract1DProfile
    ẋ :: Abstract1DProfile
    ẍ :: Abstract1DProfile
end
function ConstantVelocityDOF(ẋ0)
    f(t) = ẋ0*t
    x = FunctionProfile(f)
    ẋ = d_dt(x)
    ẍ = d_dt(ẋ)
    ConstantVelocityDOF(ẋ0,x,ẋ,ẍ)
end

function show(io::IO, p::ConstantVelocityDOF)
  print(io, "Constant velocity kinematics (velocity = $(p.ẋ0))")
end


"""
    SmoothRampDOF(x0,ẋ0,Δx,t0[;ramp=EldredgeRamp(11.0)]) <: AbstractPrescribedDOFKinematics

Kinematics describing a smooth ramp motion starting at time `t0` with nominal rate `ẋ0`.
The initial value is `x0`, and the ramp proceeds up to new value `x0 + Δx`.
The optional ramp argument is assumed to be
given by the smooth ramp `EldredgeRamp` with a smoothness factor of 11 (larger values
lead to sharper transitions on/off the ramp), but this
can be replaced by another Eldredge ramp with a different value or a `ColoniusRamp`.
"""
struct SmoothRampDOF <: AbstractPrescribedDOFKinematics
    x0 :: Float64
    ẋ0 :: Float64
    Δx :: Float64
    t0 :: Float64

    x :: Abstract1DProfile
    ẋ :: Abstract1DProfile
    ẍ :: Abstract1DProfile
end

function SmoothRampDOF(x0,ẋ0,Δx,t0; ramp = EldredgeRamp(11.0))
    Δt = Δx/ẋ0
    x = ConstantProfile(x0) + ẋ0*((ramp >> t0) - (ramp >> (t0 + Δt)))
    ẋ = d_dt(x)
    ẍ = d_dt(ẋ)
    SmoothRampDOF(x0,ẋ0,Δx,t0,x,ẋ,ẍ)
end

function show(io::IO, p::SmoothRampDOF)
  print(io, "Smooth ramp kinematics (initial state = $(p.x0), nominal rate = $(p.ẋ0), amplitude = $(p.Δx), nominal time = $(p.t0))")
end



"""
    OscillatoryDOF(amp,angfreq,phase,x0,ẋ0) <: AbstractPrescribedDOFKinematics

Set sinusoidal kinematics with amplitude `amp`, angular frequency `angfreq`,
phase `phase`, mean value `x0` (without `ẋ0`), and mean velocity `ẋ0`.
The function it provides is `x(t) = x0 + ẋ0*t + amp*sin(angfreq*t+phase)`.
"""
struct OscillatoryDOF <: AbstractPrescribedDOFKinematics
    "Amplitude"
    A :: Float64

    "Angular frequency"
    Ω :: Float64

    "Phase"
    ϕ :: Float64

    "Mean value"
    x0 :: Float64

    "Mean velocity"
    ẋ0 :: Float64

    x :: Abstract1DProfile
    ẋ :: Abstract1DProfile
    ẍ :: Abstract1DProfile
end

function OscillatoryDOF(A,Ω,ϕ,x0,ẋ0)
    f(t) = ẋ0*t
    x = ConstantProfile(x0) + FunctionProfile(f) + A*(Sinusoid(Ω) >> (ϕ/Ω))
    ẋ = d_dt(x)
    ẍ = d_dt(ẋ)
    OscillatoryDOF(A,Ω,ϕ,x0,ẋ0,x,ẋ,ẍ)
end

function show(io::IO, p::OscillatoryDOF)
  print(io, "Oscillatory kinematics (amplitude = $(p.A), ang freq = $(p.Ω), phase = $(p.ϕ), initial state = $(p.x0), mean velocity = $(p.ẋ0))")
end

"""
    CustomDOF(f::Function) <: AbstractPrescribedDOFKinematics

Set custom kinematics for a degree of freedom with a function `f`
that specifies its value at any given time.
"""
struct CustomDOF <: AbstractPrescribedDOFKinematics
    f :: Function

    x :: Abstract1DProfile
    ẋ :: Abstract1DProfile
    ẍ :: Abstract1DProfile
end
function CustomDOF(f::Function)
    x = FunctionProfile(f)
    ẋ = d_dt(x)
    ẍ = d_dt(ẋ)
    CustomDOF(f,x,ẋ,ẍ)
end

function show(io::IO, p::CustomDOF)
  print(io, "Custom kinematics")
end

"""
    ExogenousDOF() <: AbstractDOFKinematics

Sets a DOF as constrained, but with its behavior set by an exogenous process
at every instant. For such a DOF, one must provide a vector  `[x,ẋ,ẍ]`.
"""
struct ExogenousDOF <: AbstractDOFKinematics
end

function show(io::IO, p::ExogenousDOF)
  print(io, "Exogeneously-specified DOF")
end

"""
    UnconstrainedDOF([f::Function]) <: AbstractDOFKinematics

Sets a DOF as unconstrained, so that its behavior is either completely free
or determined by a given force response (e.g., spring and/or damper). This force
response is set by the optional input function `f`. The signature of `f` must be `f(x,xdot,t)`,
where `x` and `xdot` are the state of the dof and its derivative,
respectively, and `t` is the current time. It must return a single scalar,
serving as a force or torque for that DOF.
"""
struct UnconstrainedDOF <: AbstractDOFKinematics
  f :: Function
end
UnconstrainedDOF() = UnconstrainedDOF((x,xdot,t) -> zero(x))

function show(io::IO, p::UnconstrainedDOF)
  print(io, "Unconstrained DOF")
end


#####

# Joints




#####


"""
An abstract type for types that takes in time and returns `KinematicData(t,c, ċ, c̈, α, α̇, α̈)`.
It is important to note that `c` and `α` only provide the values relative to
their respective initial values at `t=0`.
"""
abstract type AbstractKinematics end

struct KinematicData
  t :: Float64
  c :: ComplexF64
  ċ :: ComplexF64
  c̈ :: ComplexF64
  α :: Float64
  α̇ :: Float64
  α̈ :: Float64
end

# APIs
"""
    complex_translational_position(k::KinematicData;[inertial=true]) -> ComplexF64

Return the complex translational position (relative to the initial
  position) of kinematic data `k`, expressed in inertial coordinates (if `inertial=true`,
    the default) or in comoving coordinates if `inertial=false`.
"""
complex_translational_position(k::KinematicData;inertial::Bool=true) =
        _complex_translation_position(k,Val(inertial))

_complex_translation_position(k,::Val{true}) = k.c
_complex_translation_position(k,::Val{false}) = k.c*exp(-im*k.α)


"""
    complex_translational_velocity(k::KinematicData;[inertial=true]) -> ComplexF64

Return the complex translational velocity of kinematic data `k`, relative
to inertial reference frame, expressed in inertial coordinates (if `inertial=true`,
  the default) or in comoving coordinates if `inertial=false`.
"""
complex_translational_velocity(k::KinematicData;inertial::Bool=true) =
        _complex_translation_velocity(k,Val(inertial))

_complex_translation_velocity(k,::Val{true}) = k.ċ
_complex_translation_velocity(k,::Val{false}) = k.ċ*exp(-im*k.α)


"""
    complex_translational_acceleration(k::KinematicData;[inertial=true]) -> ComplexF64

Return the complex translational acceleration of kinematic data `k`, relative
to inertial reference frame, expressed in inertial coordinates (if `inertial=true`,
  the default) or in comoving coordinates if `inertial=false`.
"""
complex_translational_acceleration(k::KinematicData;inertial::Bool=true) =
        _complex_translation_acceleration(k,Val(inertial))

_complex_translation_acceleration(k,::Val{true}) = k.c̈
_complex_translation_acceleration(k,::Val{false}) = k.c̈*exp(-im*k.α)

"""
    angular_position(k::KinematicData) -> Float64

Return the angular orientation of kinematic data `k`, relative to
inital value, in the inertial reference frame.
"""
angular_position(k::KinematicData) = k.α

"""
    angular_velocity(k::KinematicData) -> Float64

Return the angular velocity of kinematic data `k`, relative to
inertial reference frame.
"""
angular_velocity(k::KinematicData) = k.α̇

"""
    angular_acceleration(k::KinematicData) -> Float64

Return the angular acceleration of kinematic data `k`, relative to
inertial reference frame.
"""
angular_acceleration(k::KinematicData) = k.α̈

"""
    translational_position(k::KinematicData;[inertial=true]) -> Tuple

Return the translational position of kinematic data `k` (relative to the initial
  position), expressed as a Tuple in inertial coordinates (if `inertial=true`,
    the default) or in comoving coordinates if `inertial=false`.
"""
translational_position(k::KinematicData;kwargs...) =
          reim(complex_translational_position(k;kwargs...))
"""
    translational_velocity(k::KinematicData;[inertial=true]) -> Tuple

Return the translational velocity of kinematic data `k`, relative
to inertial reference frame, expressed as a Tuple in inertial coordinates (if `inertial=true`,
  the default) or in comoving coordinates if `inertial=false`.
"""
translational_velocity(k::KinematicData;kwargs...) =
          reim(complex_translational_velocity(k;kwargs...))

"""
    translational_acceleration(k::KinematicData;[inertial=true]) -> Tuple

Return the translational acceleration of kinematic data `k`, relative
to inertial reference frame, expressed as a Tuple in inertial coordinates (if `inertial=true`,
  the default) or in comoving coordinates if `inertial=false`.
"""
translational_acceleration(k::KinematicData;kwargs...) =
          reim(complex_translational_acceleration(k;kwargs...))


####

"""
    Kinematics(apk::AbstractDOFKinematics,xpk::AbstractDOFKinematics,ykp::AbstractDOFKinematics[;pivot=(0.0,0.0)])

Set the full 2-d kinematics of a rigid body, specifying the kinematics of the
the angle ``\\alpha`` (with `apk`) and ``x`` and ``y`` coordinates of the pivot point with `xpk` and `ypk`, respectively.
The pivot point is specified (in body-fixed coordinates)
with the optional argument `pivot`.
"""
struct Kinematics{AK<:AbstractDOFKinematics,XK<:AbstractDOFKinematics,YK<:AbstractDOFKinematics} <: AbstractKinematics
   x̃p :: Float64
   ỹp :: Float64
   apk :: AK
   xpk :: XK
   ypk :: YK
   Kinematics(x̃p::Real,ỹp::Real,apk,xpk,ypk) = new{typeof(apk),typeof(xpk),typeof(ypk)}(convert(Float64,x̃p),convert(Float64,ỹp),apk,xpk,ypk)
end

Kinematics(apk::AbstractDOFKinematics,xpk::AbstractDOFKinematics,ypk::AbstractDOFKinematics;pivot=(0.0,0.0)) =
                Kinematics(pivot...,apk,xpk,ypk)


function (kin::Kinematics)(t)
    z̃p = kin.x̃p + im*kin.ỹp
    apkd = kin.apk(t)
    xpkd = kin.xpk(t)
    ypkd = kin.ypk(t)

    zp = dof_state(xpkd) + im*dof_state(ypkd)
    żp = dof_velocity(xpkd) + im*dof_velocity(ypkd)
    z̈p = dof_acceleration(xpkd) + im*dof_acceleration(ypkd)
    α = dof_state(apkd)
    α̇ = dof_velocity(apkd)
    α̈ = dof_acceleration(apkd)

    dzp = z̃p*exp(im*α)
    c = zp - dzp
    ċ = żp - im*α̇*dzp
    c̈ = z̈p - (im*α̈ - α̇^2)*dzp
    KinematicData(t,c,ċ,c̈,α,α̇,α̈)
end

"""
    Constant(Up,Ω;[pivot = (0.0,0.0)])

Set constant translational Up and angular velocity Ω with respect to a specified point P.
This point is established by its initial position `pivot` (note that the initial angle
is assumed to be zero). By default, this initial position is (0,0).
`Up` and `pivot` can be specified by either Tuple or by complex value.
"""
Constant(żp::Complex, α̇;pivot=complex(0.0)) =
      Kinematics(ConstantVelocityDOF(α̇),
                 ConstantVelocityDOF(real(żp)),
                 ConstantVelocityDOF(imag(żp)); pivot=reim(pivot))

Constant(żp::Tuple,α̇;pivot=(0.0,0.0)) = Constant(complex(żp...),α̇;pivot=complex(pivot...))

"""
    Pitchup(U₀,a,K,α₀,t₀,Δα,ramp=EldredgeRamp(11.0)) <: AbstractKinematics

Kinematics describing a pitch-ramp motion (horizontal translation with rotation)
starting at time ``t_0`` about an axis at `a` (expressed relative to the centroid, in the ``\\tilde{x}``
  direction in the body-fixed coordinate system), with translational velocity `U₀`
in the inertial ``x`` direction, initial angle ``\\alpha_0``, dimensionless angular
velocity ``K = \\dot{\\alpha}_0c/2U_0``, and angular
change ``\\Delta\\alpha``. The optional ramp argument is assumed to be
given by the smooth ramp `EldredgeRamp` with a smoothness factor of 11 (larger values
lead to sharper transitions on/off the ramp), but this
can be replaced by another Eldredge ramp with a different value or a `ColoniusRamp`.
"""
Pitchup(U₀,a,K,α₀,t₀,Δα;ramp=EldredgeRamp(11.0)) =
            Kinematics(SmoothRampDOF(α₀,2K,Δα,t₀;ramp=ramp),
                       ConstantVelocityDOF(U₀),
                       ConstantStateDOF(0.0); pivot=(a,0.0))



"""
    Oscillation(Ux,Uy,α̇₀,ax,ay,Ω,Ax,Ay,ϕx,ϕy,α₀,Δα,ϕα) <: AbstractKinematics

Set 2-d oscillatory kinematics. This general constructor sets up motion of a
rotational axis, located at `ax`, `ay` (expressed relative to the body centroid, in
  a body-fixed coordinate system). The rotational axis motion is described by

``
x(t) = U_x t + A_x\\sin(\\Omega t - \\phi_x), \\quad y(t) = U_y t + A_y\\sin(\\Omega t - \\phi_y),
 \\quad \\alpha(t) = \\alpha_0 + \\dot{\\alpha}_0 t +
 \\Delta\\alpha \\sin(\\Omega t - \\phi_{\\alpha})
 ``

"""
Oscillation(Ux,Uy,α̇₀,ax,ay,Ω,Ax,Ay,ϕx,ϕy,α₀,Δα,ϕα) =
        Kinematics(OscillatoryDOF(Δα,Ω,ϕα,α₀,α̇₀),
                   OscillatoryDOF(Ax,Ω,ϕx,0.0,Ux),
                   OscillatoryDOF(Ay,Ω,ϕy,0.0,Uy); pivot=(ax,ay))



"""
    PitchHeave(U₀,a,Ω,α₀,Δα,ϕp,A,ϕh)

Create oscillatory pitching and heaving kinematics of a pitch axis at location ``a`` (expressed
relative to the centroid in the ``\\tilde{x}`` direction of the body-fixed coordinate system),
of the form (in inertial coordinates)

``
x(t) = U_0 t, \\quad y(t) = A\\sin(\\Omega t - \\phi_h),
 \\quad \\alpha(t) = \\alpha_0 + \\Delta\\alpha \\sin(\\Omega t - \\phi_p)
 ``
"""
PitchHeave(U₀, a, Ω, α₀, Δα, ϕp, A, ϕh) = Oscillation(U₀, 0.0, 0.0, a, 0.0, Ω, 0.0, A, 0.0, ϕh, α₀, Δα, ϕp)


"""
    OscillationXY(Ux,Uy,Ω,Ax,ϕx,Ay,ϕy)

Set oscillatory kinematics in the ``x`` and ``y`` directions, of the form

``
x(t) = U_x t + A_x \\sin(\\Omega t - \\phi_x), \\quad y(t) = U_y t + A_y \\sin(\\Omega t - \\phi_y)
``
"""
OscillationXY(Ux,Uy,Ω,Ax,ϕx,Ay,ϕy) = Oscillation(Ux, Uy, 0, 0, 0, Ω, Ax, Ay, ϕx, ϕy, 0, 0, 0)

"""
    OscillationX(Ux,Ω,Ax,ϕx)

Set oscillatory kinematics in the ``x`` direction, of the form

``
x(t) = U_x t + A_x \\sin(\\Omega t - \\phi_x),
``
"""
OscillationX(Ux,Ω,Ax,ϕx) = Oscillation(Ux, 0, 0, 0, 0, Ω, Ax, 0, ϕx, 0, 0, 0, 0)

"""
    OscillationY(Uy,Ω,Ay,ϕy)

Set oscillatory kinematics in the ``y`` direction, of the form

``
y(t) = U_y t + A_y \\sin(\\Omega t - \\phi_y)
``
"""
OscillationY(Uy,Ω,Ay,ϕy) = Oscillation(0, Uy, 0, 0, 0, Ω, 0, Ay, 0, ϕy, 0, 0, 0)

"""
    RotationalOscillation(ax,ay,Ω,α₀,α̇₀,Δα,ϕα)

Set oscillatory rotational kinematics about an axis located at `ax`, `ay`
(expressed relative to the body centroid, in a body-fixed coordinate system), of the form

``
\\alpha(t) = \\alpha_0 + \\dot{\\alpha}_0 t +
\\Delta\\alpha \\sin(\\Omega t - \\phi_{\\alpha})
``
"""
RotationalOscillation(ax,ay,Ω,α₀,α̇₀,Δα,ϕα) = Oscillation(0, 0, α̇₀, ax, ay, Ω, 0, 0, 0, 0, α₀, Δα, ϕα)

"""
    RotationalOscillation(Ω,Δα,ϕα)

Set oscillatory rotational kinematics about the centroid of the form

``
\\alpha(t) = \\Delta\\alpha \\sin(\\Omega t - \\phi_{\\alpha})
``
"""
RotationalOscillation(Ω,Δα,ϕα) = RotationalOscillation(0,0,Ω,0,0,Δα,ϕα)


####




####




abstract type Switch end
abstract type SwitchOn <: Switch end
abstract type SwitchOff <: Switch end

"""
    SwitchedKinematics <: AbstractKinematics

Modulates a given set of kinematics between simple on/off states. The velocity
specified by the given kinematics is toggled on/off.

# Fields
$(FIELDS)
"""
struct SwitchedKinematics{S <: Switch} <: AbstractKinematics

    "time at which the kinematics should be turned on"
    t_on :: Float64

    "time at which the kinematics should be turned off"
    t_off :: Float64

    "kinematics to be followed in the on state"
    kin :: AbstractKinematics

    off :: AbstractKinematics

    SwitchedKinematics(t_on,t_off,kin) = t_on > t_off ?
            new{SwitchOn}(t_on,t_off,kin,RigidBodyMotions.Constant(0,0)) :
            new{SwitchOff}(t_on,t_off,kin,RigidBodyMotions.Constant(0,0))
end

# note that these do not introduce impulsive changes into the derivatives
(p::SwitchedKinematics{SwitchOn})(t) = t <= p.t_on ? p.off(t) : p.kin(t-p.t_on)

(p::SwitchedKinematics{SwitchOff})(t) = t <= p.t_off ? p.kin(t-p.t_on) : p.off(t)
