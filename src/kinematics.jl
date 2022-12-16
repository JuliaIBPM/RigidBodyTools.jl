#=
Kinematics

Any set of kinematics added here must be accompanied by a function
(kin::Kinematics)(t) that returns kinematic data of type `KinematicData`
holding t, c, ċ, c̈, α, α̇, α̈: the evaluation time, the reference
point position, velocity, acceleration (all complex), and angle,
angular velocity, and angular acceleration
=#
using DocStringExtensions
import ForwardDiff
import Base: +, *, -, >>, <<, show

using SpaceTimeFields
import SpaceTimeFields: Abstract1DProfile, >>, ConstantProfile, d_dt

"""
An abstract type for types that takes in time and returns `KinematicData(t,c, ċ, c̈, α, α̇, α̈)`.
"""
abstract type Kinematics end

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
    Constant(Up,Ω;[pivot = (0.0,0.0)])

Set constant translational Up and angular velocity Ω with respect to a specified point P.
This point is established by its initial position `pivot` (note that the initial angle
is assumed to be zero). By default, this initial position is (0,0).
`Up` and `pivot` can be specified by either Tuple or by complex value.
"""
struct Constant{C <: Complex, A <: Real} <: Kinematics
    żp::C
    α̇::A
    zp0::C
end
Constant(żp, α̇;pivot=complex(0.0)) = Constant(complex(żp...), α̇, complex(pivot...))
function (p::Constant{C})(t) where C
  α = p.α̇*t
  # z̃pc is position of pivot rel to centroid c, in comoving coordinates
  # (initial c0 and α0 are both assumed to be zero)
  z̃pc = p.zp0
  zp = p.zp0 + p.żp*t
  zpc = exp(im*α)*z̃pc
  c = zp - zpc
  ċ = p.żp - im*p.α̇*zpc
  c̈ = p.α̇^2*zpc
  KinematicData(t, c, ċ, c̈, α, p.α̇, 0.0)
end
function show(io::IO, p::Constant)
  println(io, "Constant rotation (żp = $(p.żp), α̇ = $(p.α̇)) ")
  print(io, "   about z̃p = $(p.zp0) (in comoving coordinates) ")
end

"""
    Pitchup(U₀,a,K,α₀,t₀,Δα,ramp=EldredgeRamp(11.0)) <: Kinematics

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
struct Pitchup <: Kinematics
    "Freestream velocity"
    U₀::Float64
    "Axis of rotation, relative to the plate centroid"
    a::Float64

    "Non-dimensional pitch rate ``K = \\dot{\\alpha}_0\\frac{c}{2U_0}``"
    K::Float64

    "Initial angle of attack"
    α₀::Float64
    "Nominal start of pitch up"
    t₀::Float64

    "Total pitching angle"
    Δα::Float64

    α::Abstract1DProfile
    α̇::Abstract1DProfile
    α̈::Abstract1DProfile
end

function Pitchup(U₀, a, K, α₀, t₀, Δα, ramp=EldredgeRamp(11.0))
    Δt = 0.5Δα/K
    p = ConstantProfile(α₀) + 2K*((ramp >> t₀) - (ramp >> (t₀ + Δt)))
    ṗ = d_dt(p)
    p̈ = d_dt(ṗ)

    Pitchup(U₀, a, K, α₀, t₀, Δα, p, ṗ, p̈)
end

function (p::Pitchup)(t)
    α = p.α(t)
    α̇ = p.α̇(t)
    α̈ = p.α̈(t)

    c = p.U₀*t - p.a*exp(im*α)
    ċ = p.U₀ - p.a*im*α̇*exp(im*α)
    if (t - p.t₀) > p.Δα/p.K
        c̈ = 0.0im
    else
        c̈ = p.a*exp(im*α)*(α̇^2 - im*α̈)
    end

    return KinematicData(t, c, ċ, c̈, α, α̇, α̈)
end

function show(io::IO, p::Pitchup)
    print(io, "Pitch-up kinematics with rate K = $(p.K)")
end

"""
    Oscillation(Ux,Uy,α̇₀,ax,ay,Ω,Ax,Ay,ϕx,ϕy,α₀,Δα,ϕα) <: Kinematics

Set 2-d oscillatory kinematics. This general constructor sets up motion of a
rotational axis, located at `ax`, `ay` (expressed relative to the body centroid, in
  a body-fixed coordinate system). The rotational axis motion is described by

``
x(t) = U_x t + A_x\\sin(\\Omega t - \\phi_x), \\quad y(t) = U_y t + A_y\\sin(\\Omega t - \\phi_y),
 \\quad \\alpha(t) = \\alpha_0 + \\dot{\\alpha}_0 t +
 \\Delta\\alpha \\sin(\\Omega t - \\phi_{\\alpha})
 ``

"""
struct Oscillation <: Kinematics
    "Translational velocity in x"
    Ux::Float64

    "Translational velocity in y"
    Uy::Float64

    "Constant rotation rate"
    α̇₀::Float64

    "Axis of pitch rotation, in body coordinates"
    ax::Float64
    ay::Float64

    "Angular frequency"
    Ω::Float64

    "Amplitude of oscillatory motion in x"
    Ax::Float64

    "Amplitude of oscillatory motion in y"
    Ay::Float64

    "Phase of x motion (in radians)"
    ϕx::Float64

    "Phase of y motion (in radians)"
    ϕy::Float64

    "Mean angle of attack"
    α₀::Float64

    "Amplitude of oscillatory rotation (in radians)"
    Δα::Float64

    "Phase of rotational motion (in radians)"
    ϕα::Float64


    px::Abstract1DProfile
    ṗx::Abstract1DProfile
    p̈x::Abstract1DProfile

    py::Abstract1DProfile
    ṗy::Abstract1DProfile
    p̈y::Abstract1DProfile

    α::Abstract1DProfile
    α̇::Abstract1DProfile
    α̈::Abstract1DProfile
end

function Oscillation(Ux, Uy, α̇₀, ax, ay, Ω, Ax, Ay, ϕx, ϕy, α₀, Δα, ϕa)
    Ω > 0.0 || error("Frequency must be positive and non-zero")
    px = Ax*(Sinusoid(Ω) >> (ϕx/Ω))
    ṗx = d_dt(px)
    p̈x = d_dt(ṗx)
    py = Ay*(Sinusoid(Ω) >> (ϕy/Ω))
    ṗy = d_dt(py)
    p̈y = d_dt(ṗy)
    α = ConstantProfile(α₀) + Δα*(Sinusoid(Ω) >> (ϕa/Ω))
    α̇ = d_dt(α)
    α̈ = d_dt(α̇)
    Oscillation(Ux, Uy, α̇₀, ax, ay, Ω, Ax, Ay, ϕx, ϕy, α₀, Δα, ϕa, px, ṗx, p̈x, py, ṗy, p̈y, α, α̇, α̈)
end

function (p::Oscillation)(t)
    α = p.α̇₀*t + p.α(t)
    α̇ = p.α̇₀   + p.α̇(t)
    α̈ = p.α̈(t)

    a = complex(p.ax,p.ay)
    U = complex(p.Ux,p.Uy)
    aeiα = a*exp(im*α)
    c = U*t + p.px(t) + im*p.py(t) - aeiα
    ċ = U + p.ṗx(t) + im*p.ṗy(t) - im*α̇*aeiα
    c̈ = p.p̈x(t) + im*p.p̈y(t) + aeiα*(α̇^2 - im*α̈)

    return KinematicData(t, c, ċ, c̈, α, α̇, α̈)
end

function show(io::IO, p::Oscillation)
    println(io, "Oscillatory kinematics with")
    println(io, "     Steady velocity U = ($(p.Ux),$(p.Uy))")
    println(io, "     Ref angle α₀ = $(p.α₀)")
    println(io, "     Mean rotation rate α̇₀ = $(p.α̇₀)")
    println(io, "     Pitch axis (rel. to centroid) a = ($(p.ax),$(p.ay))")
    println(io, "     Frequency Ω = $(p.Ω)")
    println(io, "     x amplitude Ax, phase lag ϕx = ($(p.Ax), $(p.ϕx))")
    println(io, "     y amplitude Ay, phase lag ϕy = ($(p.Ay), $(p.ϕy))")
    println(io, "     α amplitude Δα, phase lag ϕα = ($(p.Δα), $(p.ϕα))")
end


"""
    PitchHeave(U₀,a,Ω,α₀,Δα,ϕp,A,ϕh)

Create oscillatory pitching and heaving kinematics of a pitch axis at location ``a`` (expressed
relative to the centroid in the ``\\tilde{x}`` direction of the body-fixed coordinate system),
of the form (in inertial coordinates)

``
x(t) = U_0, \\quad y(t) = A\\sin(\\Omega t - \\phi_h),
 \\quad \\alpha(t) = \\alpha_0 + \\Delta\\alpha \\sin(\\Omega t - \\phi_p)
 ``
"""
PitchHeave(U₀, a, Ω, α₀, Δα, ϕp, A, ϕh) = Oscillation(U₀, 0, 0, a, 0, Ω, 0, A, 0, ϕh, α₀, Δα, ϕp)


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



abstract type Switch end
abstract type SwitchOn <: Switch end
abstract type SwitchOff <: Switch end

"""
    SwitchedKinematics <: Kinematics

Modulates a given set of kinematics between simple on/off states. The velocity
specified by the given kinematics is toggled on/off.

# Fields
$(FIELDS)
"""
struct SwitchedKinematics{S <: Switch} <: Kinematics

    "time at which the kinematics should be turned on"
    t_on :: Float64

    "time at which the kinematics should be turned off"
    t_off :: Float64

    "kinematics to be followed in the on state"
    kin :: Kinematics

    off :: Kinematics

    SwitchedKinematics(t_on,t_off,kin) = t_on > t_off ?
            new{SwitchOn}(t_on,t_off,kin,RigidBodyMotions.Constant(0,0)) :
            new{SwitchOff}(t_on,t_off,kin,RigidBodyMotions.Constant(0,0))
end

# note that these do not introduce impulsive changes into the derivatives
(p::SwitchedKinematics{SwitchOn})(t) = t <= p.t_on ? p.off(t) : p.kin(t-p.t_on)

(p::SwitchedKinematics{SwitchOff})(t) = t <= p.t_off ? p.kin(t-p.t_on) : p.off(t)
