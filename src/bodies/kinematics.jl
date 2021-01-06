#=
Kinematics
=#


struct Constant{C <: Complex, A <: Real} <: Kinematics
    ċ::C
    α̇::A
end
Constant(ċ, α̇) = Constant(complex(ċ), α̇)
(c::Constant{C})(t) where C = zero(C), c.ċ, zero(C), 0.0, c.α̇, 0.0
show(io::IO, c::Constant) = print(io, "Constant (ċ = $(c.ċ), α̇ = $(c.α̇))")

"""
    Pitchup <: Kinematics

Kinematics describing a pitchup motion (horizontal translation with rotation)

# Constructors
# Fields
$(FIELDS)
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

    α::Profile
    α̇::Profile
    α̈::Profile
end

function Pitchup(U₀, a, K, α₀, t₀, Δα, ramp)
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

    return c, ċ, c̈, α, α̇, α̈
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
 \\quad \\alpha(t) = \\alpha_0 + \\dot{\\alpha}_0 t + \\Delta\\alpha \\sin(\\Omega t - \\phi_{\\alpha})
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


    px::Profile
    ṗx::Profile
    p̈x::Profile

    py::Profile
    ṗy::Profile
    p̈y::Profile

    α::Profile
    α̇::Profile
    α̈::Profile
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

    return c, ċ, c̈, α, α̇, α̈
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
    PitchHeave <: Kinematics

Kinematics describing an oscillatory pitching and heaving (i.e. plunging) motion

# Constructors
# Fields
$(FIELDS)
"""
struct PitchHeave <: Kinematics
    "Freestream velocity"
    U₀::Float64

    "Axis of pitch rotation, relative to the plate centroid"
    a::Float64

    "Reduced frequency ``K = \\frac{\\Omega c}{2U_0}``"
    K::Float64

    "Phase of pitch (in radians)"
    ϕp::Float64

    "Phase of heave (in radians)"
    ϕh::Float64

    "Mean angle of attack"
    α₀::Float64

    "Amplitude of pitching"
    Δα::Float64

    "Amplitude of translational heaving"
    A::Float64

    Y::Profile
    Ẏ::Profile
    Ÿ::Profile

    α::Profile
    α̇::Profile
    α̈::Profile
end

function PitchHeave(U₀, a, K, ϕp, α₀, Δα, A, ϕh)
    p = A*(Sinusoid(2K) >> (ϕh/(2K)))
    ṗ = d_dt(p)
    p̈ = d_dt(ṗ)
    α = ConstantProfile(α₀) + Δα*(Sinusoid(2K) >> (ϕp/(2K)))
    α̇ = d_dt(α)
    α̈ = d_dt(α̇)
    PitchHeave(U₀, a, K, ϕp, ϕh, α₀, Δα, A, p, ṗ, p̈, α, α̇, α̈)
end

function (p::PitchHeave)(t)
    α = p.α(t)
    α̇ = p.α̇(t)
    α̈ = p.α̈(t)

    c = p.U₀*t + im*p.Y(t) - p.a*exp(im*α)
    ċ = p.U₀ + im*p.Ẏ(t) - p.a*im*α̇*exp(im*α)
    c̈ = im*p.Ÿ(t) + p.a*exp(im*α)*(α̇^2 - im*α̈)

    return c, ċ, c̈, α, α̇, α̈
end

function show(io::IO, p::PitchHeave)
    println(io, "Oscillatory pitch-heave kinematics with")
    println(io, "     Reduced frequency K = $(p.K)")
    println(io, "     Heaving amplitude A = $(p.A)")
    println(io, "     Pitching amplitude Δα = $(p.Δα)")
    println(io, "     Pitch lag ϕp = $(p.ϕp)")
    println(io, "     Heave lag ϕh = $(p.ϕh)")
end

"""
    OscillationXY(Ω,Ax,ϕx,Ay,ϕy) <: Kinematics

Set oscillatory kinematics in the ``x`` and ``y`` directions, of the form

``
x(t) = A_x \\sin(\\Omega t + \\phi_x), \\quad y(t) = A_y \\sin(\\Omega t + \\phi_y)
``
"""
struct OscillationXY <: Kinematics
    "Angular frequency"
    Ω :: Float64

    "Amplitude x direction"
    Ax:: Float64

    "Phase in x direction."
    ϕx :: Float64

    "Amplitude y direction"
    Ay:: Float64

    "Phase in y direction."
    ϕy :: Float64

    cx::Profile
    ċx::Profile
    c̈x::Profile

    cy::Profile
    ċy::Profile
    c̈y::Profile

end

function OscillationXY(Ω,Ax,ϕx,Ay,ϕy)
    Δtx = ϕx/Ω
    px = Ax*(Sinusoid(Ω) << Δtx)
    ṗx = d_dt(px)
    p̈x = d_dt(ṗx)

    Δty = ϕy/Ω
    py = Ay*(Sinusoid(Ω) << Δty)
    ṗy = d_dt(py)
    p̈y = d_dt(ṗy)
    OscillationXY(Ω, Ax, ϕx, Ay, ϕy, px, ṗx, p̈x, py, ṗy, p̈y)
end

function (p::OscillationXY)(t)
    α = 0.0
    α̇ = 0.0
    α̈ = 0.0

    c = ComplexF64(p.cx(t)) + im*ComplexF64(p.cy(t))
    ċ = ComplexF64(p.ċx(t)) + im*ComplexF64(p.ċy(t))
    c̈ = ComplexF64(p.c̈x(t)) + im*ComplexF64(p.c̈y(t))
    return c, ċ, c̈, α, α̇, α̈

    #return [p.ċ(t),0.0], [p.c̈(t),0.0], α̇
end
function show(io::IO, p::OscillationXY)
    println(io, "Oscillatory xy kinematics with")
    println(io, "     Frequency K = $(p.Ω)")
    println(io, "     x amplitude Ax = $(p.Ax)")
    println(io, "     x phase ϕx = $(p.ϕx)")
    println(io, "     y amplitude Ay = $(p.Ay)")
    println(io, "     y phase ϕy = = $(p.ϕy)")
end

struct OscilX <: Kinematics
    "Angular frequency"
    Ω :: Float64

    "Mean velocity"
    Umean :: Float64

    "Velocity amplitude"
    Ux:: Float64

    "Velocity phase"
    ϕx :: Float64

    cx::Profile
    ċx::Profile
    c̈x::Profile

end

function OscilX(Ω,Umean,Ux,ϕx)
    Δtx = ϕx/Ω
    px = ConstantProfile(0.0)
    ṗx = ConstantProfile(Umean) + Ux*(Sinusoid(Ω) << Δtx)
    p̈x = d_dt(ṗx)

    OscilX(Ω, Umean, Ux, ϕx, px, ṗx, p̈x)
end

function (p::OscilX)(t)
    α = 0.0
    α̇ = 0.0
    α̈ = 0.0

    c = ComplexF64(p.cx(t))
    ċ = ComplexF64(p.ċx(t))
    c̈ = ComplexF64(p.c̈x(t))
    return c, ċ, c̈, α, α̇, α̈

    #return [p.ċ(t),0.0], [p.c̈(t),0.0], α̇
end

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
