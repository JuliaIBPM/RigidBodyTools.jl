Ux = rand()
Uy = rand()
α₀ = rand()
α̇₀ = rand()
ax = rand()
ay = rand()
Ω = 2π
Ax = rand()
ϕx = rand()
Ay = rand()
ϕy = rand()
Δα = rand()
ϕα = rand()

@testset "Motions" begin


  oscil = RigidBodyTools.Oscillation(Ux,Uy,α̇₀,ax,ay,Ω,Ax,Ay,ϕx,ϕy,α₀,Δα,ϕα)

  t = rand()
  c, ċ, c̈, α,α̇,α̈ = oscil(t)

  @test α ≈ α₀ + α̇₀*t + Δα*sin(Ω*t-ϕα)
  @test α̇ ≈ α̇₀ + Δα*Ω*cos(Ω*t-ϕα)
  @test α̈ ≈ -Δα*Ω^2*sin(Ω*t-ϕα)

  a = complex(ax,ay)
  @test real(c) ≈ Ux*t + Ax*sin(Ω*t-ϕx) - ax*cos(α) + ay*sin(α)
  @test imag(c) ≈ Uy*t + Ay*sin(Ω*t-ϕy) - ay*cos(α) - ax*sin(α)
  @test real(ċ) ≈ Ux + Ax*Ω*cos(Ω*t-ϕx) + α̇*(ax*sin(α) + ay*cos(α))
  @test imag(ċ) ≈ Uy + Ay*Ω*cos(Ω*t-ϕy) + α̇*(ay*sin(α) - ax*cos(α))
  @test real(c̈) ≈ -Ax*Ω^2*sin(Ω*t-ϕx) + α̈*(ax*sin(α) + ay*cos(α)) + α̇^2*(ax*cos(α) - ay*sin(α))
  @test imag(c̈) ≈ -Ay*Ω^2*sin(Ω*t-ϕy) + α̈*(ay*sin(α) - ax*cos(α)) + α̇^2*(ay*cos(α) + ax*sin(α))


  ph = RigidBodyTools.PitchHeave(Ux,ax,Ω,α₀,Δα,ϕα,Ay,ϕy)

  t = rand()
  c, ċ, c̈, α,α̇,α̈ = ph(t)

  @test α ≈ α₀ + Δα*sin(Ω*t-ϕα)
  @test α̇ ≈ Δα*Ω*cos(Ω*t-ϕα)
  @test α̈ ≈ -Δα*Ω^2*sin(Ω*t-ϕα)

  @test real(c) ≈ Ux*t - ax*cos(α)
  @test imag(c) ≈ Ay*sin(Ω*t-ϕy) - ax*sin(α)
  @test real(ċ) ≈ Ux + α̇*ax*sin(α)
  @test imag(ċ) ≈ Ay*Ω*cos(Ω*t-ϕy) - α̇*ax*cos(α)
  @test real(c̈) ≈ α̈*ax*sin(α) + α̇^2*ax*cos(α)
  @test imag(c̈) ≈ -Ay*Ω^2*sin(Ω*t-ϕy) - α̈*ax*cos(α) + α̇^2*ax*sin(α)

  ox = RigidBodyTools.OscillationX(Ux,Ω,Ax,ϕx)

  t = rand()
  c, ċ, c̈, α,α̇,α̈ = ox(t)
  @test α ≈ 0.0
  @test α̇ ≈ 0.0
  @test α̈ ≈ 0.0
  @test real(c) ≈ Ux*t + Ax*sin(Ω*t-ϕx)
  @test imag(c) ≈ 0.0
  @test real(ċ) ≈ Ux + Ax*Ω*cos(Ω*t-ϕx)
  @test imag(ċ) ≈ 0.0
  @test real(c̈) ≈ -Ax*Ω^2*sin(Ω*t-ϕx)
  @test imag(c̈) ≈ 0.0

  oy = RigidBodyTools.OscillationY(Uy,Ω,Ay,ϕy)

  t = rand()
  c, ċ, c̈, α,α̇,α̈ = oy(t)
  @test α ≈ 0.0
  @test α̇ ≈ 0.0
  @test α̈ ≈ 0.0
  @test real(c) ≈ 0.0
  @test imag(c) ≈ Uy*t + Ay*sin(Ω*t-ϕy)
  @test real(ċ) ≈ 0.0
  @test imag(ċ) ≈ Uy + Ay*Ω*cos(Ω*t-ϕy)
  @test real(c̈) ≈ 0.0
  @test imag(c̈) ≈ -Ay*Ω^2*sin(Ω*t-ϕy)

  ro = RigidBodyTools.RotationalOscillation(Ω,Δα,ϕα)

  t = rand()
  c, ċ, c̈, α,α̇,α̈ = ro(t)

  @test α ≈ Δα*sin(Ω*t-ϕα)
  @test α̇ ≈ Δα*Ω*cos(Ω*t-ϕα)
  @test α̈ ≈ -Δα*Ω^2*sin(Ω*t-ϕα)
  @test c ≈ 0.0
  @test ċ ≈ 0.0
  @test c̈ ≈ 0.0



end
