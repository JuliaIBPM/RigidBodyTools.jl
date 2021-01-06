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



end
