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
  #c, ċ, c̈, α,α̇,α̈ = oscil(t)
  k = oscil(t)


  @test angular_position(k) ≈ α₀ + α̇₀*t + Δα*sin(Ω*t-ϕα)
  @test angular_velocity(k) ≈ α̇₀ + Δα*Ω*cos(Ω*t-ϕα)
  @test angular_acceleration(k) ≈ -Δα*Ω^2*sin(Ω*t-ϕα)

  a = complex(ax,ay)
  α = angular_position(k)
  α̇ = angular_velocity(k)
  α̈ = angular_acceleration(k)
  c = complex_translational_position(k)
  ċ = complex_translational_velocity(k)
  c̈ = complex_translational_acceleration(k)
  @test real(c) ≈ Ux*t + Ax*sin(Ω*t-ϕx) - ax*cos(α) + ay*sin(α)
  @test imag(c) ≈ Uy*t + Ay*sin(Ω*t-ϕy) - ay*cos(α) - ax*sin(α)
  @test real(ċ) ≈ Ux + Ax*Ω*cos(Ω*t-ϕx) + α̇*(ax*sin(α) + ay*cos(α))
  @test imag(ċ) ≈ Uy + Ay*Ω*cos(Ω*t-ϕy) + α̇*(ay*sin(α) - ax*cos(α))
  @test real(c̈) ≈ -Ax*Ω^2*sin(Ω*t-ϕx) + α̈*(ax*sin(α) + ay*cos(α)) + α̇^2*(ax*cos(α) - ay*sin(α))
  @test imag(c̈) ≈ -Ay*Ω^2*sin(Ω*t-ϕy) + α̈*(ay*sin(α) - ax*cos(α)) + α̇^2*(ay*cos(α) + ax*sin(α))
  @test isa(translational_position(k),Tuple)
  @test isa(translational_velocity(k),Tuple)
  @test isa(translational_acceleration(k),Tuple)
  @test translational_position(k) == reim(c)
  @test translational_velocity(k) == reim(ċ)
  @test translational_acceleration(k) == reim(c̈)


  ph = RigidBodyTools.PitchHeave(Ux,ax,Ω,α₀,Δα,ϕα,Ay,ϕy)

  t = rand()
  #c, ċ, c̈, α,α̇,α̈ = ph(t)
  k = ph(t)

  α = angular_position(k)
  α̇ = angular_velocity(k)
  α̈ = angular_acceleration(k)
  c = complex_translational_position(k)
  ċ = complex_translational_velocity(k)
  c̈ = complex_translational_acceleration(k)

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
  #c, ċ, c̈, α,α̇,α̈ = ox(t)
  k = ox(t)

  @test angular_position(k) ≈ 0.0
  @test angular_velocity(k) ≈ 0.0
  @test angular_acceleration(k) ≈ 0.0
  xc, yc = translational_position(k)
  uc, vc = translational_velocity(k)
  axc, ayc = translational_acceleration(k)
  @test xc ≈ Ux*t + Ax*sin(Ω*t-ϕx)
  @test yc ≈ 0.0
  @test uc ≈ Ux + Ax*Ω*cos(Ω*t-ϕx)
  @test vc ≈ 0.0
  @test axc ≈ -Ax*Ω^2*sin(Ω*t-ϕx)
  @test ayc ≈ 0.0

  oy = RigidBodyTools.OscillationY(Uy,Ω,Ay,ϕy)

  t = rand()
  #c, ċ, c̈, α,α̇,α̈ = oy(t)
  k = oy(t)

  @test angular_position(k) ≈ 0.0
  @test angular_velocity(k) ≈ 0.0
  @test angular_acceleration(k) ≈ 0.0
  xc, yc = translational_position(k)
  uc, vc = translational_velocity(k)
  axc, ayc = translational_acceleration(k)
  @test xc ≈ 0.0
  @test yc ≈ Uy*t + Ay*sin(Ω*t-ϕy)
  @test uc ≈ 0.0
  @test vc ≈ Uy + Ay*Ω*cos(Ω*t-ϕy)
  @test axc ≈ 0.0
  @test ayc ≈ -Ay*Ω^2*sin(Ω*t-ϕy)

  ro = RigidBodyTools.RotationalOscillation(Ω,Δα,ϕα)

  t = rand()
  #c, ċ, c̈, α,α̇,α̈ = ro(t)
  k = ro(t)

  @test angular_position(k) ≈ Δα*sin(Ω*t-ϕα)
  @test angular_velocity(k) ≈ Δα*Ω*cos(Ω*t-ϕα)
  @test angular_acceleration(k) ≈ -Δα*Ω^2*sin(Ω*t-ϕα)
  @test complex_translational_position(k) ≈ 0.0
  @test complex_translational_velocity(k) ≈ 0.0
  @test complex_translational_acceleration(k) ≈ 0.0

  U₀ = 0.0
  a = 0.5
  K = 0.2
  t₀ = 0.5
  pu = Pitchup(U₀,a,K,α₀,t₀,Δα,EldredgeRamp(20.0))
  t = rand()
  k = pu(t)
  α̇ = angular_velocity(k)
  α̈ = angular_acceleration(k)

  @test complex_translational_position(k,inertial=false) ≈ -a
  @test complex_translational_velocity(k,inertial=false) ≈ -a*im*α̇
  @test complex_translational_acceleration(k,inertial=false) ≈ a*(α̇^2 - im*α̈)

end

@testset "Direct motions" begin

  u, v = rand(5), rand(5)
  m = RigidBodyTools.ConstantDeformationMotion(u,v)

  ml = RigidBodyTools.MotionList([m])

  push!(ml,m)


end

@testset "Motion velocities and states" begin

  b1 = Circle(1.0,100)
  b2 = Rectangle(2.0,0.5,400)
  bl = BodyList([b1,b2])

  T1 = RigidTransform((rand(),rand()),rand())
  T2 = RigidTransform((rand(),rand()),rand())
  tl = RigidTransformList([T1,T2])

  tl(bl)

  kin = Pitchup(1.0,0.5,0.2,0.0,0.5,π/4,EldredgeRamp(20.0))
  m1 = RigidBodyMotion(kin)
  m2 = ConstantDeformationMotion(one.(b2.x),zero.(b2.y))

  ml = MotionList([m1,m2])

  x0 = motion_state(bl,ml)

  @test x0[1:3] == vec(T1)[1:3]
  @test x0[4:end] == vcat(b2.x̃end,b2.ỹend)

  t = rand()
  u = motion_velocity(bl,ml,t)
  @test u[1:3] == motion_velocity(b1,m1,t)
  @test u[4:end] == motion_velocity(b2,m2,t)


  bl2 = deepcopy(bl)
  update_body!(bl2,x0,ml) # This should not change the bodies
  for i in 1:length(bl)
    @test bl2[i].cent == bl[i].cent && bl2[i].α == bl[i].α
    @test bl2[i].x̃ == bl[i].x̃ && bl2[i].ỹ == bl[i].ỹ
    @test bl2[i].x == bl[i].x && bl2[i].y == bl[i].y
  end

  u, v = zero(b1.x),zero(b1.y)
  surface_velocity!(u,v,b1,m1,0.0)
  maxvelocity(b1,m1)

  m = RigidAndDeformingMotion(m1,m2)
  x0 = motion_state(b2,m)

  @test x0[1:3] == vec(T2)[1:3]
  @test x0[4:end] == vcat(b2.x̃end,b2.ỹend)

  t = rand()
  u = motion_velocity(b2,m,t)
  @test u[1:3] == motion_velocity(b2,m1,t)
  @test u[4:end] == motion_velocity(b2,m2,t)



end
