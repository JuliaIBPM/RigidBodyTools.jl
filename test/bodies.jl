using Statistics

@testset "Bodies" begin

  p = Plate(1,101)
  dx, dy = diff(p)
  @test maximum(dx) ≈ minimum(dx) ≈ 0.01
  @test maximum(dy) ≈ minimum(dy) ≈ 0.0

  dxc, dyc = centraldiff(p)
  @test dxc[1] ≈ 0.005 && dxc[101] ≈ 0.005 && maximum(dxc[2:100]) ≈ 0.01
  @test maximum(dyc) ≈ minimum(dyc) ≈ 0.0
  @test sum(dlength(p)) ≈ sum(dlengthmid(p))


  c = Rectangle(1,2,101)
  dx, dy = diff(c)
  @test length(dx) == 600
  @test sum(dlength(c)) ≈ 12.0

  c = Square(1,0.01)
  @test isapprox(mean(dlength(c)),0.01,atol=1e-4)

  c = Ellipse(1,2,0.01)
  @test isapprox(mean(dlength(c)),0.01,atol=1e-4)

  c = Ellipse(1,2,100)
  nx, ny = normalmid(c)
  @test nx[1] == ny[26] == -nx[51] == -ny[76] == 1.0
  @test abs(sum(nx)) < 1000.0*eps(1.0)
  @test abs(sum(ny)) < 1000.0*eps(1.0)


end

@testset "Lists" begin

  n1 = 101
  n2 = 201
  p = Plate(1,n1)
  c = Rectangle(1,2,n2)

  bl = BodyList([p,c])
  @test length(bl) == 2

  @test eltype(bl) == Body

  m1 = RigidBodyMotion(complex(0.0),0.0)
  m2 = RigidBodyMotion(RigidBodyTools.PitchHeave(1.0, 11.0, 0.2, 0.0, 0.0, 0.5, 1.0, 0.0))

  ml = RigidMotionList([m1,m2])
  @test length(ml) == 2
  @test eltype(ml) == RigidBodyMotion

  t1 = RigidTransform((0.0,0.0),0.0)
  t2 = RigidTransform((1.0,0.0),π/2)

  tl = RigidTransformList([t1,t2])
  push!(tl,t1)
  @test length(tl) == 3
  @test eltype(tl) == RigidTransform

  v = rand(numpts(bl))
  @test numpts(bl) == length(p)+length(c)
  @test v[1:length(p)] == view(v,bl,1)
  @test v[(length(p)+1):(length(p)+length(c))] == view(v,bl,2)

  tl = RigidTransformList([t1,t2])
  tl(bl)
  @test bl[1].cent == tl[1].trans
  @test bl[2].cent == tl[2].trans
  @test bl[1].α == tl[1].α
  @test bl[2].α == tl[2].α


end

@testset "Assign velocity" begin
    b = Circle(1.0,100)
    T = RigidTransform((rand(),rand()),0.0)
    T(b)

    ċ = rand(ComplexF64)
    m = RigidBodyMotion(ċ,1.0)
    u,v = assign_velocity(b,m,0.0)
    @test u[26] == -1.0+real(ċ)
    @test v[51] == -1.0+imag(ċ)

    u2, v2 = m(0.0,b)
    @test u == u2 && v == v2

    b2 = Circle(1.0,100)
    T2 = RigidTransform((rand(),rand()),0.0)
    T2(b2)
    ċ2 = rand(ComplexF64)
    m2 = RigidBodyMotion(ċ2,1.0)

    bl = BodyList([b,b2])
    ml = RigidMotionList([m,m2])
    u, v = assign_velocity(bl,ml,0.0)
    @test u[26] == -1.0+real(ċ)
    @test u[126] == -1.0+real(ċ2)
    @test v[151] == -1.0+imag(ċ2)

    u2, v2 = ml(0.0,bl)
    @test u2 == u && v2 == v

    vel = rigidbodyvelocity(ml,0.0)
    @test vel[3] == 1.0
    @test vel[1]+im*vel[2] == ċ
    @test vel[4]+im*vel[5] == ċ2
    @test vel[6] == 1.0

end