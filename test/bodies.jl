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



end
