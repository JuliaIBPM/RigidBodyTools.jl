using Statistics

const MYEPS = 20*eps()

@testset "Bodies" begin

  bn = nothing
  @test numpts(bn) == 0

  p = Plate(1,100)
  dx, dy = diff(p)
  @test maximum(dx) ≈ minimum(dx) ≈ 0.01
  @test maximum(dy) ≈ minimum(dy) ≈ 0.0

  @test sum(dlength(p)) ≈ sum(dlengthmid(p))


  c = Rectangle(1,2,400)
  dx, dy = diff(c)
  @test length(dx) == 400
  @test sum(dlength(c)) ≈ 12.0

  s = arccoord(c)
  @test s[end] ≈ 12.0 ≈ arclength(c)
  smid = arccoordmid(c)
  @test smid[1] > 0.0 && smid[end] < s[end]


  c = Square(1,0.01)
  @test isapprox(mean(dlength(c)),0.01,atol=1e-4)

  c = Ellipse(1,2,0.01)
  @test isapprox(mean(dlength(c)),0.01,atol=1e-4)

  c = Ellipse(1,2,100)
  nx, ny = normalmid(c)
  #@test nx[1] == ny[26] == -nx[51] == -ny[76] == 1.0
  @test abs(sum(nx)) < 1000.0*eps(1.0)
  @test abs(sum(ny)) < 1000.0*eps(1.0)


  T = MotionTransform((1.0,-1.0),π/4)
  update_body!(c,T)
  nx2, ny2 = normalmid(c,ref=true)
  @test sum(abs.(nx .- nx2)) < 1000.0*eps(1.0)
  @test sum(abs.(ny .- ny2)) < 1000.0*eps(1.0)


end

@testset "Rectangle" begin

  b = Rectangle(0.5,1.0,60)
  nx0, ny0 = normalmid(b)

  @test all(isapprox.(nx0[1:10],0.0,atol=MYEPS)) &&
        all(isapprox.(nx0[11:30],1.0,atol=MYEPS)) &&
        all(isapprox.(nx0[31:40],0.0,atol=MYEPS)) &&
        all(isapprox.(nx0[41:60],-1.0,atol=MYEPS))
  @test all(isapprox.(ny0[1:10],-1.0,atol=MYEPS)) &&
        all(isapprox.(ny0[11:30],0.0,atol=MYEPS)) &&
        all(isapprox.(ny0[31:40],1.0,atol=MYEPS)) &&
        all(isapprox.(ny0[41:60],0.0,atol=MYEPS))

  ds = dlengthmid(b)
  @test all(ds .≈ 0.1)

  T = MotionTransform((1.0,-1.0),π/2)
  update_body!(b,T)

  T2 = MotionTransform((1,-1),π/2)
  update_body!(b,T2)

  nx, ny = normalmid(b)
  @test all(isapprox.(nx[1:10],1.0,atol=MYEPS)) &&
        all(isapprox.(nx[11:30],0.0,atol=MYEPS)) &&
        all(isapprox.(nx[31:40],-1.0,atol=MYEPS)) &&
        all(isapprox.(nx[41:60],0.0,atol=MYEPS))
  @test all(isapprox.(ny[1:10],0.0,atol=MYEPS)) &&
        all(isapprox.(ny[11:30],1.0,atol=MYEPS)) &&
        all(isapprox.(ny[31:40],0.0,atol=MYEPS)) &&
        all(isapprox.(ny[41:60],-1.0,atol=MYEPS))

  ds = dlengthmid(b)
  @test all(ds .≈ 0.1)

  T = MotionTransform((1.2,-1.5),π/4)
  update_body!(b,T)

  nx2, ny2 = normalmid(b,ref=true)
  @test sum(abs.(nx0 .- nx2)) < 1000.0*eps(1.0)
  @test sum(abs.(ny0 .- ny2)) < 1000.0*eps(1.0)


end

@testset "Polygons" begin

  b = Polygon([1.0,1.0,0.0,0.0],[0.0,0.5,0.5,0.0],0.02)
  nx0, ny0 = normalmid(b)

  @test b.side[1] == 1:25 && b.side[2] == 26:75 && b.side[3] == 76:100 && b.side[4] == 101:150


  @test all(isapprox.(nx0[b.side[1]],1.0,atol=MYEPS)) &&
        all(isapprox.(nx0[b.side[2]],0.0,atol=MYEPS)) &&
        all(isapprox.(nx0[b.side[3]],-1.0,atol=MYEPS)) &&
        all(isapprox.(nx0[b.side[4]],0.0,atol=MYEPS))
  @test all(isapprox.(ny0[b.side[1]],0.0,atol=MYEPS)) &&
        all(isapprox.(ny0[b.side[2]],1.0,atol=MYEPS)) &&
        all(isapprox.(ny0[b.side[3]],0.0,atol=MYEPS)) &&
        all(isapprox.(ny0[b.side[4]],-1.0,atol=MYEPS))

  ds = dlengthmid(b)
  @test all(ds .≈ 0.02)


end

@testset "Lists" begin

  bl = BodyList()
  @test numpts(bl) == 0

  n1 = 101
  n2 = 201
  p = Plate(1,n1)
  c = Rectangle(1,2,n2)

  bl = BodyList([p,c])
  @test length(bl) == 2

  @test eltype(bl) == Body

  sl = arccoord(bl)
  @test sl[1:numpts(p)] == arccoord(p)
  @test sl[numpts(p)+1:numpts(bl)] == arccoord(c)


  t1 = MotionTransform((0.0,0.0),0.0)
  t2 = MotionTransform((1.0,0.0),π/2)

  tl = MotionTransformList([t1,t2])
  push!(tl,t1)
  @test length(tl) == 3
  @test eltype(tl) == MotionTransform

  v = rand(numpts(bl))
  @test numpts(bl) == numpts(p)+numpts(c)
  @test v[1:numpts(p)] == view(v,bl,1)
  @test v[(numpts(p)+1):(numpts(p)+numpts(c))] == view(v,bl,2)

  tl = MotionTransformList([t1,t2])
  update_body!(bl,tl)

  @test [bl[1].cent...] == translation(tl[1])
  @test [bl[2].cent...] == translation(tl[2])
  @test bl[1].α == RigidBodyTools._get_angle_of_2d_transform(tl[1])
  @test bl[2].α == RigidBodyTools._get_angle_of_2d_transform(tl[2])

  b1 = Rectangle(0.5,1.0,60)
  b2 = Rectangle(0.5,1.0,60)

  bl = BodyList()
  push!(bl,b1)
  push!(bl,b2)
  nx, ny = normalmid(bl)
  nxb1, nyb1 = normalmid(b1)
  nxb2, nyb2 = normalmid(b2)
  @test nx[1:60] == nxb1 && ny[1:60] == nyb1 && nx[61:120] == nxb2 && ny[61:120] == nyb2


end
