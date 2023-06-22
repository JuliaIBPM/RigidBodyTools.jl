using LinearAlgebra

@testset "Rotations" begin

  θ = rand()
  r = [0,0,1]
  R1 = rotation_about_z(θ)
  R2 = rotation_about_axis(θ,r)
  @test R1 ≈ R2

  r = [0,1,0]
  R1 = rotation_about_y(θ)
  R2 = rotation_about_axis(θ,r)
  @test R1 ≈ R2

  r = [1,0,0]
  R1 = rotation_about_x(θ)
  R2 = rotation_about_axis(θ,r)
  @test R1 ≈ R2


end

@testset "Transforms" begin

  x = rand(3)
  @test cross_vector(cross_matrix(x)) ≈ x

  r = rand(2)
  Θ = 0.0
  TM = MotionTransform(r,Θ)
  vA = rand(3)

  vB = TM*vA

  # Should give velocity based at origin of B
  @test vB ≈ vA .+ [0.0,-vA[1]*r[2],vA[1]*r[1]]

  TF = ForceTransform(r,Θ)
  fA = rand(3)
  fB = TF*fA

  # Should give moment about origin of B
  @test fB ≈ fA .- [r[1]*fA[3]-r[2]*fA[2],0.0,0.0]

  x2 = rand(2)
  Θ = rand()
  R2 = [cos(-Θ) -sin(-Θ); sin(-Θ) cos(-Θ)]
  TM2_1 = MotionTransform(x2,R2)
  TM2_2 = MotionTransform(x2,Θ)
  @test TM2_1 isa MotionTransform{2}
  @test TM2_1.matrix ≈ TM2_2.matrix

  TM3_1 = MotionTransform(rand(3),rotation_about_axis(rand(),rand(3)))
  TM3_2 = MotionTransform(rand(3),rotation_about_axis(rand(),rand(3)))
  @test TM3_1 isa MotionTransform{3}

  TM3_12 = TM3_1*TM3_2
  @test TM3_12 isa MotionTransform
  @test MotionTransform(TM3_12.x,TM3_12.R).matrix ≈ TM3_12.matrix

  TM3_t = transpose(TM3_1)
  @test ForceTransform(TM3_t.x,TM3_t.R).matrix ≈ TM3_t.matrix

  TM2_1 = MotionTransform(rand(2),rotation_about_z(rand()))
  TM2_2 = MotionTransform(rand(2),rotation_about_z(rand()))
  @test TM2_1 isa MotionTransform{2}

  TM2_12 = TM2_1*TM2_2
  @test TM2_12 isa MotionTransform
  @test MotionTransform{2}(TM2_12.x,TM2_12.R).matrix ≈ TM2_12.matrix

  TM2_t = transpose(TM2_1)
  @test ForceTransform{2}(TM2_t.x,TM2_t.R).matrix ≈ TM2_t.matrix

  TMinv = inv(TM)
  TMiTM = TMinv*TM
  @test isapprox(TMiTM.x,zeros(3),atol=1e-15)
  @test isapprox(TMiTM.R,I,atol=1e-15)
  @test isapprox(TMiTM.matrix,I,atol=1e-15)

  TMiTM = TM*TMinv
  @test isapprox(TMiTM.x,zeros(3),atol=1e-15)
  @test isapprox(TMiTM.R,I,atol=1e-15)
  @test isapprox(TMiTM.matrix,I,atol=1e-15)

  x = (rand(),rand())
  θ = rand()
  TM = MotionTransform(x,θ)
  TR = RigidTransform(x,θ)

  x̃ = (rand(),rand())
  @test TM(x̃...) == TR(x̃...)

end

@testset "Plucker vectors" begin

  vA = PluckerMotion(rand(3))

  vAr = angular_only(vA)
  @test vA.angular == vAr.angular
  @test all(vAr.linear .== 0.0)

  vAl = linear_only(vA)
  @test vA.linear == vAl.linear
  @test all(vAl.angular .== 0.0)

  vA = PluckerMotion(rand(6))

  vAr = angular_only(vA)
  @test vA.angular == vAr.angular
  @test all(vAr.linear .== 0.0)

  vAl = linear_only(vA)
  @test vA.linear == vAl.linear
  @test all(vAl.angular .== 0.0)

  @test vAr + vAl == vA

  x = (rand(),rand())
  θ = rand()
  TM = MotionTransform(x,θ)
  vA = PluckerMotion(rand(3))
  TM*vA

  vA = PluckerForce(rand(3))

  vAr = angular_only(vA)
  @test vA.angular == vAr.angular
  @test all(vAr.linear .== 0.0)

  vAl = linear_only(vA)
  @test vA.linear == vAl.linear
  @test all(vAl.angular .== 0.0)

  vA = PluckerForce(rand(6))

  vAr = angular_only(vA)
  @test vA.angular == vAr.angular
  @test all(vAr.linear .== 0.0)

  vAl = linear_only(vA)
  @test vA.linear == vAl.linear
  @test all(vAl.angular .== 0.0)

  @test vAr + vAl == vA

  x = (rand(),rand())
  θ = rand()
  TM = ForceTransform(x,θ)
  vA = PluckerForce(rand(3))
  TM*vA

end
