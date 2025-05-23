using LinearAlgebra
using JLD2

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
  @test RigidBodyTools._get_angular_part(vA) === RigidBodyTools._get_angular_part(vAr)

  vAl = linear_only(vA)
  @test RigidBodyTools._get_linear_part(vA) === RigidBodyTools._get_linear_part(vAl)

  vA = PluckerMotion(rand(6))

  vAr = angular_only(vA)
  @test RigidBodyTools._get_angular_part(vA) === RigidBodyTools._get_angular_part(vAr)

  vAl = linear_only(vA)
  @test RigidBodyTools._get_linear_part(vA) === RigidBodyTools._get_linear_part(vAl)

  x = (rand(),rand())
  θ = rand()
  TM = MotionTransform(x,θ)
  vA = PluckerMotion(rand(3))
  TM*vA

  vA = PluckerForce(rand(3))

  vAr = angular_only(vA)
  @test RigidBodyTools._get_angular_part(vA) === RigidBodyTools._get_angular_part(vAr)

  vAl = linear_only(vA)
  @test RigidBodyTools._get_linear_part(vA) === RigidBodyTools._get_linear_part(vAl)

  vA = PluckerForce(rand(6))

  vAr = angular_only(vA)
  @test RigidBodyTools._get_angular_part(vA) === RigidBodyTools._get_angular_part(vAr)

  vAl = linear_only(vA)
  @test RigidBodyTools._get_linear_part(vA) === RigidBodyTools._get_linear_part(vAl)

  x = (rand(),rand())
  θ = rand()
  TM = ForceTransform(x,θ)
  vA = PluckerForce(rand(3))
  TM*vA

  TM*angular_only(vA)

  TM*linear_only(vA)


  TM = MotionTransform(x,θ)
  vA = PluckerMotion{2}(angular=1)
  vB = PluckerMotion([1,2,3])
  @test TM*angular_only(vA) == TM*angular_only(vB)

  vA = PluckerMotion{2}(linear=[2,3])
  @test TM*linear_only(vA) == TM*linear_only(vB)

  TM = MotionTransform(rand(3),rotation_about_axis(rand(),rand(3)))
  vA = PluckerMotion{3}(angular=[1,2,3])
  vB = PluckerMotion([1,2,3,4,5,6])
  @test TM*angular_only(vA) == TM*angular_only(vB)

  vA = PluckerMotion{3}(linear=[4,5,6])
  @test TM*linear_only(vA) == TM*linear_only(vB)

  vA = PluckerMotion(rand(3))
  fA = PluckerForce(rand(3))

  @test dot(angular_only(fA),vA) == dot(fA,angular_only(vA)) == dot(angular_only(vA),fA) == dot(vA,angular_only(fA))
  @test dot(linear_only(fA),vA) == dot(fA,linear_only(vA)) == dot(linear_only(vA),fA) == dot(vA,linear_only(fA))

  vA = PluckerMotion(rand(6))
  fA = PluckerForce(rand(6))

  @test dot(angular_only(fA),vA) == dot(fA,angular_only(vA)) == dot(angular_only(vA),fA) == dot(vA,angular_only(fA))
  @test dot(linear_only(fA),vA) == dot(fA,linear_only(vA)) == dot(linear_only(vA),fA) == dot(vA,linear_only(fA))


end

@testset "Linked systems" begin
  Xp_to_j1 = MotionTransform(0.0,0.0,0.0)
  Xch_to_j1 = MotionTransform(-0.5,0.0,0.0)

  Xp_to_j2 = MotionTransform(1.02,0.0,0.0)
  Xch_to_j2 = MotionTransform(-1.02,0.0,0.0)

  Xp_to_j3 = MotionTransform(-5.0,0.0,0.0)
  Xch_to_j3 = MotionTransform(-0.5,0.0,0.0)

  dofs1 = [OscillatoryDOF(π/4,2π,0.0,0.0),ConstantPositionDOF(0.0),OscillatoryDOF(1.0,2π,-π/2,0.0)]
  #dofs11 = [OscillatoryDOF(π/4,2π,0.0,0.0),UnconstrainedDOF(),ExogenousDOF()]
  #dofs12 = [OscillatoryDOF(π/4,2π,π/3,0.0),ConstantPositionDOF(0.0),ExogenousDOF()]


  dofs2 = [OscillatoryDOF(π/4,2π,π/4,0.0)]
  dofs3 = [OscillatoryDOF(π/4,2π,π/2,0.0)]
  dofs4 = [ConstantVelocityDOF(0)]


  joint1 = Joint(FreeJoint2d,0,Xp_to_j1,1,Xch_to_j1,dofs1)
  joint2 = Joint(RevoluteJoint,1,Xp_to_j2,2,Xch_to_j2,dofs2)
  joint3 = Joint(RevoluteJoint,2,Xp_to_j2,3,Xch_to_j2,dofs3)
  joint4 = Joint(RevoluteJoint,2,Xp_to_j2,4,Xch_to_j2,dofs4)
  joint5 = Joint(FreeJoint2d,0,Xp_to_j3,5,Xch_to_j3,dofs1)
  joint6 = Joint(RevoluteJoint,5,Xp_to_j2,6,Xch_to_j2,dofs4)

  @test ismoving(joint3)
  @test !ismoving(joint4)

  joints = [joint2,joint6,joint3,joint4,joint5,joint1]


  body1 = Ellipse(1.0,0.2,200)
  body2 = deepcopy(body1)
  body3 = deepcopy(body1)
  body4 = deepcopy(body1)
  body5 = deepcopy(body1)
  body6 = deepcopy(body1)
  bl = BodyList([body1,body2,body3,body4,body5,body6])

  ls = RigidBodyMotion(joints,bl)

  lsid = 1 # this system should be internally in motion
  @test is_system_in_relative_motion(lsid,ls)

  lsid = 2 # this system should be internally fixed
  @test !is_system_in_relative_motion(lsid,ls)

  # test saving and loading in JLD2
  filen = "testmotion.jld2"
  save(filen,"motion",ls)
  data = load(filen)
  ls_new = data["motion"]

  lsid = 1
  @test is_system_in_relative_motion(lsid,ls_new)
  lsid = 2 # this system should be internally fixed
  @test !is_system_in_relative_motion(lsid,ls_new)
  @test ismoving(ls_new.joints[3])



  ufcn(x,y,t) = 0.25*2π*x*y*cos(2π*t)
  vfcn(x,y,t) = 0.25*2π*(x^2-y^2)*cos(2π*t)
  def = DeformationMotion(ufcn,vfcn)
  X = MotionTransform([0,0],0)
  joint = Joint(X)
  ls = RigidBodyMotion(joint,body1,def)

  @test ismoving(ls)

  defs = AbstractDeformationMotion[NullDeformationMotion(),
                                   NullDeformationMotion(),
                                   NullDeformationMotion(),
                                   NullDeformationMotion(),
                                   def,
                                   NullDeformationMotion()]
  ls = RigidBodyMotion(joints,bl,defs)
  lsid = 2 # now this system is no longer internally fixed, because body 5 deforms
  @test is_system_in_relative_motion(lsid,ls)

  x = init_motion_state(bl,ls)
  Xl = body_transforms(x,ls)
  X1_to_2 = rebase_from_inertial_to_reference(Xl[2],x,ls,1)
  X0_to_2 = X1_to_2*Xl[1]
  @test X0_to_2.x ≈ Xl[2].x && X0_to_2.R ≈ Xl[2].R

  X1_to_4 = rebase_from_inertial_to_reference(Xl[4],x,ls,1)
  X0_to_4 = X1_to_4*Xl[1]
  @test X0_to_4.x ≈ Xl[4].x && X0_to_4.R ≈ Xl[4].R

end
