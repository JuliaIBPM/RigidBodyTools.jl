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

  r = rand(2)
  Θ = 0.0
  TM = motion_transform_matrix_2d(r,Θ)
  vA = rand(3)

  vB = TM*vA

  # Should give velocity based at origin of B
  @test vB ≈ vA .+ [0.0,-vA[1]*r[2],vA[1]*r[1]]

  TF = force_transform_matrix_2d(r,Θ)
  fA = rand(3)
  fB = TF*fA

  # Should give moment about origin of B
  @test fB ≈ fA .- [r[1]*fA[3]-r[2]*fA[2],0.0,0.0]




end
