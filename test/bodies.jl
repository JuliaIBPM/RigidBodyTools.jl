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




end
