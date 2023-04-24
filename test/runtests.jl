using LinearAlgebra, Test
using Ballomatic9000
const GRAVITY = Ballomatic9000.GRAVITY

@testset "Ballomatic9000" begin

@testset "Check for brainfarts" begin
  for _ in 1:10
    x1, x2 = randn(2), randn(2)
    v1, v2 = randn(2), randn(2)
    @test dot(x1 - x2, v1 - v2) == dot(x2 - x1, v2 - v1)
  end
end

@testset "Collisions" begin
  @testset "Conservation" begin
    for _ in 1:10
      x1, x2 = rand(2), rand(2)
      @assert !iszero(norm(x1 - x2)) "This doesnt tes degenerate positions"
      a = Ball(x1, rand(2), radius=1.0, mass=rand())
      b = Ball(x2, rand(2), radius=1.0, mass=rand())
      momentum_before = momentum(a) + momentum(b)
      energy_before = energy(a) + energy(b)
      collide!(a, b)
      momentum_after = momentum(a) + momentum(b)
      energy_after = energy(a) + energy(b)
      @test momentum_before ≈ momentum_after
      @test energy_before ≈ energy_after
    end
  end

  @testset "Conservation" begin
    x = rand(2)
    a = Ball(x, rand(2), rand(), rand())
    b = Ball(x, rand(2), rand(), rand())
    @test_throws DomainError collide!(a, b)
  end

  @testset "Collide with domain wall" begin
    d = Domain(10, 20)
    a = Ball([5.0, 5.0], [0.0, 0.0], 1.0, 1.0)
    # fall at 9.81
    #s = x + ut + 1/2 a t^2
    #sqrt(8/a) = t
    g = -GRAVITY[2]
    @test sqrt(8 / g) ≈ Ballomatic9000.timeofcollision(a, 0.0, d)
    #20-1 = 5 + 1/2 a t^2
    @test sqrt(28 / g) ≈ Ballomatic9000.timeofcollision(a, [0.0, 2g], d)
  end

  @testset "Collide with eachother" begin
    a = Ball([0.0, 0.0], [1.0, 0.0], 1.0, 1.0)
    b = Ball([4.0, 0.0], [-1.0, 0.0], 1.0, 1.0)
    # a starts a 0, b starts at 4, they move together at speed 1 each and collide
    # when radii touch after travelling 1 space unit each in 1 time unit
    #s = x + ut + 1/2 a t^2
    @test 1 ≈ Ballomatic9000.timeofcollision(a, b, 0, 0)
    a = Ball([0.0, 0.0], [1.0, 0.0], 1.0, 1.0)
    b = Ball([8.0, 0.0], [-1.0, 0.0], 1.0, 1.0)
    @test 3 ≈ Ballomatic9000.timeofcollision(a, b, 0, 0)
    a = Ball([0.0, 0.0], [2.0, 0.0], 1.0, 1.0)
    b = Ball([8.0, 0.0], [-2.0, 0.0], 1.0, 1.0)
    @test 1.5 ≈ Ballomatic9000.timeofcollision(a, b, 0, 0)
  end


end

end
