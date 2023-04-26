using Ballomatic9000, Plots, LinearAlgebra

import Logging
Logging.disable_logging(Logging.Info)

function circle(x, y, r; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end



function run(;L = 100.0, nballs = 100, vth = 1.0)
  domain = Domain(L, L)
  plotkwargs = (aspect_ratio=:equal, legend=false, line=nothing,
                grid=false, xlim=(0, L), ylim=(0,L))
  balls = Vector{Ball}()
  sizehint!(balls, nballs)
  for i in 1:nballs
    while true
      ball = Ball(rand(2) .* L, vth*randn(2), radius=1.0, mass=rand())
      overlaps = false
      for a in balls
        # access member variable = naughty
        overlaps |= norm(a.position - ball.position) < a.radius + ball.radius
        overlaps && break
      end
      if !overlaps
        push!(balls, ball)
        break
      end
    end
  end
  dt = 1 / 100 # 100 Hz
  forces = [zeros(2) for _ in 1:nballs]
  @gif for _ in 1:2^9 # iteration
    simulate!(balls, forces, domain, dt)
    # access member variable = naughty
    plot([circle(a.position[1], a.position[2], a.radius) for a in balls];
         plotkwargs...)
  end every 8
end

run(;vth=10)
