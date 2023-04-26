using Ballomatic9000, Plots

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

  balls = [Ball(rand(2) .* L, vth*randn(2), radius=1.0, mass=rand()) for _ in 1:nballs]
  dt = 1 / 100 # 100 Hz
  forces = [zeros(2) for _ in 1:nballs]
  @gif for _ in 1:2^10 # iteration
    simulate!(balls, forces, domain, dt)
    plot([circle(a.position[1], a.position[2], a.radius) for a in balls];
         plotkwargs...)
  end every 8
end

run(;vth=4)
