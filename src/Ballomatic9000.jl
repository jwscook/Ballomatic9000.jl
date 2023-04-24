module Ballomatic9000

using LinearAlgebra, PolynomialRoots

const GRAVITY = [0.0, -9.81]

struct Ball{T}
  position::T
  velocity::T
  radius::Float64
  mass::Float64
  function Ball(x::T, v::T, r, m) where T
    m <= 0 && throw(ArgumentError("Ball mass must be > 0, but is $m"))
    r <= 0 && throw(ArgumentError("Ball radius must be > 0, but is $r"))
    all(isfinite, x) || throw(ArgumentError("Ball position must finite, but is $x"))
    all(isfinite, v) || throw(ArgumentError("Ball velocity must finite, but is $v"))
    length(x) == 2 || throw(ArgumentError("Ball position must be 2D; x = $x"))
    length(v) == 2 || throw(ArgumentError("Ball velocity must be 2D; v = $v"))
    return new{T}(x, v, r, m)
  end
end
Ball(x, v; radius, mass) = Ball(x, v, radius, mass)


area(b::Ball) = π * b.radius^2
speed(b::Ball) = norm(b.velocity)
momentum(b::Ball) = b.mass * b.velocity
energy(b::Ball) = b.mass/2 * sum(abs2, b.velocity)

function move!(b::Ball, force, dt)
  @. b.position += (b.velocity + (GRAVITY + force / b.mass) / 2 * dt) * dt
end

function touch(a::Ball, b::Ball, rtol=10eps())
  return isapprox(sum(abs2, a.position - b.position),
                  (a.radius + b.radius)^2, atol=0, rtol=rtol)
end

# See https://en.wikipedia.org/wiki/Elastic_collision
function collide!(a::Ball, b::Ball)
  # this function has zero allocations
  x1, x2 = a.position, b.position
  v1, v2 = a.velocity, b.velocity
  m1, m2 = a.mass, b.mass
  Δx² = ΔvΔx = zero(promote_type(eltype(v1), eltype(v2)))
  for i in eachindex(v1)
    Δxi = x1[i] - x2[i]
    ΔvΔx += (v1[i] - v2[i]) * Δxi
    Δx² += Δxi^2
  end
  # TODO: Is domain error the right kind of error to throw?
  iszero(Δx²) && throw(DomainError("The distance between two balls must be > 0"))
  commonfactor = 2 / (m1 + m2) * ΔvΔx / Δx²
  @. a.velocity -= m2 * (x1 - x2) * commonfactor
  @. b.velocity -= m1 * (x2 - x1) * commonfactor
  return nothing
end

function solvequadratic(a, b, c)
  # ax^2 + bx + c = 0
  solutions = if iszero(a) 
    (iszero(b) ? -Inf : -c / b) .* (1, 1)
  else
    sqrtarg = b^2 - 4a*c
    if sqrtarg < 0
      (-Inf, -Inf) # this is safe for the use case, but otherwise a silly thing to do
    else
      ((-b + sqrt(sqrtarg)) / 2a, (-b - sqrt(sqrtarg)) / 2a)
    end
  end
  return solutions
end
#  >>> import sympy
#  >>> x1, x2, v1, v2, a1, a2, t, r1, r2 = sympy.symbols('x1 x2 v1 v2 a1 a2 t r1 r2')
#  >>> expr = ((x1 - x2) + (v1 - v2)*t + (a1 - a2)*t**2 /2)**2 - (r1 + r2)**2
#  >>> sympy.Poly(expr, t).all_coeffs()
#  [a1**2/4 - a1*a2/2 + a2**2/4, a1*v1 - a1*v2 - a2*v1 + a2*v2, a1*x1 - a1*x2 - a2*x1 + a2*x2 + v1**2 - 2*v1*v2 + v2**2, 2*v1*x1 - 2*v1*x2 - 2*v2*x1 + 2*v2*x2, -r1**2 - 2*r1*r2 - r2**2 + x1**2 - 2*x1*x2 + x2**2]
function timeofcollision(a::Ball, b::Ball, forcea, forceb)
  x1, x2 = a.position, b.position
  v1, v2 = a.velocity, b.velocity
  a1, a2 = GRAVITY .+ forcea / a.mass, GRAVITY .+ forceb / b.mass
  r1, r2 = a.radius, b.radius
  coeffs = [dot(a1, a1)/4 - dot(a1, a2)/2 + dot(a2, a2)/4,
            dot(a1, v1) - dot(a1, v2) - dot(a2, v1) + dot(a2, v2),
            dot(a1, x1) - dot(a1, x2) - dot(a2, x1) + dot(a2, x2) + dot(v1, v1) - 2dot(v1, v2) + dot(v2, v2),
            2dot(v1, x1) - 2dot(v1, x2) - 2dot(v2, x1) + 2dot(v2, x2),
            -dot(r1, r1) - 2dot(r1, r2) - dot(r2, r2) + dot(x1, x1) - 2dot(x1, x2) + dot(x2, x2)]
  @show times_of_potential_collisions = sort(real.(roots(coeffs)))
  index_of_earliest_collision_time = findfirst(x->x > 0, times_of_potential_collisions)
  if isnothing(index_of_earliest_collision_time)
    return Inf
  else
    return times_of_potential_collisions[index_of_earliest_collision_time]
  end
end

struct Domain
  xmax::Float64
  ymax::Float64
end

function timeofcollision(a::Ball, force, d::Domain)
  r1 = a.radius
  x1 = a.position
  v1 = a.velocity
  a1 = GRAVITY .+ force / a.mass
  # factor of 1/2 neceause x0 - x1 + ut + 1/2t^2 = 0
  tx0 = solvequadratic(a1[1] / 2, v1[1], x1[1] - r1) # collision at x=0
  ty0 = solvequadratic(a1[2] / 2, v1[2], x1[2] - r1) # collision at y=0
  tx1 = solvequadratic(a1[1] / 2, v1[1], x1[1] - d.xmax + r1) # collision at x=max
  ty1 = solvequadratic(a1[2] / 2, v1[2], x1[2] - d.ymax + r1) # collision at y=max
  return min(Inf, max(tx0[1], ty0[1], tx1[1], ty1[1],
                      tx0[2], ty0[2], tx1[2], ty1[2]))
end

export Ball, collide!, move!, area, speed, momentum, energy, Domain

end # module Ballomatic9000
