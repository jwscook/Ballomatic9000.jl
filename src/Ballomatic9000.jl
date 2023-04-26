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

momentum(b::Ball) = b.mass * b.velocity
energy(b::Ball) = b.mass/2 * sum(abs2, b.velocity)

function move!(b::Ball, force, dt)
  @. b.position += (b.velocity + (GRAVITY + force / b.mass) / 2 * dt) * dt
end


"""
    calculate_time_from_coeffs(coeffs)

Coeffs of a polynomial of form ∑ᵢ₌₀ᴺ cᵢtⁱ

Return the smallest non-zero positive root corresonding the time of collision.

...
# Arguments
- `coeffs`: Vector-like container of numbers
...

# Example
```julia
```
"""
function calculate_time_from_coeffs(coeffs)
  all(iszero, coeffs) && return Inf
  rts = roots(coeffs)
  # if a solution has an imaginary component then subtract in from the real part
  rts .= remove_small_imaginary_part.(rts)
  rts .-= (abs.(imag.(rts)) .> 0) * Inf
  # now only take the real solutions (bad solutions collided at infinity in the past)
  times_of_potential_collisions = sort(real.(rts))
  # any tiny times until collisions, e.g. 1e-13, will be just after the last collision
  times_of_potential_collisions .= remove_small_timestep.(times_of_potential_collisions)
  # retrieve the soonest collision in the future
  index_of_earliest_collision_time = findfirst(x->x > 0, times_of_potential_collisions)
  if isnothing(index_of_earliest_collision_time)
    return Inf
  else
    return times_of_potential_collisions[index_of_earliest_collision_time]
  end
end


"""
    collide!(a::Ball,b::Ball)

Collide two balls off one another in an elastic collision that
alters their velocities. No temporal advance is taken.

# See https://en.wikipedia.org/wiki/Elastic_collision

...
# Arguments
- `a::Ball`: 
- `b::Ball`: 
...

# Example
```julia
```
"""
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
  iszero(Δx²) && throw(DomainError("The distance between two balls cannot be = 0"))
  commonfactor = 2 / (m1 + m2) * ΔvΔx / Δx²
  @. a.velocity -= m2 * (x1 - x2) * commonfactor
  @. b.velocity -= m1 * (x2 - x1) * commonfactor
  return nothing
end


function remove_small_imaginary_part(z)
  iszero(z) && return z
  nz = norm(z)
  normz = z ./ nz
  normz = trunc.(real.(normz), digits=16) .+ im * trunc.(imag.(normz), digits=16)
  output = normz .* nz
  @assert output ≈ z "$output is not ≈ $z"
  return output
end

remove_small_timestep(dt) = dt < 1e-8 ? zero(dt) : dt

"""
    timeofcollision(a::Ball,b::Ball,forcea,forceb)

Calculate all the collisions that the two accelerating balls may undertake and
return the soonest collision time.

Position at time t for ball a, x1' = x1 + v1t + a1/2 * t^2
Position at time t for ball b, x2' = x2 + v2t + a2/2 * t^2

Difference in positions at time t, (x1' - x2') = (x1 - x2) + (v1 - v2)t + (a1 - a2)/2 * t^2
Square of distance between balls at time t, (x1' - x2')^2 = ((x1 - x2) + (v1 - v2)t + (a1 - a2)/2 * t^2)^2
At the time of collision that squared distance is (r1 + r2)^2 = (x1' - x2')^2

  >>> import sympy
  >>> x1, x2, v1, v2, a1, a2, t, r1, r2 = sympy.symbols('x1 x2 v1 v2 a1 a2 t r1 r2')
  >>> expr = ((x1 - x2) + (v1 - v2)*t + (a1 - a2)*t**2 / 2)**2 - (r1 + r2)**2
  >>> sympy.Poly(expr, t).all_coeffs()
  [a1**2/4 - a1*a2/2 + a2**2/4, a1*v1 - a1*v2 - a2*v1 + a2*v2, a1*x1 - a1*x2 - a2*x1 + a2*x2 + v1**2 - 2*v1*v2 + v2**2, 2*v1*x1 - 2*v1*x2 - 2*v2*x1 + 2*v2*x2, -r1**2 - 2*r1*r2 - r2**2 + x1**2 - 2*x1*x2 + x2**2]
...
# Arguments
- `a::Ball`: ball to collide
- `b::Ball`: ball to collide
- `forcea`: force on ball a
- `forceb`: force on ball b
...

# Example
```julia
```
"""
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
  reverse!(coeffs) # reverse, julia lib expects opposite order to python lib :(
  return calculate_time_from_coeffs(coeffs)
end

"""
The Domain representing the locations of the walls, which are at:
 - x = 0
 - y = 0
 - x = xmax
 - y = ymax
"""
struct Domain
  xmax::Float64
  ymax::Float64
end

"""
    timeofcollision(a::Ball,force,d::Domain)

Calculate the time that the ball will collide with the walls of the domain

...
# Arguments
- `a::Ball`: the ball
- `force`: the force acting on the ball
- `d::Domain`: the domain in which the ball is bouncing
...

# Example
```julia
```
"""
function timeofcollision(a::Ball, force, d::Domain)
  r = a.radius
  x, y = a.position
  vx, vy = a.velocity
  ax, ay = GRAVITY .+ force / a.mass
  # factor of 1/2 because x0 - x1 + ut + a/2t^2 = 0
  tx0 = calculate_time_from_coeffs([x - r, vx, ax /2])
  ty0 = calculate_time_from_coeffs([y - r, vy, ay /2])
  tx1 = calculate_time_from_coeffs([x - d.xmax + r, vx, ax /2])
  ty1 = calculate_time_from_coeffs([y - d.ymax + r, vy, ay /2])
  return min(tx0, ty0, tx1, ty1)
end

"""
    collide!(a::Ball,d::Domain,atol=sqrt(eps()),rtol=sqrt(eps()))

Collide a ball with the wall of the Domain

...
# Arguments
- `a::Ball`: ball that could collide with the domain wall
- `d::Domain`: domain defining positions of the wall
- `atol=sqrt(eps())`: absolute tolerance to check position of ball wrt walls
- `rtol=sqrt(eps())`: as above but relative tolerance
...

# Example
```julia
```
"""
function collide!(a::Ball, d::Domain, atol=sqrt(eps()), rtol=sqrt(eps()))
  # so ugly!
  if isapprox(a.position[1], 0.0 + a.radius, atol=atol, rtol=rtol)
    a.velocity[1] -= 2a.velocity[1]
  elseif isapprox(a.position[1], d.xmax - a.radius, atol=atol, rtol=rtol)
    a.velocity[1] -= 2a.velocity[1]
  elseif isapprox(a.position[2], 0.0 + a.radius, atol=atol, rtol=rtol)
    a.velocity[2] -= 2a.velocity[2]
  elseif isapprox(a.position[2], d.ymax - a.radius, atol=atol, rtol=rtol)
    a.velocity[2] -= 2a.velocity[2]
  end
end

function simulate!(balls, forces, domain::Domain, dt)
  length(forces) == length(balls) || throw(ArgumentError("Length of forces must be length of balls container"))
  substep = dt
  time = 0.0
  collisionpair = [0, 0]
  # some balls may collide within a time interval dt.
  # so find earliest collision and integrate all balls to that time
  # then keep going until reach dt in the future
  while time < dt
    for (i, a) in enumerate(balls) # while some balls havent yet been integrated through a full dt
      for j in i+1:length(balls) # all other balls, not good n^2, put in a tree
        i == j && continue # don't collide last ball with itself (or any ball with itself for that matter)
        b = balls[j]
        collisiontime = timeofcollision(a, b, forces[i], forces[j]) # get collision time if at all
        @assert collisiontime > 0
        if collisiontime < substep # is this collision soonest?
          collisionpair .= [i, j] # ball i collides with ball j
          substep = min(dt, collisiontime)
        end
      end
      # what if it collided with a wall?
      wallcollisiontime = timeofcollision(a, forces[i], domain)
      if wallcollisiontime < substep
        substep = min(dt, wallcollisiontime)
        collisionpair .= [i, 0] # special value of 0 denotes a wall collision
      end
    end
    substep = min(substep, dt - time) # don't overshoot dt
    iszero(substep) && break
    @assert substep > 0 "substep is $substep"
    for i in eachindex(balls)
      move!(balls[i], forces[i], substep) # move all balls to substep
      if i == collisionpair[1] && iszero(collisionpair[2]) # then do collision with wall
        @info "Ball at position $(balls[i].position) collided with the wall"
        collide!(balls[i], domain) # potentially a collision with a wall
      elseif i == collisionpair[1] # it was a ball-ball collision
        j = collisionpair[2]
        @info "Balls at positions $(balls[i].position) and $(balls[j].position) collided"
        collide!(balls[i], balls[j])
      end
    end
    time += substep # add substep to the time
    fill!(collisionpair, 0) # zero out collision index pair
  end
end

export Ball, collide!, move!, momentum, energy, Domain, simulate!

end # module Ballomatic9000
