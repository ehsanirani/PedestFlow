using DrWatson
@quickactivate "PedestFlow"

using Agents, LinearAlgebra, StaticArrays, Random
using LazySets

@enum Destination b_exit d_exit
@agent Person ContinuousAgent{2} begin
    σ::Float64
    m::Float64
    speed::Float64
    f::NTuple{2,Float64}
    dest::Destination
    Rv::Float64
    θ::Float64
end
Person(; σ=0.5, m=1.0, speed=1.0, f=(0.0, 0.0),
         dest=:b_exit, Rv=1.0, θ=π) = Person(σ, m, speed, f, dest, Rv, θ)

function initialize_model(;
    n_persons=1,
    extent=(30, 30),
    offset=2.5,
    spacing=5.0,
    ϕ=0.1,
    dt=0.1,
    β=1.0,
    γ = 1.0,
    σ=0.5,
    m=1.0,
    θ=π,
    Rv=2.0,
    seed=12345
)
#    dest_lines = Dict(
#        :b_exit => LazySets.Line(from=[12.5, 27.5], to=[17.5, 27.5]),
#        :d_exit => LazySets.Line(from=[27.5, 12.5], to=[27.5, 17.5])
#    )
    walls = Dict(
        :w1 => [(-5.0, 10.0), (10.0, 10.0)],  # l2
        :w2 => [(-5.0, 15.0), (10.0, 15.0)],  # l1
        :w3 => [(15.0, 10.0), (25.0, 10.0)], # r2
        :w4 => [(15.0, 15.0), (25.0, 15.0)], # r1
        :w5 => [(10.0, -5.0), (10.0, 10.0)],  # b1
        :w6 => [(15.0, -5.0), (15.0, 10.0)],  # b2
        :w7 => [(10.0, 15.0), (10.0, 25.0)], # t1
        :w8 => [(15.0, 15.0), (15.0, 25.0)]  # t2
     )
     for (k, v) in walls
        walls[k][1] = walls[k][1] .+ offset
        walls[k][2] = walls[k][2] .+ offset
     end
    space2d = ContinuousSpace(extent, spacing)
    model = ABM(Person, space2d, #scheduler = Schedulers.Randomly(), 
        properties=Dict(
            :ϕ => ϕ,
            :dt => dt,
            :γ => γ,
            :β => β,
            :σ => σ,
            :m => m,
            :Rv => Rv,
            :θ => θ,
            :buffer_size => offset,
            #:dest_lines => dest_lines,
            :walls => walls),
            rng = MersenneTwister(seed))
    for i in 1:n_persons
        r = (rand(model.rng) * 4.8 + 10.1 + offset, 0.0)
        v = (0.0, 1.0)
        speed = rand(model.rng) * 0.2 + 1.1#1.0,
        add_agent!(r, model, v, σ, m, speed, (0.0, 0.0), b_exit, Rv, θ)
    end
    return model
end

"""
Distance from a line segment:
    function dist2wall(p, a, b, model)

p: point
a, b: two end of the line segment
model: ABM model

returns distance and the vector connecting p to a-b: dist, NTuple
"""
function dist2wall(p, a, b)
    # if p on a-b

    ab = b .- a
    ab_mag = norm(ab)

    ap = p .- a
    # x has the shortest distance to p
    x = a .+ (dot(ab, ap) .* ab) ./ ab_mag^2 

    # if x ∈ ab:
    line = LazySets.LineSegment([a...], [b...])
    if [x...] ∈ line
        xp = p .- x
        return norm(xp), xp
    end
    # if x ∉ ab:
    bp = p .- b
    ap_mag = norm(ap)
    bp_mag = norm(bp)
    if ap_mag < bp_mag
        return ap_mag, ap
    else
        return bp_mag, bp
    end
end

"""
Interaction with walls: A harmonic repulsion.
"""
function interact_wall!(p::Person, model::ABM)
    rcut = p.σ * 1
    k = 100.0
    #    t1 t2 
    # l1       r1
    # l2       r2
    #    b1 b2
    for (ks, w) in model.walls
        dik, rik = dist2wall(p.pos, w[1], w[2])
        if dik > rcut
            continue
        end
        p.f = p.f .+ k .* (p.σ .- dik) .*(rik./dik)
    end

end

"""
Calculate the predicted collision point and returns 
the corresponding unit vector pointing to it.
"""
function get_pij(pi, pj, rij, model)
    dij = 1.5*pi.σ
    τ_max = 100.0
    vij = pi.vel .- pj.vel
    vij_mag = norm(vij)
    a = vij_mag^2
    b = 2 * rij * vij_mag
    c = rij^2 - dij^2

    b2_4ac = b^2 - 4 * a * c
    if b2_4ac < 0.0
        return (0.0, 0.0)
    end
    tcoll1 = (-b + sqrt(b2_4ac))/(2a)
    tcoll2 = (-b - sqrt(b2_4ac))/(2a)
    tcoll = min(tcoll1, tcoll2)
    if tcoll<0
        tcoll = max(tcoll1, tcoll2)
    end
    if tcoll < 0 || tcoll > τ_max
        return (0.0, 0.0)
    end
    pij = (pi.pos .+ pi.vel .* tcoll) .- (pj.pos .+ pj.vel .* tcoll)
    return pij ./ norm(pij)
end

"""
Advancing the agent!
Integrating equations of motion takes place here
"""
function agent_step!(person, model)
    # interactios + driving
    f_inter = person.f
    # additional damping, if desired
    f_damp = -model.γ .* person.vel
    # random force
    ζ = (randn(model.rng, 2)..., )
    f_rand = model.β .* ζ
    # total force
    f_tot = (f_inter .+ f_damp .+ f_rand)
    a = f_tot ./ person.m
    
    # limit the acceleration to 2 m/s^2
    a_mag = norm(a)
    if a_mag > 2.0
        a = (a ./ a_mag) .* 2.0
    end
    
    # integrate
    Δv = a .* model.dt
    person.vel = person.vel .+ Δv
    Δr = person.vel .* model.dt
    walk!(person, Δr, model)
    return person.pos
end

"""
Calculate and add the relevant force to drive
the person toward its destination.
"""
function agent_target!(person, model)
    α = 1.0

    # Driving force to the destination
    if person.dest == b_exit
        da = (11.0, 26.0) .+ model.buffer_size
        db = (14.0, 26.0) .+ model.buffer_size
    else
        da = (26.0, 11.0) .+ model.buffer_size
        db = (26.0, 14.0) .+ model.buffer_size
    end
    d_idest, r_idest = dist2wall(person.pos, da, db)

    # Desired velocity, tuned by the prescribed speed
    v_des = person.speed .* r_idest ./ -d_idest
    ai0 = α .* (v_des .- person.vel)
    person.f = person.f .+ ai0
end

"""
Interaction with other people is calculated here
"""
function agent_interact!(person, model)
    k = 1000.0
    # Interactions with other people
    # Obtain the ids of neighbors
    neighbor_ids = nearby_ids(person, model, person.Rv)  
    for id in neighbor_ids
        pj = model[id]
        rij =  person.pos .- pj.pos
        dij = norm(rij)

        aij = (0.0, 0.0)

        cos_γij = sum(person.vel .* (-1 .* rij)) / (norm(person.vel) * norm(rij))
        # bound cos(γij) to [-1, 1]
        if cos_γij > 1
            cos_γij = 1.0
        elseif cos_γij<-1
            cos_γij = -1.0
        end
        γij = acos(cos_γij)
        # interaction matters only if neighbor pj is in the view field of the person
        if (person.θ / 2) < (γij)
            continue
        end
        # get the unit vector in the direction of the predicted collision
        pij = get_pij(person, pj, dij, model)

        aij = aij .+ (k/(dij^2)) .* pij 
       
        person.f = person.f .+ aij
    end

end

"""
Time evolution of the model.
It calls force calculations for pedestrians,
and add new people with rate of ϕ,
and remove poeple who arrived their destination.
"""
function model_step!(model)
    extent = model.space.extent

    # remove people who arrived
    remove_people!(model)

    # reset forces and interact with walls
    for p in allagents(model)
        p.f = (0.0, 0.0)
        agent_target!(p, model)
        agent_interact!(p, model)
        interact_wall!(p, model)
    end

    # add new people with flux ϕ
    add_people!(model)
end

"""
Returns true if person is arrived to the destination.
"""
function is_arrived(p::Person, model::ABM)
    if p.dest == b_exit
        if p.pos[2] > 25.0 + model.buffer_size
            return true
        end
    else
        if p.pos[1] > 25.0 + model.buffer_size
            return true
        end
    end
    return false
end

"""
Remove people from the simulation if arrived 
to their destination.
"""
function remove_people!(model)
    if model.agents.count == 0
        return
    end
    for p in allagents(model)
        if is_arrived(p, model)
            kill_agent!(p, model)
        end
    end
end

"""
Add new people to the simulation with
the rate of ϕ
"""
function add_people!(model)
    ϕ = model.ϕ
    prob = model.dt * ϕ
    if rand(model.rng) < 0.5*prob
        r = (rand(model.rng) * 4. + 13.0, 0.0)
        v = (0.0, 1.0)
        speed = rand(model.rng) * 0.2 + 1.1#1.0,
        add_agent!(r, model, v, model.σ, model.m, speed, (0.0, 0.0), b_exit, model.Rv, model.θ)
    end
    if rand(model.rng) < 0.5*prob
        r = (0.0, rand(model.rng) * 4. + 13.0)
        v = (1.0, 0.0)
        speed = rand(model.rng) * 0.2 + 1.1#1.0,
        add_agent!(r, model, v, model.σ, model.m, speed, (0.0, 0.0), d_exit, model.Rv, model.θ)
    end
end

## Visulaization helpers
"""
unit_vec(x::NTuple{2, Real})

returns the unit vector in the direction of x
"""
function unit_vec(x)
    return x ./ norm(x)
end

using InteractiveDynamics, GLMakie
const bird_polygon = Polygon(Point2f[(-0.25, -0.25), (0.25, 0), (-0.25, 0.25)])
function bird_marker(b::Person)
    φ = atan(b.vel[2], b.vel[1]) #+ π/2 + π
    scale(rotate2D(bird_polygon, φ), 1)
end

function ac(p::Person)
if p.dest == b_exit
    return :blue
end
return :green
end

function disk_marker(p::Person)
    t = LinRange(0, 2π, 50)
    x = cos.(t)
    y = sin.(t)
    disk = permutedims([x y])
    coords = [Point2f(x, y) for (x, y) in zip(disk[1, :], disk[2, :])]

    φ = atan(p.vel[2], p.vel[1])
    scale(rotate2D(Polygon(coords), φ), 0.25)
end

function static_preplot!(ax, model)
    offset = model.buffer_size
    lines!(ax, [10, 10, NaN, 10, 10] .+ offset, [0, 10, NaN, 15, 25] .+ offset, color=:black)
    lines!(ax, [15, 15, NaN, 15, 15] .+ offset, [0, 10, NaN, 15, 25] .+ offset, color=:black)
    lines!(ax, [0, 10, NaN, 15, 25] .+ offset, [10, 10, NaN, 10, 10] .+ offset, color=:black)
    lines!(ax, [0, 10, NaN, 15, 25] .+ offset, [15, 15, NaN, 15, 15] .+ offset, color=:black)
    hidedecorations!(ax) # hide tick labels etc.
end

  
  
