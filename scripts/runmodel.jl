using DrWatson
@quickactivate "PedestFlow"

include(srcdir("ehsan.jl"))

function init_model(;n_persons=1, ϕ=0.2, dt=0.005, σ=0.5, m=1, θ=3π/4, Rv=4.0, β=0.0, seed=1)
    model = initialize_model(n_persons=n_persons, ϕ=ϕ, dt=dt, σ=σ, m=m, θ=θ, Rv=Rv, β=β, seed=seed)
    return model
end




