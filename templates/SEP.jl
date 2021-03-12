#=
    The template for 2-level systems
=#

mutable struct SEP{}
    size::Int
    particles::Int
    c::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    escapeRate::Float64
end

"""
    Spin(N, c)

Return a two-level spin structure with N sites and c bias.
"""
function SEP(N, m, c)
    return SEP(N, m, c, zeros(N), zeros(N-1), 0.0)
end


"""
    setState!(template, state)

Change the state in the model
"""
function setState!(template, state)
    template.state = copy(state)
end


"""
    transition!(template, idx)

Transition the state in the model according to the unique identifier idx.
"""
function transition!(template, idx)
    template.state[idx] = 1 - template.state[idx]
    template.state[idx+1] = 1 - template.state[idx+1]
end
