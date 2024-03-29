#=
    The template for 2-level systems
=#

mutable struct Spin{}
    size::Int
    c::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    escapeRate::Float64
end

"""
    Spin(N, c)

Return a two-level spin structure with N sites and c bias.
"""
function Spin(N, c)
    return Spin(N, c, zeros(N), zeros(N), 0.0)
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
end
