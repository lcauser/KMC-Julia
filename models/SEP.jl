#=
    This creates the type for a SEP with N sites.
=#

include("../templates/SEP.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateTransitionRates!(template, idxMin, idxMax)
    for idx = idxMin:idxMax
        constraint = template.state[idx] != template.state[idx+1]
        rate = constraint * ((1-template.c) * state[idx] + template.c * state[idx+1])
        template.transitionRates[idx] = rate
    end
end


function updateTransitionRates!(template)
    updateTransitionRates!(template, 1, template.size-1)
end

function updateTransitionRates!(template, idx)
    updateTransitionRates!(template, max(1, idx-1), min(template.size-1, idx+1))
end


"""
    equilibriumState(template)

Generates a configuration from equilibrium.

### NEEDS CREATING ###
"""
function equilibriumState(template)
    state = zeros(Bool, template.size)
    for i = 1:template.size
        state[i] = rand(Float64) < template.c ? 1 : 0
    end
    return state
end
initialState(template) = equilibriumState(template)
