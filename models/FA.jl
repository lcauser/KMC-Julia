#=
    This creates the type for a 2-level spin system with N sites.
=#

include("../templates/spin.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateTransitionRates!(template, idxMin, idxMax)
    for i = idxMin:idxMax
        constraint = Float64(0)
        if i != 1
            constraint += template.state[i-1]
        end
        if i != template.size
            constraint += template.state[i+1]
        end

        if template.state[i] == 1
            template.transitionRates[i] = constraint * (1 - template.c)
        else
            template.transitionRates[i] = constraint * template.c
        end
    end
end

function updateTransitionRates!(template, idx)
    idxMin = max(1, idx-1)
    idxMax = min(template.size, idx+1)

    updateTransitionRates!(template, idxMin, idxMax)
end

function updateTransitionRates!(template)
    updateTransitionRates!(template, 1, template.size)
end


"""
    equilibriumState(template)

Generates a configuration from equilibrium.
"""
function equilibriumState(template)
    state = zeros(Bool, template.size)
    for i = 1:template.size
        state[i] = rand(Float64) < template.c ? 1 : 0
    end
    return state
end
initialState(template) = equilibriumState(template)
