#=
    This creates the type for a 2-level spin system with N sites.
=#

include("../templates/spin.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateTransitionRates!(template)
    for i = 1:template.size
        if template.state[i] == 1
            template.transitionRates[i] = 1 - template.c
        else
            template.transitionRates[i] = template.c
        end
    end
end

function updateTransitionRates!(template, idx)
    if template.state[idx] == 1
        template.transitionRates[idx] = 1 - template.c
    else
        template.transitionRates[idx] = template.c
    end
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
