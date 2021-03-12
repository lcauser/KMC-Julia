#=
    This creates the type for a 2-level spin system with N sites.
=#

using ITensors
include("../templates/SEP.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateTransitionRates!(template, idxMin, idxMax)
    for idx = idxMin:idxMax
        if idx == 1 || idx == template.size - 1
            constraint = 0.0
        else
            constraint = convert(Float64, template.state[idx] != template.state[idx+1])
            constraint *= (1 + template.state[idx-1] - template.state[idx+2])
        end

        rate = constraint * ((1-template.c) * template.state[idx] + template.c * template.state[idx+1])
        template.transitionRates[idx] = rate
    end
end


function updateTransitionRates!(template)
    updateTransitionRates!(template, 1, template.size-1)
end

function updateTransitionRates!(template, idx)
    updateTransitionRates!(template, max(1, idx-2), min(template.size-1, idx+2))
end


"""
    equilibriumState(template)

Generates a configuration from equilibrium for c = 0.5

"""
function equilibriumState(template)
    state = zeros(Bool, template.size)
    height = 0
    for i = 1:template.size
        denom = (2*(height+1)*(template.size - (i-1)))
        p1 = (height + 2)*(template.size-(i-1)-height) / denom
        choice = rand() < p1
        height += 2*convert(Int64, choice)-1
        state[i] = choice
    end

    return state
end
initialState(template) = equilibriumState(template)


"""
    sampleState(MPS)

Samples a configuration from an MPS
"""
function sampleState(template, MPS)
    config = sample(MPS)
    config = [site == 1 for site = config]
    return config
end
