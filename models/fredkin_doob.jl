#=
    This creates the type for a 2-level spin system with N sites.
=#

using ITensors
include("../templates/SEP_doob.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateOriginalTransitionRates!(template, idxMin, idxMax)
    for idx = idxMin:idxMax
        if idx == 1 || idx == template.size - 1
            constraint = 0.0
        else
            constraint = convert(Float64, template.state[idx] != template.state[idx+1])
            constraint *= (1 + template.state[idx-1] - template.state[idx+2])
        end

        rate = constraint * ((1-template.c) * template.state[idx] + template.c * template.state[idx+1])
        template.originalTransitionRates[idx] = rate
    end
end


function updateOriginalTransitionRates!(template)
    updateOriginalTransitionRates!(template, 1, template.size-1)
end

function updateOriginalTransitionRates!(template, idx)
    updateOriginalTransitionRates!(template, max(1, idx-2), min(template.size-1, idx+2))
end
