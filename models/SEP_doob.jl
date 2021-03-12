#=
    This creates the type for a SEP with N sites.
=#

include("../templates/SEP_doob.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateOriginalTransitionRates!(template, idxMin, idxMax)
    for idx = idxMin:idxMax
        constraint = convert(Float64, template.state[idx] != template.state[idx+1])
        rate = constraint * ((1-template.c) * template.state[idx] + template.c * template.state[idx+1])
        template.originalTransitionRates[idx] = rate
    end
end


function updateOriginalTransitionRates!(template)
    updateOriginalTransitionRates!(template, 1, template.size-1)
end

function updateOriginalTransitionRates!(template, idx)
    #updateOriginalTransitionRates!(template, max(1, idx-1), min(template.size-1, idx+1))#
    updateOriginalTransitionRates!(template, 1, template.size-1)
end
