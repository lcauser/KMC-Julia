#=
    This creates the type for a 2-level spin system with N sites.
=#

include("../templates/spin_doob.jl")

"""
    updateOriginalTransitionRates!(template, [idx,])

Updates the original transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateOriginalTransitionRates!(template, idxMin, idxMax)
    for i = idxMin:idxMax
        constraint = Float64(0)
        if i != 1
            constraint += template.state[i-1]
        end
        if i != template.size
            constraint += template.state[i+1]
        end

        if template.state[i] == 1
            template.originalTransitionRates[i] = constraint * (1 - template.c)
        else
            template.originalTransitionRates[i] = constraint * template.c
        end
    end
end

function updateOriginalTransitionRates!(template, idx)
    idxMin = max(1, idx-1)
    idxMax = min(template.size, idx+1)

    updateOriginalTransitionRates!(template, idxMin, idxMax)
end

function updateOriginalTransitionRates!(template)
    updateOriginalTransitionRates!(template, 1, template.size)
end
