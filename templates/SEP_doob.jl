#=
    The template for interacting hard-core bosons
=#

using ITensors

mutable struct SEPDoob{}
    size::Int
    particles::Int
    c::Float64
    s::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    originalTransitionRates::Array{Float64}
    escapeRate::Float64
    originalEscapeRate::Float64
    psi::MPS
    leftComponent::Float64
    leftBlocks::Array{ITensor}
    rightBlocks::Array{ITensor}
end

"""
    SpinDoob(N, c, s, psi)

Return a two-level spin structure with N sites and c, s bias. Requires an MPS
probability vector psi.
"""
function SEPDoob(N, m, c, s, psi)
    return SEPDoob(N, m, c, s, zeros(N-1), zeros(N-1), zeros(N-1), 0.0, 0.0,
                 psi, 0, Array{ITensor}(undef, N), Array{ITensor}(undef, N))
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



function updateLeftBlock!(template, idx)
    # Get the state and the psi block
    state = template.state[idx]
    A = template.psi[idx]
    s = siteind(template.psi, idx)
    projA = ITensor(s)
    projA[s[2 - state]] = 1

    # Fetch the previous block and multiply by A and projA
    block = idx == 1 ? 1 : copy(template.leftBlocks[idx-1])
    block *= A
    block *= dag(projA)
    template.leftBlocks[idx] = block
end

function updateRightBlock!(template, idx)
    # Get the state and the psi block
    state = template.state[idx]
    A = template.psi[idx]
    s = siteind(template.psi, idx)
    projA = ITensor(s)
    projA[s[2 - state]] = 1

    # Fetch the previous block and multiply by A and projA
    block = idx == template.size ? 1 : copy(template.rightBlocks[idx+1])
    block *= A
    block *= dag(projA)
    template.rightBlocks[idx] = block
end

function psiToLeft(template, state)
    # Needs doing for SEP
    return 1.0
end

function leftComponent(template, idx, current=false)
    # Get the psi component at the sites and construct projectors
    state1 = current ? template.state[idx] : 1 - template.state[idx]
    state2 = current ? template.state[idx+1] : 1 - template.state[idx+1]
    A1 = template.psi[idx]
    A2 = template.psi[idx+1]
    s1 = siteind(template.psi, idx)
    s2 = siteind(template.psi, idx+1)
    projA1 = ITensor(s1)
    projA1[s1[2 - state1]] = 1
    projA2 = ITensor(s2)
    projA2[s2[2 - state2]] = 1

    # Get the left and right blocks
    leftBlock = idx == 1 ? 1 : template.leftBlocks[idx-1]
    rightBlock = idx == template.size - 1 ? 1 : template.rightBlocks[idx+2]

    # Multiply together
    prod = leftBlock * A1
    prod *= dag(projA1)
    prod *= A2
    prod *= dag(projA2)
    prod *= rightBlock
    prod = real(scalar(prod))

    # Get conversion factor
    state = copy(template.state)
    state[idx] = current ? state[idx] : 1-state[idx]
    conversion = psiToLeft(template, state)
    return abs(prod) + 10^-18
end


function updateTransitionRates!(template, blockLeft, blockRight)
    # Find the first and last sites which can flip
    idxLeft = findfirst([x != 0 for x = template.originalTransitionRates])
    idxRight = findlast([x != 0 for x = template.originalTransitionRates])

    # Build up the left and right blocks
    for i = blockLeft:idxRight
        updateLeftBlock!(template, i)
    end
    for i = blockRight:-1:idxLeft
        updateRightBlock!(template, i)
    end

    # Update the left component
    template.leftComponent = leftComponent(template, idxLeft, true)

    # Loop through each site and calculate the transition rates
    transitionRates = copy(template.originalTransitionRates)
    for i = 1:template.size-1
        if transitionRates[i] != 0
            left = leftComponent(template, i)
            transitionRates[i] *= exp(-template.s) * left / template.leftComponent
        end
    end

    template.transitionRates = transitionRates
end


"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateTransitionRates!(template, idx)
    updateOriginalTransitionRates!(template, idx)
    updateTransitionRates!(template, idx, min(idx+1, template.size))
end

function updateTransitionRates!(template)
    updateOriginalTransitionRates!(template)
    updateTransitionRates!(template, 1, template.size)
end


"""
    equilibriumState(template)

Generates a configuration from equilibrium.
"""
function equilibriumState(template)
    state = 2 .- sample(template.psi)
    return convert(Array{Bool}, state)
end
initialState(template) = equilibriumState(template)
