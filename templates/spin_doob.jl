#=
    The template for 2-level systems
=#

using ITensors

mutable struct SpinDoob{}
    size::Int
    c::Float64
    s::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    originalTransitionRates::Array{Float64}
    escapeRate::Float64
    originalEscapeRate::Float64
    psis
    weights
    leftComponent::Float64
    leftBlocks::Array{Array{ITensor}}
    rightBlocks::Array{Array{ITensor}}
end

"""
    SpinDoob(N, c, s, psi)

Return a two-level spin structure with N sites and c, s bias. Requires an MPS
probability vector psi.
"""
function SpinDoob(N, c, s, psi)
    arr = Array{Array{ITensor}}
    push!(arr, Array{ITensor}(undef, N))
    return SpinDoob(N, c, s, zeros(N), zeros(N), zeros(N), 0.0, 0.0, [psi], [1], 0,
                deepcopy(arr), deepcopy(arr))
end

function SpinDoob(N, c, s, psis)
    arr = Array{Array{ITensor}}(undef, 0)
    for i = 1:size(psis)[1]
        push!(arr, Array{ITensor}(undef, N))
    end
    return SpinDoob(N, c, s, zeros(N), zeros(N), zeros(N), psis, ones(size(psis[1])), 0,
                deepcopy(arr), deepcopy(arr))
end

function SpinDoob(N, c, s, psis, weights)
    arr = Array{Array{ITensor}}(undef, 0)
    for i = 1:size(psis)[1]
        push!(arr, Array{ITensor}(undef, N))
    end
    return SpinDoob(N, c, s, zeros(N), zeros(N), zeros(N), psis, weights, 0,
                deepcopy(arr), deepcopy(arr))
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


function updateLeftBlock!(template, idx)
    for i = 1:size(template.psis)[1]
        # Get the state and the psi block
        state = template.state[idx]
        A = template.psis[i][idx]
        s = siteind(template.psis[i], idx)
        projA = ITensor(s)
        projA[s[2 - state]] = 1

        # Fetch the previous block and multiply by A and projA
        block = idx == 1 ? 1 : copy(template.leftBlocks[i][idx-1])
        block *= A
        block *= dag(projA)
        template.leftBlocks[i][idx] = block
    end

end

function updateRightBlock!(template, idx)
    for i = 1:size(template.psis)[1]
        # Get the state and the psi block
        state = template.state[idx]
        A = template.psis[i][idx]
        s = siteind(template.psis[i], idx)
        projA = ITensor(s)
        projA[s[2 - state]] = 1

        # Fetch the previous block and multiply by A and projA
        block = idx == template.size ? 1 : copy(template.rightBlocks[i][idx+1])
        block *= A
        block *= dag(projA)
        template.rightBlocks[i][idx] = block
    end
end

function psiToLeft(template, state)
    # Count the number of ones and zeros in the state
    ones = sum(state)
    zeros = size(state)[1] - ones

    return (sqrt(template.c^(-1)))^(ones) * (sqrt((1-template.c)^(-1)))^(zeros)
end

function leftComponent(template, idx, current=false)
    prodTotal = 0.0
    for i = 1:size(template.psis)[1]
        # Get the psi component at the site and construct projector
        state = current ? template.state[idx] : 1 - template.state[idx]
        A = template.psis[i][idx]
        s = siteind(template.psis[i], idx)
        projA = ITensor(s)
        projA[s[2 - state]] = 1

        # Get the left and right blocks
        leftBlock = idx == 1 ? 1 : copy(template.leftBlocks[i][idx-1])
        rightBlock = idx == template.size ? 1 : copy(template.rightBlocks[i][idx+1])

        # Multiply together
        prod = leftBlock * A
        prod *= dag(projA)
        prod *= rightBlock
        prod = abs(real(scalar(prod)))
        prodTotal += template.weights[i] * prod
    end
    state = copy(template.state)
    state[idx] = current ? state[idx] : 1-state[idx]

    # Get conversion factor
    conversion = psiToLeft(template, state)
    if (prodTotal == 0.0)
        println(prodTotal)
    end
    return prodTotal*conversion
end

function leftComponent2(template, idx, current=false)
    # Get the state
    state = copy(template.state)

    if current==false
        state[idx] = 1 - state[idx]
    end
    mpsState = [x ? "Up" : "Dn" for x = state]

    sites = siteinds(template.psi)
    mpsState = productMPS(sites, mpsState)

    lc = dot(template.psi, mpsState)
    return lc
end


function updateTransitionRates!(template, blockLeft, blockRight)
    # Find the first and last sites which can flip
    idxLeft = findfirst([x != 0 for x = template.originalTransitionRates])
    idxRight = findlast([x != 0 for x = template.originalTransitionRates])

    # Build up the left and right blocks
    for i = blockLeft:idxRight-1
        updateLeftBlock!(template, i)
    end
    for i = blockRight:-1:idxLeft+1
        updateRightBlock!(template, i)
    end

    # Update the left component
    template.leftComponent = leftComponent(template, idxLeft, true)

    # Loop through each site and calculate the transition rates
    transitionRates = copy(template.originalTransitionRates)
    for idx = 1:template.size
        if transitionRates[idx] != 0
            left = leftComponent(template, idx)
            transitionRates[idx] *= exp(-template.s) * left / template.leftComponent
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
    updateTransitionRates!(template, idx, idx)
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
    # Add psis together
    psi = template.weights[1] * copy(template.psis[1])
    for i = 1:size(template.psis)[1]-1
        psi = add(psi, template.weights[i+1]*template.psis[i+1])
    end
    psi = psi * (1/sqrt(dot(psi, psi)))
    ITensors.orthogonalize!(psi, 1)

    # Sample
    state = 2 .- sample(psi)
    return convert(Array{Bool}, state)
end
initialState(template) = equilibriumState(template)
