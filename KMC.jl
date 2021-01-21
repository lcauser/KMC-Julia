#=
    Code to create Kinetic Monte Carlo simulations.

    Author: Luke Causer
=#

function randomChoice(weights)
    cs = cumsum(weights)
    r = rand(Float64) * cs[end]
    idx = 1
    for val in cs
        if val < r
            idx += 1
        end
    end

    return idx
end


function simulation(model, initial, maxTime)
    # Set the initial state and variables
    setState!(model, initial)
    updateTransitionRates!(model)
    time = 0
    idxs = []
    times = []

    # Loop through until the max time has passed
    while time < maxTime
        # Find a configuration to transition into and what time
        idx = randomChoice(model.transitionRates)
        escapeRate = sum(model.transitionRates)
        transitionTime = -log(rand(Float64)) / escapeRate
        time += transitionTime

        if time < maxTime
            # Do the transition
            transition!(model, idx)
            updateTransitionRates!(model, idx)

            # Store information
            append!(idxs, idx)
            append!(times, transitionTime)
        end
    end

    return idxs, times
end


function reconstruct(model, initial, idxs, times)
    # Preallocate memory
    sz = size(idxs)[1]
    realTimes = Array{Float64}(undef, sz+1)
    states = Array{typeof(model.state)}(undef, sz+1)

    # Set the initial state and variables
    time = 0
    setState!(model, initial)

    # Set starting information
    realTimes[1] = time
    states[1] = initial

    # Loop through changing the configuration
    for i = 1:sz
        transition!(model, idxs[i])
        time += times[i]
        realTimes[i+1] = time
        states[i+1] = copy(model.state)
    end

    return states, realTimes

end
