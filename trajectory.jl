#=
    Trajectories contain all the information needed to reconstruct a
    stochastic trajectory.

    Author: Luke Causer
=#


"""
    trajectory

A stochastic trajectory. Contains the initial state, intermediate states and
the transition times.
"""
mutable struct trajectory{}
    initial
    states
    times::Array{Float64}
    maxTime::Float64
    observers::Array{Any}
end


"""
    trajectory(model, initial, idxs, times, maxTime, observers)

Construct a trajectory with minimal information needed (initial, transition
identifiers, transition times, the max time, and observer information.)
"""
function trajectory(model, initial, idxs, times, maxTime, observers)
    configs, realTimes = reconstruct(model, initial, idxs, times)
    return trajectory(initial, configs, realTimes, maxTime, observers)
end


"""
    trajectory(model, initial, idxs, times, maxTime)

Construct a trajectory with minimal information needed (initial, transition
identifiers, transition times and the max time)
"""
function trajectory(model, initial, idxs, times, maxTime)
    configs, realTimes = reconstruct(model, initial, idxs, times)
    return trajectory(initial, configs, realTimes, maxTime, [])
end


"""
    splitTrajectory(traj::trajectory, splitTime::Float64)

Splits a trajectory into two at the given time.
"""
function splitTrajectory(traj, splitTime::Float64)
    # Find the index where to split
    idx = gidx(traj.times, splitTime)

    # First partial trajectory
    endIdx = idx == 0 ? size(traj.times)[1] : max(1, idx - 1)
    partialInitial1 = copy(traj.initial)
    partialTimes1 = copy(traj.times[1:endIdx])
    partialStates1 = copy(traj.states[1:endIdx])
    partialObservers1 = []
    for observer in traj.observers
        push!(partialObservers1, copy(observer[1:endIdx]))
    end
    partialMaxTime1 = splitTime
    partialTrajectory1 = trajectory(partialInitial1, partialStates1,
                                    partialTimes1, partialMaxTime1,
                                    partialObservers1)

    # Second partial trajectory
    startIdx = idx == 0 ? size(traj.times)[1] : max(1, idx - 1)
    partialInitial2 = copy(traj.states[startIdx])
    partialTimes2 = copy(traj.times[startIdx:end]) .- splitTime
    partialTimes2[1] = 0
    partialStates2 = copy(traj.states[startIdx:end])
    partialObservers2 = []
    for observer in traj.observers
        push!(partialObservers2, copy(observer[startIdx:end]))
    end
    partialMaxTime2 = copy(traj.maxTime) - splitTime
    partialTrajectory2 = trajectory(partialInitial2, partialStates2,
                                    partialTimes2, partialMaxTime2,
                                    partialObservers2)

    return partialTrajectory1, partialTrajectory2
end


"""
    reverseTrajectory(traj::trajectory)

Reverses a trajectory
"""
function reverseTrajectory(traj)
    states = reverse(traj.states)
    initial = states[1]
    observers = []
    for observer in traj.observers
        push!(observers, reverse(observer))
    end

    times = [0.0]
    sz = size(traj.times)[1]
    for idx in 1:sz-1
        push!(times, traj.maxTime - traj.times[sz+1-idx])
    end

    return trajectory(initial, states, times, traj.maxTime, observers)
end


"""
    joinTrajectory(traj1::trajectory, traj2::trajectory)

Merges two trajectories.
"""
function joinTrajectory(traj1, traj2)
    initial = copy(traj1.initial)
    maxTime = traj1.maxTime + traj2.maxTime
    if size(traj2.states)[1] > 1
        states = vcat(traj1.states, traj2.states[2:end])
        times = vcat(traj1.times, traj2.times[2:end] .+ traj1.maxTime)
        observers = []
        for i = 1:size(traj1.observers)[1]
            push!(observers, vcat(traj1.observers[i], traj2.observers[i][2:end]))
        end
    else
        states = copy(traj1.states)
        times = copy(traj1.times)
        observers = []
        for i = 1:size(traj1.observers)[1]
            push!(observers, copy(traj1.observers[i]))
        end
    end

    return trajectory(initial, states, times, maxTime, observers)
end
