#=
    Observers for the spin template
=#

# Occupations
"""
    occupationsMeasure(model)

Calculate the occupations of the system.
"""
function occupationsMeasure(model)
    return convert(Array{Float64}, copy(model.state))
end

"""
    occupationsObserver(N::Int64)

Create an observer to store occupations.
"""
function occupationsObserver(N)
    return observer("occupations", "configuration", occupationsMeasure, 0, zeros(N), zeros(N))
end

# Activity
"""
    activityMeasure(traj::trajectory)

Calculate the activity of a trajectory.
"""
function activityMeasure(traj)
    return size(traj.times)[1]-1
end


"""
    occupationsObserver()

Create an observer to store the activity.
"""
function activityObserver()
    return observer("activity", "trajectory", activityMeasure, 0, 0, 0)
end


# Time dependant occupations
"""
    timeOccupationsMeasure(traj::trajectory, times)

Calculate the time-occupations of a trajectory.
"""
function timeOccupationsMeasure(traj, times)
    # Loop through each time
    occupations = []
    for time = times
        # Find the last time smaller or eq to time
        idx = findlast([t <= time for t = traj.times])
        push!(occupations, convert(Array{Float64}, copy(traj.states[idx])))
    end

    return occupations
end


"""
    timeOccupationsObserver(N::Int64, times)

Create an observer to store time-occupations.
"""
function timeOccupationsObserver(N, times)
    measure(traj) = timeOccupationsMeasure(traj, times)
    initial = []
    for i = 1:size(times)[1]
        push!(initial, zeros(N))
    end
    return observer("time occupations", "trajectory", measure, 0, initial, initial)
end


# Auto correlation function
"""
    autoCorrelationMeasure(traj::trajectory, times)

Calculate the auto correlator of a trajectory.
"""
function autoCorrelationMeasure(traj, times)
    # Loop through each time
    correlations = []
    for time = times
        # Find the last time smaller or eq to time
        idx = findlast([t <= time for t = traj.times])
        push!(correlations, dot(traj.initial, traj.states[idx]))
    end

    return convert(Array{Float64}, correlations)
end


"""
    autoCorrelationObserver(N::Int64, times)

Create an observer to store auto correlations.
"""
function autoCorrelationObserver(times)
    measure(traj) = autoCorrelationMeasure(traj, times)
    initial = zeros(Float64, size(times))
    return observer("time occupations", "trajectory", measure, 0, initial, initial)
end


# Timescale of AC
"""
    timescaleMeasure(traj::trajectory)

Measure the time-integrated timescale.
"""
function timescaleMeasure(traj)
    timescale = 0.0
    # Loop through each state
    for i = 2:size(traj.states)[1]
        correl = dot(traj.initial, traj.states[i-1])
        timescale += correl * (traj.times[i] - traj.times[i-1])
    end
    correl = dot(traj.initial, traj.states[end])
    timescale += correl * (traj.maxTime - traj.times[end])

    return timescale
end


"""
    timescaleObserver(N::Int64, times)

Create an observer to store the time-integrated timescale.
"""
function timescaleObserver()
    measure(traj) = timescaleMeasure(traj)
    return observer("timescale", "trajectory", measure, 0, 0.0, 0.0)
end
