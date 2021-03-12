#=
    The observer struct is the base for observers
=#
"""
    observer

An observer (configuration or trajectory observable) measures an observable on
stochastic trajectories and stores the cumulative sum.
"""
mutable struct observer{}
    name::String
    type::String
    func
    numSims::Int
    measure
    currentMeasure
end


"""
    measureObserver(observer, model)

Measure an configuration observable given the model.
"""
function measureObserver(observer, model)
    return observer.func(model)
end


"""
    measureObserver(observer, traj::trajectory)

Measure a trajectory observerable on a trajectory.
"""
function measureObserver(observer, traj::trajectory)
    return observer.func(traj)
end


"""
    updateObserver!(observer, measure)

Update the cumulative measure of an observable.
"""
function updateObserver!(observer, measure)
    observer.numSims += 1
    observer.measure += measure
end


"""
    storeMeasure!(observer, measure)

Store the current measure of an observable.
"""
function storeMeasure!(observer, measure)
    observer.currentMeasure = measure
end


"""
    getMeasure!(observer)

Get the average measure of an observable.
"""
function getMeasure(observer)
    return observer.measure / observer.numSims
end
