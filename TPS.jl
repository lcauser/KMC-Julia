#=
    Run TPS simulations
=#

include("KMC.jl")

"""
    propose(model, traj::trajectory, configObservers, portion)

Propose a new trajectory for the model given a current trajectory and the
measured observables with the shifting method. You must specify which portion
to keep.
"""
function propose(model, traj, configObservers, portion, rate=1.0)
    # How to cut the trajectory
    splitTime = rand(Float64)*traj.maxTime
    portion = rand(Bool) # 0 first 1 second
    #portion = false
    splitTime = portion ? splitTime*rate : traj.maxTime-rate*splitTime

    # Split the trajectory and just keep one
    pTraj1, pTraj2 = splitTrajectory(traj, splitTime)
    partialTraj1 = portion ? pTraj2 : pTraj1

    # Propose a new partial trajectory
    pMaxTime = portion ? splitTime : traj.maxTime - splitTime
    pInitial = portion ? partialTraj1.states[end] : partialTraj1.initial
    partialTraj2 = simulation(model, pInitial, pMaxTime, configObservers)
    partialTraj2 = portion ? partialTraj2 : reverseTrajectory(partialTraj2)

    # Join the two trajectories
    newTrajectory = portion ? joinTrajectory(partialTraj1, partialTraj2) :
                              joinTrajectory(partialTraj2, partialTraj1)

    return newTrajectory

end

function gidx(lst, val)
    idx = 1
    for item = lst
        if item < val
            idx += 1
        end
    end

    if idx > size(lst)[1]
        idx = 0
    end

    return idx
end


"""
    metropolis(crit1::Float64, crit2::Float64)

Find the acceptance rate according to metropolis criterion.
"""
function metropolis(crit1, crit2)
    r = rand(Float64)
    return r < exp(crit2 - crit1) ? true : false
end


"""
    TPS(model, traj::Trajectory, trajMeasures, numSims::Int, maxTime:Float64,
    observers, warmUp::Int=0)

Run TPS simulations for the model. You must define the criterion function.
"""
function TPS(model, traj, trajMeasures, numSims, maxTime, observers, warmUp=0,
             quiet=false, flipRate = 0.5, ratio = 1)
    # Get initial probability
    crit = criterion(traj, trajMeasures)

    # Split observers into configuration based and trajectory based
    configObservers = []
    trajectoryObservers = []
    for i = 1:size(observers)[1]
        if observers[i].type == "configuration"
            push!(configObservers, i)
        else
            push!(trajectoryObservers, i)
        end
    end

    acceptance = 0
    portion = true
    for sim = 1:(numSims+warmUp)
        # Propose a new trajectory and calculate observables
        newTraj = propose(model, traj, observers[configObservers], portion, ratio)
        newTrajMeasures = []
        for idx in trajectoryObservers
            observers[idx].currentMeasure = measureObserver(observers[idx], newTraj)
            push!(newTrajMeasures, observers[idx].currentMeasure)
        end
        newCrit = criterion(newTraj, newTrajMeasures)

        # Accept or reject
        accept = metropolis(crit, newCrit)
        if accept
            traj = newTraj
            trajMeasures = newTrajMeasures
            crit = newCrit
            acceptance += 1
        end

        r = rand()
        if r < flipRate
            portion = portion ? false : true
        end

        # Measure observers
        if sim > warmUp
            i = 1
            for idx in configObservers
                measure = timeIntegrate(traj.times, traj.observers[i], maxTime)
                updateObserver!(observers[idx], measure)
                i += 1
            end
            for idx in trajectoryObservers
                updateObserver!(observers[idx], trajMeasures[idx])
            end
        end

        if !quiet
            println(string("Simulation ", sim, "/", numSims+warmUp, " completed."))
        end
    end
    if !quiet
        println(string("Acceptance rate ", acceptance/(numSims+warmUp), " completed."))
    end

    return traj, trajMeasures, acceptance
end
