#=
    Observers for the doob dynamics
=#

# Left components
"""
    leftMeasure(model)

Calculate the left-eigenvector component
"""
function leftMeasure(model)
    return copy(model.leftComponent)
end


"""
    leftObserver())

Create an observer to left components.
"""
function leftObserver()
    return observer("left", "configuration", leftMeasure, 0, 0, 0)
end


# Escape rate
"""
    escapeRateMeasure(model)

Calculate the escape rate of a model.
"""
function escapeRateMeasure(model)
    return copy(sum(model.transitionRates))
end


"""
    escapeRateObserver())

Create an observer to store the escape rates.
"""
function escapeRateObserver()
    return observer("ER", "configuration", escapeRateMeasure, 0, 0, 0)
end


# Original escape rate
"""
    originalEscapeRateMeasure(model)

Calculate the escape rate of the original dynamics.
"""
function originalEscapeRateMeasure(model)
    return copy(sum(model.originalTransitionRates))
end


"""
    originalEscapeRateObserver())

Create an observer to store the original dynamics escape rates.
"""
function originalEscapeRateObserver()
    return observer("OER", "configuration", originalEscapeRateMeasure, 0, 0, 0)
end
