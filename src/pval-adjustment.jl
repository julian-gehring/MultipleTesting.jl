## p-value adjustment methods ##

type Bonferroni <: PValueAdjustmentMethod
end

adjust(pvals, method::Bonferroni) = bonferroni(pvals)
adjust{T}(pvals::Vector{T}, ws::MultipleTestingWeights{T}, method::Bonferroni) = bonferroni(pvals./ws)

function bonferroni{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    return min(pValues * length(pValues), 1.)
end

type BenjaminiHochberg <: PValueAdjustmentMethod
end

adjust(pvals, method::BenjaminiHochberg) = benjamini_hochberg(pvals)
adjust{T}(pvals::Vector{T}, ws::MultipleTestingWeights{T}, method::BenjaminiHochberg) = benjamini_hochberg(pvals./ws)

function benjamini_hochberg{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, bejamini_hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

bejamini_hochberg_multiplier(i::Int, n::Int) = n/(n-i)

type BenjaminiHochbergAdaptive <: PValueAdjustmentMethod
    pi0estimator::Pi0Estimator
end

## default to BenjaminiHochberg
BenjaminiHochbergAdaptive(π0::AbstractFloat) = BenjaminiHochbergAdaptive(Oracle(π0))

BenjaminiHochbergAdaptive() = BenjaminiHochbergAdaptive(1.0)

function adjust(pvals, method::BenjaminiHochbergAdaptive)
  π0 = estimate_pi0(pvals, method.pi0estimator)
  benjamini_hochberg(pvals, π0)
end

function adjust{T}(pvals::Vector{T}, ws::MultipleTestingWeights{T}, ::BenjaminiHochbergAdaptive)
  π0 = estimate_pi0(pvals, ws, method.pi0estimator)
  benjamini_hochberg(pvals, ws, π0)
end

function benjamini_hochberg{T<:AbstractFloat}(pValues::Vector{T}, pi0::T)
    validPValues([pi0])
    benjamini_hochberg(pValues) .* pi0
end


type BenjaminiYekutieli <: PValueAdjustmentMethod
end

adjust(pvals, method::BenjaminiYekutieli) = benjamini_yekutieli(pvals)

function benjamini_yekutieli{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, benjamini_yekutieli_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

function benjamini_yekutieli_multiplier(i::Int, n::Int)
    c = sum([1/i for i in 1:n])
    return ((n*c)/(n-i))
end

type BenjaminiLiu <: PValueAdjustmentMethod
end

adjust(pvals, method::BenjaminiLiu) = benjamini_liu(pvals)

function benjamini_liu{T<:AbstractFloat}(pvalues::Vector{T})
    validPValues(pvalues)
    n = length(pvalues)
    if n <= 1
        return pvalues
    end
    sortedIndexes, originalOrder = reorder(pvalues)
    sortedPValues = pvalues[sortedIndexes]
    general_stepdown!(sortedPValues, benjamini_liu_step, n)
    min(sortedPValues[originalOrder], 1.)
end

function benjamini_liu_step{T<:AbstractFloat}(p::T, i::Int, n::Int)
    # a bit more involved because cutoffs at significance α have the form:
    # P_(i) <= 1- [1 - min(1, m/(m-i+1)α)]^{1/(m-i+1)}
    adjp = (1-(1-p)^(n-i+1))*(n-i+1)/n
    if adjp*n/(n-i+1) > 1
        return 0.0
    else
        return adjp
    end
end

type Hochberg <: PValueAdjustmentMethod
end

adjust(pvals, method::Hochberg) = hochberg(pvals)

function hochberg{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, hochberg_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

hochberg_multiplier(i::Int, n::Int) = (i+1)


type Holm <: PValueAdjustmentMethod
end

adjust(pvals, method::Holm) = holm(pvals)

function holm{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, holm_multiplier, n)
    min(sortedPValues[originalOrder], 1.)
end

holm_multiplier(i::Int, n::Int) = (n-i+1)


type Hommel <: PValueAdjustmentMethod
end

adjust(pvals, method::Hommel) = hommel(pvals)

function hommel{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    n = length(pValues)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    q = fill(minimum(n .* pValues./[1:n; ]), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* sortedPValues[i2]./([2:j; ]))
        q[ij] = min(j .* sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max(pa, q)
    end
    max(pa, sortedPValues)[originalOrder]
end


type Sidak <: PValueAdjustmentMethod
end

adjust(pvals, method::Sidak) = sidak(pvals)

function sidak{T<:AbstractFloat}(pValues::Vector{T})
    validPValues(pValues)
    return min(1-(1-pValues).^length(pValues), 1.)
end

type ForwardStop <: PValueAdjustmentMethod
end

adjust(pvals, method::ForwardStop) = forwardstop(pvals)

function forwardstop{T<:AbstractFloat}(pvalues::Vector{T})
    validPValues(pvalues)
    n = length(pvalues)
    logsums = - cumsum(log(1-pvalues))
    stepup!(logsums, forwardstop_multiplier, n)
    max(min(logsums, 1.),0.)
end

forwardstop_multiplier(i::Int, n::Int) = 1/(n-i)
