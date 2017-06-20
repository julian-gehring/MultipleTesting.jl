## p-value adjustment methods ##

# promotion from float vectors to PValues type

adjust{T<:AbstractFloat, M<:PValueAdjustment}(pvals::Vector{T}, method::M) = adjust(PValues(pvals), method)

adjust{T<:AbstractFloat, M<:PValueAdjustment}(pvals::Vector{T}, n::Integer, method::M) = adjust(PValues(pvals), n, method)


# Bonferroni

immutable Bonferroni <: PValueAdjustment
end

adjust{P<:GeneralizedPValues}(pvals::P, method::Bonferroni) = adjust(pvals, length(pvals), method)

adjust{P<:GeneralizedPValues}(pvals::P, n::Integer, method::Bonferroni) = bonferroni(pvals, n)

function adjust(pvals::PValues, weights::PriorityWeights, method::Bonferroni)
    adjust( WeightedPValues(pvals, weights), method)
end

function bonferroni{P<:GeneralizedPValues}(pValues::P, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    return min.(pValues * n, 1)
end


# Benjamini-Hochberg

immutable BenjaminiHochberg <: PValueAdjustment
end

adjust{P<:GeneralizedPValues}(pvals::P, method::BenjaminiHochberg) = adjust(pvals, length(pvals), method)

adjust{P<:GeneralizedPValues}(pvals::P, n::Integer, method::BenjaminiHochberg) = benjamini_hochberg(pvals, n)

function adjust(pvals::PValues, weights::PriorityWeights, method::BenjaminiHochberg)
    adjust( WeightedPValues(pvals, weights), method)
end

function benjamini_hochberg{P<:GeneralizedPValues}(pValues::P, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, bejamini_hochberg_step, k, n)
    return min.(sortedPValues[originalOrder], 1)
end

bejamini_hochberg_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * n/(k-i)


# Benjamini-Hochberg Adaptive

immutable BenjaminiHochbergAdaptive <: PValueAdjustment
    pi0estimator::Pi0Estimator
end

## default to BenjaminiHochberg
BenjaminiHochbergAdaptive(π0::AbstractFloat) = BenjaminiHochbergAdaptive(Oracle(π0))

BenjaminiHochbergAdaptive() = BenjaminiHochbergAdaptive(1.0)

adjust(pvals::PValues, method::BenjaminiHochbergAdaptive) = adjust(pvals, length(pvals), method)

function adjust(pvals::PValues, n::Integer, method::BenjaminiHochbergAdaptive)
    π0 = estimate_pi0(pvals, method.pi0estimator)
    return benjamini_hochberg(pvals, n) * π0
end


# Benjamini-Yekutieli

immutable BenjaminiYekutieli <: PValueAdjustment
end

adjust(pvals::PValues, method::BenjaminiYekutieli) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::BenjaminiYekutieli) = benjamini_yekutieli(pvals, n)

function benjamini_yekutieli(pValues::PValues, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, benjamini_yekutieli_step, k, n)
    return min.(sortedPValues[originalOrder], 1)
end

benjamini_yekutieli_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * harmonic_number(n) * n/(k-i)


# Benjamini-Liu

immutable BenjaminiLiu <: PValueAdjustment
end

adjust(pvals::PValues, method::BenjaminiLiu) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::BenjaminiLiu) = benjamini_liu(pvals, n)

function benjamini_liu(pValues::PValues, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, benjamini_liu_step, k, n)
    return min.(sortedPValues[originalOrder], 1)
end

function benjamini_liu_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer)
    # a bit more involved because cutoffs at significance α have the form:
    # P_(i) <= 1- [1 - min(1, m/(m-i+1)α)]^{1/(m-i+1)}
    s = n-i+1
    return (1 - (1-p)^s) * s / n
end


# Hochberg

immutable Hochberg <: PValueAdjustment
end

adjust(pvals::PValues, method::Hochberg) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::Hochberg) = hochberg(pvals, n)

function hochberg(pValues::PValues, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepup!(sortedPValues, hochberg_step, k, n)
    return min.(sortedPValues[originalOrder], 1)
end

hochberg_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * (n-k+i+1)


# Holm

immutable Holm <: PValueAdjustment
end

adjust(pvals::PValues, method::Holm) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::Holm) = holm(pvals, n)

function holm(pValues::PValues, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    stepdown!(sortedPValues, holm_step, k, n)
    return min.(sortedPValues[originalOrder], 1)
end

holm_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * (n-i+1)


# Hommel

immutable Hommel <: PValueAdjustment
end

adjust(pvals::PValues, method::Hommel) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::Hommel) = hommel(pvals, n)

function hommel(pValues::PValues, n::Integer)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    pValues = vcat(pValues, ones(n-k))  # TODO avoid sorting of ones
    sortedIndexes, originalOrder = reorder(pValues)
    sortedPValues = pValues[sortedIndexes]
    q = fill(minimum(n .* pValues./(1:n)), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* sortedPValues[i2]./((2:j)))
        q[ij] = min.(j .* sortedPValues[ij], q1)
        q[i2] = q[n-j+1]
        pa = max.(pa, q)
    end
    return max.(pa, sortedPValues)[originalOrder[1:k]]
end


# Sidak

immutable Sidak <: PValueAdjustment
end

adjust(pvals::PValues, method::Sidak) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::Sidak) = sidak(pvals, n)

function sidak(pValues::PValues, n::Integer)
    check_number_tests(length(pValues), n)
    return min.(1-(1-pValues).^n, 1)
end


# Forward Stop

immutable ForwardStop <: PValueAdjustment
end

adjust(pvals::PValues, method::ForwardStop) = adjust(pvals, length(pvals), method)

adjust(pvals::PValues, n::Integer, method::ForwardStop) = forwardstop(pvals, n)

function forwardstop(pvalues::PValues, n::Integer)
    k = length(pvalues)
    check_number_tests(k, n)
    logsums = -cumsum(log.(1-pvalues))
    stepup!(logsums, forwardstop_step, k, n)
    return max.(min.(logsums, 1), 0)
end

forwardstop_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p * 1/(k-i)


# Barber-Candès
immutable BarberCandes <: PValueAdjustment
end

function adjust(pvals::PValues, method::BarberCandes)
    n = length(pvals)

    if n <= 1
        return ones(pvals) # unlike other p-adjust methods
    end

    sorted_indexes, original_order = reorder(pvals)
    estimated_fdrs = pvals[sorted_indexes]

    Rt = 1 # current number of discoveries
    Vt = 1 # estimated false discoveries at t

    left_pv = estimated_fdrs[1]
    right_pv = estimated_fdrs[n]

    while (left_pv < 0.5)
        while ((1-right_pv <= left_pv) & (right_pv >= 0.5) & (Vt+Rt <= n))
            estimated_fdrs[n-Vt+1] = 1.0;
            right_pv = estimated_fdrs[n-Vt]
            Vt +=1
        end
        estimated_fdrs[Rt] = Vt/Rt
        Rt += 1
        left_pv = estimated_fdrs[Rt]
    end

    while ((right_pv >= 0.5) & (Vt + Rt <= n+1))
      estimated_fdrs[n-Vt+1] = 1.0;
      right_pv = (Vt + Rt <= n)?estimated_fdrs[n-Vt]:0.0
      Vt +=1
    end

    stepup!(estimated_fdrs, identity_step, n, n) # just monotonize, no multiplier needed
    return min.(estimated_fdrs[original_order], 1)
end

# as test, inefficient implementation
function barber_candes_brute_force{T<:AbstractFloat}(pvals::AbstractVector{T})
    n = length(pvals)
    sorted_indexes, original_order = reorder(pvals)
    sorted_pvals = pvals[sorted_indexes]
    estimated_fdrs = ones(pvals)
    for (i,pv) in enumerate(sorted_pvals)
        if pv >= 0.5
            break
        else
            estimated_fdrs[i] = (sum(1-pvals .<= pv)+1)/i
        end
    end
    stepup!(estimated_fdrs, identity_step, n, n) # just monotonize, no multiplier needed
    return min.(estimated_fdrs[original_order], 1)
end

identity_step(p::AbstractFloat, i::Integer, k::Integer, n::Integer) = p


## internal ##

# step-up / step-down

function stepup!{T<:AbstractFloat}(sortedPValues::AbstractVector{T}, stepfun::Function, k::Integer, n::Integer)
    sortedPValues[k] = stepfun(sortedPValues[k], 0, k, n)
    for i in 1:(k-1)
        sortedPValues[k-i] = min.(sortedPValues[k-i+1], stepfun(sortedPValues[k-i], i, k, n))
    end
    return sortedPValues
end

function stepdown!{T<:AbstractFloat}(sortedPValues::AbstractVector{T}, stepfun::Function, k::Integer, n::Integer)
    sortedPValues[1] = stepfun(sortedPValues[1], 1, k, n)
    for i in 2:k
        sortedPValues[i] = max.(sortedPValues[i-1], stepfun(sortedPValues[i], i, k, n))
    end
    return sortedPValues
end


function check_number_tests(k::Integer, n::Integer)
    if k > n
        msg = "Number of total tests ($n) is smaller than number of p-values ($k)"
        throw(ArgumentError(msg))
    end
end


harmonic_number(n::Integer) =  digamma(n+1) + γ
