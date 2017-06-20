### Combination methods for p-values ###

combine{T<:AbstractFloat, M<:PValueCombination}(pValues::AbstractVector{T}, method::M) = combine(PValues(pValues), method)

combine{T<:AbstractFloat, M<:PValueCombination}(pValues::AbstractVector{T}, weights::PriorityWeights, method::M) = combine(PValues(pValues), weights, method)

combine{T<:AbstractFloat, R<:Real, M<:PValueCombination}(pValues::AbstractVector{T}, weights::AbstractVector{R}, method::M) = combine(PValues(pValues), weights, method)


## Fisher combination ##

immutable FisherCombination <: PValueCombination
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::FisherCombination)
    fisher_combination(pValues)
end

function fisher_combination{T<:AbstractFloat}(pValues::PValues{T})
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    if minimum(pValues) == 0.0
        return NaN
    end
    x = -2 * sum(log.(pValues))
    p = ccdf(Chisq(2n), x)
    return p
end


## Logit combination ##

immutable LogitCombination <: PValueCombination
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::LogitCombination)
    logit_combination(pValues)
end

function logit_combination{T<:AbstractFloat}(pValues::PValues{T})
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    c = sqrt( (5n+2)*n*pi^2 / ((5n+4)*3) )
    x = -sum(log.(pValues./(1-pValues))) / c
    p = ccdf(TDist(5n+4), x)
    return p
end


## Stouffer combination ##

immutable StoufferCombination <: PValueCombination
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::StoufferCombination)
    stouffer_combination(pValues)
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, weights::Vector{T}, method::StoufferCombination)
    stouffer_combination(pValues, weights)
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, weights::PriorityWeights, method::StoufferCombination)
    stouffer_combination(pValues, values(weights))
end

function stouffer_combination{T<:AbstractFloat}(pValues::PValues{T})
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    z = cquantile(Normal(), pValues)
    z = sum(z) ./ sqrt(n)
    p = ccdf(Normal(), z)
    return p
end

function stouffer_combination{T<:AbstractFloat}(pValues::PValues{T}, weights::Vector{T})
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pmin, pmax = extrema(pValues)
    if pmin == 0.0 || pmax == 1.0
        return NaN
    end
    z = cquantile(Normal(), pValues) .* weights
    z = sum(z) ./ sqrt(sum(weights.^2))
    p = ccdf(Normal(), z)
    return p
end


## Tippett combination ##

immutable TippettCombination <: PValueCombination
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::TippettCombination)
    tippett_combination(pValues)
end

function tippett_combination{T<:AbstractFloat}(pValues::PValues{T})
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    p = 1.0 - (1.0 - minimum(pValues))^n
    return p
end


## Simes combination ##

immutable SimesCombination <: PValueCombination
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::SimesCombination)
    simes_combination(pValues)
end

function simes_combination{T<:AbstractFloat}(pValues::PValues{T})
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    pValues = sort(pValues)  # faster than `sortperm`
    p = n * minimum(pValues./(1:n))
    return p
end


## Wilkinson combination ##

immutable WilkinsonCombination <: PValueCombination
    rank::Integer

    function WilkinsonCombination(rank)
        if rank < 1
            throw(ArgumentError("Rank must be positive."))
        end
        return new(rank)
    end
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::WilkinsonCombination)
    wilkinson_combination(pValues, method.rank)
end

function wilkinson_combination{T<:AbstractFloat}(pValues::PValues{T}, rank::Integer)
    n = length(pValues)
    if rank < 1 || rank > n
        throw(ArgumentError("Rank must be in 1,..,$(n)"))
    end
    if n == 1
        return pValues[1]
    end
    p_rank = sort(pValues)[rank]
    p = cdf(Beta(rank, n-rank+1), p_rank)
    return p
end


## Generalised minimum combination ##

immutable MinimumCombination <: PValueCombination
    method::PValueAdjustment
end

function combine{T<:AbstractFloat}(pValues::PValues{T}, method::MinimumCombination)
    minimum_combination(pValues, method.method)
end

function minimum_combination{T<:AbstractFloat}(pValues::PValues{T}, pAdjustMethod::PValueAdjustment)
    n = length(pValues)
    if n == 1
        return pValues[1]
    end
    padj = adjust(pValues, pAdjustMethod)
    p = minimum(padj)
    return p
end
