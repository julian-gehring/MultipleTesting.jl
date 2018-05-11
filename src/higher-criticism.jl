## Higher criticism

"""
Higher criticism scores

**Examples**

```jldoctest
HigherCriticismScores()
```

**References**

Donoho, D., and Jin, J. (2008). Higher criticism thresholding: Optimal feature
selection when useful features are rare and weak. PNAS 105, 14790–14795.

Klaus, B., and Strimmer, K. (2013). Signal identification for rare and weak
features: higher criticism or false discovery rates? Biostatistics 14, 129–143.

"""
struct HigherCriticismScores
end

function estimate(pValues::PValues{T}, method::HigherCriticismScores) where T<:AbstractFloat
    n = length(pValues)
    F = (n+1 .- competerank(-pValues)) ./ n  # ECDF
    denom = F .* (one(T) .- F) ./ n
    # avoid denominator of 0 for last value
    idx0 = denom .== 0
    denom[idx0] = minimum(denom[.!idx0]) .+ eps()  # conservative
    hcs = abs.(F .- pValues) ./ sqrt.(denom)
    return hcs
end


"""
Higher criticism threshold

**Examples**

```jldoctest
HigherCriticismThreshold()
```

**References**

Donoho, D., and Jin, J. (2008). Higher criticism thresholding: Optimal feature
selection when useful features are rare and weak. PNAS 105, 14790–14795.

Klaus, B., and Strimmer, K. (2013). Signal identification for rare and weak
features: higher criticism or false discovery rates? Biostatistics 14, 129–143.

"""
struct HigherCriticismThreshold
end

function estimate(pValues::PValues{T}, method::HigherCriticismThreshold) where T<:AbstractFloat
    idx_hcv = indmax(estimate(pValues, HigherCriticismScores()))
    return pValues[idx_hcv]
end
