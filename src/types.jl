## abstract types ##

@compat abstract type Pi0Estimator end

@compat abstract type Pi0Fit end

@compat abstract type Pi0RegressionEstimator end

@compat abstract type Pi0RegressionFit end

@compat abstract type PValueAdjustment end

@compat abstract type PValueCombination end


## concrete types ##

## PValues

@compat abstract type GeneralizedPValues{T<:AbstractFloat} <:AbstractVector{T} end

immutable PValues{T<:AbstractFloat} <: GeneralizedPValues{T}
    values::AbstractVector{T}
    min::T
    max::T

    function(::Type{PValues}){T}(values::AbstractVector{T})
        min, max = extrema(values)
        if min < 0.0 || max > 1.0
            throw(DomainError())
        end
        new{T}(copy(values), min, max)
    end
end

Base.convert{T<:AbstractFloat}(::Type{PValues}, x::AbstractVector{T}) = PValues(x)

Base.size(pv::PValues) = (length(pv.values), )
@compat Base.IndexStyle{T<:PValues}(::Type{T}) = IndexLinear()
Base.getindex(pv::PValues, i::Integer) = pv.values[i]
Base.setindex!(pv::PValues, x::AbstractFloat, i::Integer) =
    throw(ErrorException("Modification of values is not permitted"))

Base.values(pv::PValues) = pv.values
Base.minimum(pv::PValues) = pv.min
Base.maximum(pv::PValues) = pv.max
Base.extrema(pv::PValues) = (minimum(pv), maximum(pv))
Base.length(pv::PValues) = length(pv.values)
