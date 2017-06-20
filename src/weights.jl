immutable PriorityWeights{S<:AbstractFloat} <: AbstractWeights{S, S, Vector{S}}
    values::Vector{S}
    sum::S
end

PriorityWeights{S<:AbstractFloat}(vs::Vector{S}) = PriorityWeights{S}(vs/sum(vs)*length(vs), length(vs))


immutable WeightedPValues{T<:AbstractFloat} <: GeneralizedPValues{T}
    pvalues::PValues{T}
    weights::PriorityWeights{T} #TODO check these are of the same size
end


Base.size(wpv::WeightedPValues) = (length(wpv.pvalues), )
@compat Base.IndexStyle{T<:WeightedPValues}(::Type{T}) = IndexLinear()

function Base.getindex{T<:AbstractFloat}(wpv::WeightedPValues{T}, i::Integer)
    pvi = wpv.pvalues[i]
    wi = wpv.weights[i]
    if pvi == zero(T)
      wpv = zero(T)
    elseif wi == zero(T)
      wpv = one(T)
    else
      wpv = pvi/wi
    end
    min(wpv, one(T))
end

Base.setindex!(wpv::WeightedPValues, x::AbstractFloat, i::Integer) =
    throw(ErrorException("Modification of values is not permitted"))


Base.length(wpv::WeightedPValues) = length(wpv.pvalues)
