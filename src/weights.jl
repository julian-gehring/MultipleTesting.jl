immutable MultipleTestingWeights{T <: AbstractFloat} <: AbstractVector{T}
    ws::Vector{T}
    function MultipleTestingWeights(ws::Vector{T})
      minimum(ws) >= 0 || throw(DomainError())
      s = mean(ws)
      sum(ws) >0 || throw(DomainError())
      new(ws./s)
    end
end

MultipleTestingWeights{T<:AbstractFloat}(ws::Vector{T}) = MultipleTestingWeights{T}(ws)

# make it work as array
Base.size{T}(ws::MultipleTestingWeights{T}) = size(ws.ws)
Base.length{T}(ws::MultipleTestingWeights{T}) = length(ws.ws)
Base.linearindexing{T}(::Type{MultipleTestingWeights{T}}) = Base.LinearFast()
Base.getindex{T}(ws::MultipleTesting.MultipleTestingWeights{T}, i::Int) = ws.ws[i]


@inline function mtpdiv{T}(p::T, w::T)
    if p == 0.0
        wp = 0.0
    else
        wp = min(one(T), p/w)
    end
    wp
end

function ./{T}(pvalue::Vector{T}, ws::MultipleTestingWeights{T})
    length(pvalue) == length(ws) || throw(DomainError())
    wpv = copy(pvalue) # make this in-place?
    @inbounds for i=1:length(pvalue)
        wpv[i] = mtpdiv(pvalue[i],ws[i])
    end
    wpv
end
