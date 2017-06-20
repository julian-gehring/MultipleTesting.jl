type IsotonicRegression{T<:Real} <: RegressionModel
    Xorig::Vector{T} # original grid points, but sorted
    Yorig::Vector{T} # original values, sorted according to x
    Xrank::Vector{Int}
    Yfitted::Vector{T}
end

function fit{T<:Real}(::Type{IsotonicRegression}, X::AbstractVector{T}, Y::AbstractVector{T},
        direction::Symbol=:isotonic)
    Xrank = ordinalrank(X) # maybe speed this up when already sorted
    Xperm = sortperm(X)
    Xorig = X[Xperm]
    Yorig = Y[Xperm]

    if direction == :isotonic
        Yfitted = MultipleTesting.isotonic_regression(Yorig)
    elseif direction == :antitonic
        Yfitted = -MultipleTesting.isotonic_regression(-Yorig)
    end

    IsotonicRegression{T}(Xorig, Yorig, Xrank, Yfitted)
end

function predict(iso::IsotonicRegression)
    iso.Yfitted[iso.Xrank]
end

function predict{T<:Real}(iso::IsotonicRegression{T}, x::T)
    n = length(iso.Xorig)
    idxl = searchsortedlast(iso.Xorig, x)
    idxr = idxl + 1
    if idxl == 0 # constant extrapolation to the left
        val = iso.Yfitted[1]
    elseif idxl == n # constant extrapolation to the right
        val =  iso.Yfitted[n]
    else # linear interpolation
        xl = iso.Xorig[idxl]
        xr = iso.Xorig[idxr]
        valr = iso.Yfitted[idxl]
        vall = iso.Yfitted[idxr]
        h = xr-xl
        val = ((x-xl)*valr + (xr-x)*vall)/h
    end
    val
end

function predict{T<:Real}(iso::IsotonicRegression{T}, xs::AbstractVector{T})
    g = zeros(xs)
    for (i,x) in enumerate(xs)
        # can be made more efficient as in StatsBase ECDF code
        g[i] = predict(iso, x)
    end
    g
end
