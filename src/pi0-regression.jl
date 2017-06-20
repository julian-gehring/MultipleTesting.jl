struct BocaLeek{T<:RegressionModel} <: Pi0RegressionEstimator
    λ::AbstractFloat
    fit_pars::Tuple
    function BocaLeek{T}(λ, fit_pars) where {T<:RegressionModel}
        isin(λ, 0., 1.) ? new{T}(λ, fit_pars) : throw(DomainError())
    end
end

BocaLeek{T<:RegressionModel}(::Type{T}, λ=0.5, fit_pars=()) = BocaLeek{T}(λ, fit_pars)

struct BocaLeekFit{T<:RegressionModel} <: Pi0RegressionFit
    λ::AbstractFloat
    model::T
end

function fit{T<:RegressionModel}(pi0est::BocaLeek{T}, pvalues, X)
    λ = pi0est.λ
    m = length(pvalues)
    Y = 1.0*(pvalues .> λ)
    model = fit(T, X, Y, pi0est.fit_pars...)
    BocaLeekFit{T}(λ, model)
end

# special-casing for Isotonic regression to avoid jump at boundary
function fit(pi0est::BocaLeek{IsotonicRegression}, pvalues, X)
    λ = pi0est.λ
    m = length(pvalues)
    min_idx = indmin(X)
    max_idx = indmax(X)
    Y = 1.0*(pvalues .> λ)
    Y[min_idx] = 1.0
    Y[max_idx] = 1.0
    model = fit(IsotonicRegression, X, Y, pi0est.fit_pars...)
    BocaLeekFit{IsotonicRegression}(λ, model)
end

function predict(pi0fit::BocaLeekFit, x...)
    λ = pi0fit.λ
    pi0_tmp = predict(pi0fit.model, x...)./(1-λ)
    max.(0.,min.(pi0_tmp, 1.))
end
