### estimators for π0 (pi0) ###

## Storey estimator ##

"""
Storey π0 estimator

**Parameters**

- λ : tuning parameter, FloatingPoint, default: 0.1

**Examples**

```julia
Storey()
Storey(0.1)
```

**References**

Storey, JD (2002). "A Direct Approach to False Discovery Rates." Journal of the
Royal Statistical Society, doi:10.1111/1467-9868.00346

"""
type Storey <: Pi0Estimator
    λ::AbstractFloat

    Storey(λ) = isin(λ, 0., 1.) ? new(λ) : throw(DomainError())
end

Storey() = Storey(0.1)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::Storey)
    storey_pi0(pValues, pi0estimator.λ)
end

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, ws::MultipleTestingWeights{T}, pi0estimator::Storey)
    storey_pi0(pValues./ws, pi0estimator.λ)
end

function storey_pi0{T<:AbstractFloat}(pValues::Vector{T}, lambda::T)
    validPValues(pValues)
    pi0 = (sum(pValues .>= lambda) / length(pValues)) / (1.-lambda)
    pi0 = min(pi0, 1.)
    return pi0
end


## Storey bootstrap estimator ##
"""
Storey closed-form bootstrap π0 estimator

StoreyBootstrap(λseq, q)

Reference: David Robinson, 2012
"""
type StoreyBootstrap <: Pi0Estimator
    ## check range of arguments
    λseq::Vector{AbstractFloat}
    q   ::AbstractFloat

    StoreyBootstrap(λseq, q) =
        isin(λseq, 0., 1.) && isin(q, 0., 1.) ? new(λseq, q) : throw(DomainError())
end

StoreyBootstrap() = StoreyBootstrap(collect(0.05:0.05:0.95), 0.1)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::StoreyBootstrap)
    bootstrap_pi0(pValues, pi0estimator.λseq, pi0estimator.q)
end

function bootstrap_pi0{T<:AbstractFloat,S<:AbstractFloat}(pValues::Vector{T}, lambda::Vector{S} = [0.05:0.05:0.95;], q::S = 0.1)
    #validPValues(pValues)
    n = length(pValues)
    w = Int[sum(pValues .>= l) for l in lambda] ## TODO: check >= or >
    pi0 = w ./ n ./ (1. - lambda)
    min_pi0 = quantile(pi0, q)
    mse = (w ./ (n.^2 .* (1. - lambda).^2 )) .* (1. - w/n) + (pi0 - min_pi0).^2
    pi0 = min(pi0[indmin(mse)], 1.)
    pi0
end


## Least SLope (LSL) estimator ##

"""
Least SLope (LSL) π0 estimator

LeastSlope()
"""
type LeastSlope <: Pi0Estimator
end

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::LeastSlope)
    lsl_pi0(pValues)
end

function lsl_pi0{T<:AbstractFloat}(pValues::Vector{T})
    n = length(pValues)
    if !issorted(pValues)
        pValues = sort(pValues)
    end
    s0 = lsl_slope(1, n, pValues)
    sx = 0.
    for i in 2:n
        s1 = lsl_slope(i, n, pValues)
        if (s1 - s0) < 0.
            sx = s1
            break
        end
        s0 = s1
    end
    pi0 = min( 1/sx + 1, n ) / n
    return(pi0)
end

function lsl_slope{T<:AbstractFloat}(i::Int, n::Int, pval::Vector{T})
    s = (1 - pval[i]) / (n - i + 1)
    return s
end

## alternative, vectorized version
## used for comparison and compactness
function lsl_pi0_vec{T<:AbstractFloat}(pValues::Vector{T})
    n = length(pValues)
    if !issorted(pValues)
        pValues = sort(pValues)
    end
    s = (1 - pValues) ./ (n - collect(1:n) + 1)
    d = diff(s) .< 0
    idx = findfirst(d) + 1
    pi0 = min( 1/s[idx] + 1, n ) / n
    return(pi0)
end


## Oracle
"""
Oracle π0

Oracle(π0)
"""
type Oracle <: Pi0Estimator
    π0::AbstractFloat

    Oracle(π0) = isin(π0, 0., 1.) ? new(π0) : throw(DomainError())
end

Oracle() = Oracle(1.0)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::Oracle)
    pi0estimator.π0
end


## Two-Step estimator: Benjamini, Krieger and Yekutieli (2006) ##

"""
Two step π0 estimator

TwoStep(α)

Reference: Benjamini, Krieger and Yekutieli, 2006
"""
type TwoStep <: Pi0Estimator
    α::AbstractFloat
    ## method::PValueAdjustmentMethod

    TwoStep(α) = isin(α, 0., 1.) ? new(α) : throw(DomainError())
end

TwoStep() = TwoStep(0.05)

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::TwoStep)
    twostep_pi0(pValues, pi0estimator.α)
end

function twostep_pi0{T<:AbstractFloat}(pValues::Vector{T}, alpha::T, method::PValueAdjustmentMethod = BenjaminiHochberg())
    padj = adjust(pValues, method)
    pi0 = sum(padj .>= (alpha/(1+alpha))) / length(padj)
    return(pi0)
end


# RightBoundary procedure as defined in Definition 2 of Liang and Nettleton 2012
# "Adaptive and dynamic adaptive procedures for false discovery rate control and estimation"
"""
Right boundary π0 estimator

RightBoundary(λseq)
"""
type RightBoundary <: Pi0Estimator
    ## check range of arguments
    λseq::Vector{Float64}

    RightBoundary(λseq) =
        isin(λseq, 0., 1.) ? new(λseq) : throw(DomainError())
end

# λseq used in Liang, Nettleton 2012
RightBoundary() = RightBoundary([collect(0.02:0.02:0.1); collect(0.15:0.05:0.95)])

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::RightBoundary)
    rightboundary_pi0(pValues, pi0estimator.λseq)
end

function rightboundary_pi0{T<:AbstractFloat}(pValues::Vector{T}, λseq)
    validPValues(pValues)
    n = length(pValues)
    if !issorted(λseq)
        λseq = sort(λseq)
    end
    # make sure we catch p-values equal to 1 despite left closure
    # use closed=:left because we have been using >= convention in this package
    # note that original paper uses > convention.
    h = fit(Histogram, pValues, [λseq; Inf], closed=:left)
    pi0_estimates = reverse(cumsum(reverse(h.weights)))./(1.-λseq)./n
    pi0_decrease = diff(pi0_estimates) .>= 0
    pi0_decrease[end] = true
    pi0 = pi0_estimates[findfirst(pi0_decrease, true) + 1]
    return(min(pi0,1))
end


## Censored BUM
"""
Censored BUM π0 estimator

CensoredBUM(γ0, λ, xtol, maxiter)
"""
type CensoredBUM <: Pi0Estimator
    γ0::Float64
    λ::Float64
    xtol::Float64
    maxiter::Int64

    function CensoredBUM(γ0, λ, xtol, maxiter)
        if isin(γ0, 0., 1.) && isin(λ, 0., 1.) && isin(xtol, 0., 1.) && maxiter > 0
            new(γ0, λ, xtol, maxiter)
        else
            throw(DomainError())
        end
    end
end

CensoredBUM() = CensoredBUM(0.5, 0.05, 1e-6, 10000)
CensoredBUM(γ0, λ) = CensoredBUM(γ0, λ, 1e-6, 10000)

type CensoredBUMFit <: Pi0Fit
    π0::Float64
    param::Vector{Float64}
    is_converged::Bool
end

function fit(pi0estimator::CensoredBUM, pValues;
             kw...)
    π0, param, is_converged = cbum_pi0(pValues, pi0estimator.γ0, pi0estimator.λ,
                                       pi0estimator.xtol, pi0estimator.maxiter;
                                       kw...)
    return CensoredBUMFit(π0, param, is_converged)
end

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::CensoredBUM)
    estimate_pi0(fit(pi0estimator, pValues))
end

function estimate_pi0(pi0fit::CensoredBUMFit)
    π0 = pi0fit.is_converged ? π0 = pi0fit.π0 : NaN
    return π0
end

function cbum_pi0{T<:AbstractFloat}(pValues::Vector{T}, γ0::T = 0.5, λ::T = 0.05,
                                    xtol::Float64 = 1e-6, maxiter::Int = 10000)
    n = length(pValues)
    idx_right = pValues .>= λ
    n2 = sum(idx_right)
    n1 = n - n2
    sz = (1-γ0)*n
    szr = (1-γ0)*n2
    szl = sz - szr
    ## compute constant values only once
    pr = pValues[idx_right]
    lpr = log(pr)
    ll = log(λ)
    zr = fill(1-γ0, n2)
    pi0_old = γ0 = α = γ = Inf
    for i in 1:maxiter
        γ = 1 - sz/n
        α = -szr / ( ll * szl + sum(zr .* lpr) )
        γ = max(min(γ, 1.), 0.)
        α = max(min(α, 1.), 0.)
        xl = (1-γ) * (λ^α)
        szl = (xl ./ (γ*λ + xl)) * n1
        xr = (1-γ) * α * pr.^(α-1)
        zr = xr ./ (γ + xr)
        szr = sum(zr)
        sz = szl + szr
        pi0_new = γ + (1-γ)*α
        if abs(pi0_new - pi0_old) <= xtol
            return pi0_new, [γ, α], true
        end
        γ0 = γ
        pi0_old = pi0_new
    end
    return NaN, [γ, α], false
end

function cbum_pi0_naive(pValues, γ0 = 0.5, λ = 0.05,
                        xtol = 1e-6, maxiter = 10000)
    n = length(pValues)
    z = fill(1-γ0, n)
    idx_left = pValues .< λ
    idx_right = !idx_left
    pi0_old = γ0 = α = γ = Inf
    ## compute constant values only once
    lpr = log(pValues[idx_right])
    ll = log(λ)
    for i in 1:maxiter
        γ = sum(1-z) / n
        α = -sum(z[idx_right])
        α = α / ( ll * sum(z[idx_left]) + sum(z[idx_right] .* lpr) )
        xl = (1-γ) * (λ^α)
        z[idx_left] = xl ./ (γ*λ + xl)
        xr = (1-γ) * α * pValues[idx_right].^(α-1)
        z[idx_right] = xr ./ (γ + xr)
        pi0_new = γ + (1-γ)*α
        if abs(pi0_new - pi0_old) <= xtol
            return pi0_new, [γ, α], true
        end
        γ0 = γ
        pi0_old = pi0_new
    end
    return NaN, [γ, α], false
end


## BUM
"""
BUM π0 estimator

BUM(γ0, xtol, maxiter)
"""
type BUM <: Pi0Estimator
    γ0::Float64
    xtol::Float64
    maxiter::Int64

    function BUM(γ0, xtol, maxiter)
        if isin(γ0, 0., 1.) && isin(xtol, 0., 1.)
            new(γ0, xtol, maxiter)
        else
            throw(DomainError())
        end
    end
end

BUM() = BUM(0.5, 1e-6, 10000)

BUM(y0::Float64) = BUM(y0, 1e-6, 10000)

type BUMFit <: Pi0Fit
    π0::Float64
    param::Vector{Float64}
    is_converged::Bool
end

function fit(pi0estimator::BUM, pValues;
             kw...)
    π0, param, is_converged = cbum_pi0(pValues, pi0estimator.γ0, eps(),
                                       pi0estimator.xtol, pi0estimator.maxiter;
                                       kw...)
    return BUMFit(π0, param, is_converged)
end

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::BUM)
    estimate_pi0(fit(pi0estimator, pValues))
end

function estimate_pi0(pi0fit::BUMFit)
    π0 = pi0fit.is_converged ? π0 = pi0fit.π0 : NaN
    return π0
end


## Longest constant interval in the Grenander estimator: Langaas et al., 2005
# Leave unexported until thoroughly tested

"""
Flat Grenander π0 estimator

FlatGrenander()

Estimates π0 by the longest constant interval in the Grenander estimator

Reference: Langaas et al., 2005: section 4.3
"""
type FlatGrenander <: Pi0Estimator
end

function estimate_pi0{T<:AbstractFloat}(pValues::Vector{T}, pi0estimator::FlatGrenander)
    flat_grenander_pi0(pValues)
end

function flat_grenander_pi0(pv::Vector{Float64})
    p, f, F = grenander(pv)
    pi0 = longest_constant_interval(p, f)
    return pi0
end

function longest_constant_interval(p::Vector{Float64}, f::Vector{Float64})
    p = [0.0; p]
    f = [Inf; f]
    i2 = length(f)
    Δp_max = Δp = 0.0
    pi0 = 1.0
    for i1 in (length(f)-1):-1:1
        if f[i2] ≈ f[i1] ## within constant interval
            Δp = p[i2] - p[i1]
        else
            if Δp >= Δp_max
                Δp_max = Δp
                pi0 = f[i2]
            end
            i2 = i1
            Δp = 0.0
        end
        if f[i1] > 1.0
            break
        end
    end
    return pi0
end
