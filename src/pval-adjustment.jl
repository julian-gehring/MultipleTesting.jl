## p-value adjustment methods ##

# promotion from float vectors to PValues type

function adjust(pValues::Vector{T}, method::M) where {T<:AbstractFloat, M<:PValueAdjustment}
    adjust(PValues(pValues), method)
end

function adjust(pValues::Vector{T}, n::Integer, method::M) where {T<:AbstractFloat, M<:PValueAdjustment}
    adjust(PValues(pValues), n, method)
end

# Bonferroni

"""
Bonferroni adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, Bonferroni())
4-element Array{Float64,1}:
 0.004
 0.04
 0.12
 1.0

julia> adjust(pvals, 6, Bonferroni())
4-element Array{Float64,1}:
 0.006
 0.06
 0.18
 1.0
```


# References

Bonferroni, C.E. (1936). Teoria statistica delle classi e calcolo delle
probabilita (Libreria internazionale Seeber).

"""
struct Bonferroni <: PValueAdjustment
end

adjust(pValues::PValues, method::Bonferroni) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Bonferroni)
    k = length(pValues)
    check_number_tests(k, n)
    return min.(pValues * n, 1)
end


# Benjamini-Hochberg

"""
Benjamini-Hochberg adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, BenjaminiHochberg())
4-element Array{Float64,1}:
 0.004
 0.02
 0.04
 0.5

julia> adjust(pvals, 6, BenjaminiHochberg())
4-element Array{Float64,1}:
 0.006
 0.03
 0.06
 0.75
```


# References

Benjamini, Y., and Hochberg, Y. (1995). Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. Journal of the Royal
Statistical Society. Series B (Methodological) 57, 289–300.

"""
struct BenjaminiHochberg <: PValueAdjustment
end

adjust(pValues::PValues, method::BenjaminiHochberg) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiHochberg)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= n ./ (1:k)
    stepup!(pAdjusted)
    pAdjusted = min.(pAdjusted[originalOrder], 1)
    return pAdjusted
end


# Benjamini-Hochberg Adaptive

"""
Adaptive Benjamini-Hochberg adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, BenjaminiHochbergAdaptive(Oracle(0.5))) # known π₀ of 0.5
4-element Array{Float64,1}:
 0.002
 0.01
 0.02
 0.25

julia> adjust(pvals, BenjaminiHochbergAdaptive(StoreyBootstrap())) # π₀ estimator
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> adjust(pvals, 6, BenjaminiHochbergAdaptive(StoreyBootstrap()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
```


# References

Benjamini, Y., and Hochberg, Y. (1995). Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. Journal of the Royal
Statistical Society. Series B (Methodological) 57, 289–300.

"""
struct BenjaminiHochbergAdaptive <: PValueAdjustment
    pi0estimator::Pi0Estimator
end

## default to BenjaminiHochberg
BenjaminiHochbergAdaptive(π0::AbstractFloat) = BenjaminiHochbergAdaptive(Oracle(π0))

BenjaminiHochbergAdaptive() = BenjaminiHochbergAdaptive(1.0)

adjust(pValues::PValues, method::BenjaminiHochbergAdaptive) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiHochbergAdaptive)
    π0 = estimate_pi0(pValues, method.pi0estimator)
    pAdjusted = adjust(pValues, n, BenjaminiHochberg()) * π0
    return pAdjusted
end


# Benjamini-Yekutieli

"""
Benjamini-Yekutieli adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, BenjaminiYekutieli())
4-element Array{Float64,1}:
 0.00833333
 0.0416667
 0.0833333
 1.0

julia> adjust(pvals, 6, BenjaminiYekutieli())
4-element Array{Float64,1}:
 0.0147
 0.0735
 0.147
 1.0
```


# References

Benjamini, Y., and Yekutieli, D. (2001). The Control of the False Discovery Rate
in Multiple Testing under Dependency. The Annals of Statistics 29, 1165–1188.

"""
struct BenjaminiYekutieli <: PValueAdjustment
end

adjust(pValues::PValues, method::BenjaminiYekutieli) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiYekutieli)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= harmonic_number(n) .* n ./ (1:k)
    stepup!(pAdjusted)
    pAdjusted = min.(pAdjusted[originalOrder], 1)
    return pAdjusted
end


# Benjamini-Liu

"""
Benjamini-Liu adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, BenjaminiLiu())
4-element Array{Float64,1}:
 0.003994
 0.0222757
 0.02955
 0.125

julia> adjust(pvals, 6, BenjaminiLiu())
4-element Array{Float64,1}:
 0.00598502
 0.0408416
 0.0764715
 0.4375
```


# References

Benjamini, Y., and Liu, W. (1999). A step-down multiple hypotheses testing
procedure that controls the false discovery rate under independence. Journal of
Statistical Planning and Inference 82, 163–170.

"""
struct BenjaminiLiu <: PValueAdjustment
end

adjust(pValues::PValues, method::BenjaminiLiu) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::BenjaminiLiu)
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    # a bit more involved because cutoffs at significance α have the form:
    # P_(i) <= 1- [1 - min(1, m/(m-i+1)α)]^{1/(m-i+1)}
    s = n .- (1:k) .+ 1
    pAdjusted = (1 .- (1 .- pAdjusted) .^ s) .* s ./ n
    stepdown!(pAdjusted)
    pAdjusted = min.(pAdjusted[originalOrder], 1)
    return pAdjusted
end


# Hochberg

"""
Hochberg adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, Hochberg())
4-element Array{Float64,1}:
 0.004
 0.03
 0.06
 0.5

julia> adjust(pvals, 6, Hochberg())
4-element Array{Float64,1}:
 0.006
 0.05
 0.12
 1.0
```


# References

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of
significance. Biometrika 75, 800–802.

"""
struct Hochberg <: PValueAdjustment
end

adjust(pValues::PValues, method::Hochberg) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Hochberg)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= (n .- (1:k) .+ 1)
    stepup!(pAdjusted)
    pAdjusted = min.(pAdjusted[originalOrder], 1)
    return pAdjusted
end


# Holm

"""
Holm adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, Holm())
4-element Array{Float64,1}:
 0.004
 0.03
 0.06
 0.5

julia> adjust(pvals, 6, Holm())
4-element Array{Float64,1}:
 0.006
 0.05
 0.12
 1.0
```


# References

Holm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure.
Scandinavian Journal of Statistics 6, 65–70.

"""
struct Holm <: PValueAdjustment
end

adjust(pValues::PValues, method::Holm) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Holm)
    k = length(pValues)
    check_number_tests(k, n)
    if n <= 1
        return pValues
    end
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    pAdjusted .*= (n .- (1:k) .+ 1)
    stepdown!(pAdjusted)
    pAdjusted = min.(pAdjusted[originalOrder], 1)
    return pAdjusted
end


# Hommel

"""
Hommel adjustment

# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, Hommel())
4-element Array{Float64,1}:
 0.004
 0.03
 0.06
 0.5

julia> adjust(pvals, 6, Hommel())
4-element Array{Float64,1}:
 0.006
 0.05
 0.12
 1.0
```


# References

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a
modified Bonferroni test. Biometrika 75, 383–386.

"""
struct Hommel <: PValueAdjustment
end

adjust(pValues::PValues, method::Hommel) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Hommel)
    k = length(pValues)
    check_number_tests(k, n)
    if k <= 1
        return pValues
    end
    pValues = vcat(pValues, fill(1.0, n-k))  # TODO avoid sorting of ones
    sortedOrder, originalOrder = reorder(pValues)
    pAdjusted = pValues[sortedOrder]
    q = fill(minimum(n .* pAdjusted./(1:n)), n)
    pa = fill(q[1], n)
    for j in (n-1):-1:2
        ij = 1:(n-j+1)
        i2 = (n-j+2):n
        q1 = minimum(j .* pAdjusted[i2]./((2:j)))
        q[ij] = min.(j .* pAdjusted[ij], q1)
        q[i2] = q[n-j+1]
        pa = max.(pa, q)
    end
    pAdjusted = max.(pa, pAdjusted)[originalOrder[1:k]]
    return pAdjusted
end


# Sidak

"""
Šidák adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, Sidak())
4-element Array{Float64,1}:
 0.003994
 0.039404
 0.114707
 0.9375

julia> adjust(pvals, 6, Sidak())
4-element Array{Float64,1}:
 0.00598502
 0.0585199
 0.167028
 0.984375
```


# References

Šidák, Z. (1967). Rectangular Confidence Regions for the Means of Multivariate
Normal Distributions. Journal of the American Statistical Association 62,
626–633.

"""
struct Sidak <: PValueAdjustment
end

adjust(pValues::PValues, method::Sidak) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::Sidak)
    check_number_tests(length(pValues), n)
    pAdjusted = min.(1 .- (1 .- pValues).^n, 1)
    return pAdjusted
end


# Forward Stop

"""
Forward-Stop adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, ForwardStop())
4-element Array{Float64,1}:
 0.0010005
 0.00552542
 0.0138367
 0.183664

julia> adjust(pvals, 6, ForwardStop())
4-element Array{Float64,1}:
 0.0010005
 0.00552542
 0.0138367
 0.183664
```


# References

G’Sell, M.G., Wager, S., Chouldechova, A., and Tibshirani, R. (2016). Sequential
selection procedures and false discovery rate control. J. R. Stat. Soc. B 78,
423–444.

"""
struct ForwardStop <: PValueAdjustment
end

adjust(pValues::PValues, method::ForwardStop) = adjust(pValues, length(pValues), method)

function adjust(pValues::PValues, n::Integer, method::ForwardStop)
    k = length(pValues)
    check_number_tests(k, n)
    sortedOrder, originalOrder = reorder(pValues)
    logsums = -cumsum(log.(1 .- pValues[sortedOrder]))
    logsums ./= (1:k)
    stepup!(logsums)
    pAdjusted = max.(min.(logsums[originalOrder], 1), 0)
    return pAdjusted
end


# Barber-Candès

"""
Barber-Candès adjustment


# Examples

```jldoctest
julia> pvals = PValues([0.001, 0.01, 0.03, 0.5]);

julia> adjust(pvals, BarberCandes())
4-element Array{Float64,1}:
 0.333333
 0.333333
 0.333333
 1.0
```


# References

Barber, R.F., and Candès, E.J. (2015). Controlling the false discovery rate via
knockoffs. Ann. Statist. 43, 2055–2085.

Arias-Castro, E., and Chen, S. (2017). Distribution-free multiple testing.
Electron. J. Statist. 11, 1983–2001.

"""
struct BarberCandes <: PValueAdjustment
end

function adjust(pValues::PValues, method::BarberCandes)
    n = length(pValues)
    if n <= 1
        return fill(1.0, size(pValues)) # unlike other p-adjust methods
    end

    sorted_indexes, original_order = reorder(pValues)
    estimated_fdrs = pValues[sorted_indexes]

    Rt = 1 # current number of discoveries
    Vt = 1 # estimated false discoveries at t

    left_pv = estimated_fdrs[1]
    right_pv = estimated_fdrs[n]

    while left_pv < 0.5
        while (1 - right_pv <= left_pv) && (right_pv >= 0.5) && (Vt + Rt <= n)
            estimated_fdrs[n-Vt+1] = 1.0
            right_pv = estimated_fdrs[n-Vt]
            Vt += 1
        end
        estimated_fdrs[Rt] = Vt/Rt
        Rt += 1
        left_pv = estimated_fdrs[Rt]
    end

    while (right_pv >= 0.5) && (Vt + Rt <= n+1)
      estimated_fdrs[n-Vt+1] = 1.0
      right_pv = (Vt + Rt <= n) ? estimated_fdrs[n-Vt] : 0.0
      Vt += 1
    end

    stepup!(estimated_fdrs)
    pAdjusted = min.(estimated_fdrs[original_order], 1)
    return pAdjusted
end

# as test, inefficient implementation
function barber_candes_brute_force(pValues::AbstractVector{T}) where T<:AbstractFloat
    n = length(pValues)
    sorted_indexes, original_order = reorder(pValues)
    sorted_pValues = pValues[sorted_indexes]
    estimated_fdrs = fill(1.0, size(pValues))
    for (i,pv) in enumerate(sorted_pValues)
        if pv >= 0.5
            break
        else
            estimated_fdrs[i] = (sum((1 .- pValues) .<= pv) + 1)/i
        end
    end
    stepup!(estimated_fdrs)
    pAdjusted = min.(estimated_fdrs[original_order], 1)
    return pAdjusted
end


## internal ##

# step-up / step-down

function stepup!(values::AbstractVector{T}) where T<:AbstractFloat
    accumulate!(min, values, reverse!(values))
    reverse!(values)
    return values
end

function stepdown!(values::AbstractVector{T}) where T<:AbstractFloat
    accumulate!(max, values, values)
    return values
end


function check_number_tests(k::Integer, n::Integer)
    if k > n
        msg = "Number of total tests ($n) is smaller than number of p-values ($k)"
        throw(ArgumentError(msg))
    end
end


harmonic_number(n::Integer) =  digamma(n+1) + γ
