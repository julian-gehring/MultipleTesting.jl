# MultipleTesting

The `MultipleTesting` package offers common algorithms for p-value adjustment and the estimation of the proportion π₀ of true null hypotheses.

![xkcd p-value guide](pvalues.png)


## Package Status

[![Package Status](http://pkg.julialang.org/badges/MultipleTesting_0.5.svg)](http://pkg.julialang.org/?pkg=MultipleTesting)
[![Linux/Mac Build Status](https://travis-ci.org/juliangehring/MultipleTesting.jl.svg?branch=master)](https://travis-ci.org/juliangehring/MultipleTesting.jl)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/1ld0ppptisirryt1/branch/master?svg=true)](https://ci.appveyor.com/project/juliangehring/multipletesting-jl/branch/master)
[![Coverage Status](https://codecov.io/gh/juliangehring/MultipleTesting.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliangehring/MultipleTesting.jl)


## Features

### Adjustment of p-values

* Bonferroni
* Benjamini-Hochberg
* Adaptive Benjamini-Hochberg with known π₀ or π₀ estimation method (see section below)
* Benjamini-Yekutieli
* Benjamini-Liu
* Hochberg
* Holm
* Hommel
* Sidak
* Forward Stop

```julia
adjust(pvals, Bonferroni())
adjust(pvals, BenjaminiHochberg())
adjust(pvals, BenjaminiHochbergAdaptive(0.9))
adjust(pvals, BenjaminiHochbergAdaptive(Storey()))
adjust(pvals, BenjaminiYekutieli())
adjust(pvals, BenjaminiLiu())
adjust(pvals, Hochberg())
adjust(pvals, Holm())
adjust(pvals, Hommel())
adjust(pvals, Sidak())
adjust(pvals, ForwardStop())
```


### Estimation of π₀

* Storey
* Storey's closed-form bootstrap
* Least Slope
* Two Step
* RightBoundary (Storey's estimate with dynamically chosen λ)
* Beta-Uniform Mixture (BUM)
* Censored BUM
* Flat Grenander
* Oracle for known π₀

```julia
estimate_pi0(pvals, Storey())
estimate_pi0(pvals, StoreyBootstrap())
estimate_pi0(pvals, LeastSlope())
estimate_pi0(pvals, TwoStep())
estimate_pi0(pvals, TwoStep(0.05))
estimate_pi0(pvals, TwoStep(0.05, BenjaminiHochbergAdaptive(0.9))
estimate_pi0(pvals, RightBoundary())
estimate_pi0(pvals, CensoredBUM())
estimate_pi0(pvals, BUM())
estimate_pi0(pvals, FlatGrenander())
estimate_pi0(pvals, Oracle(0.9))
```
