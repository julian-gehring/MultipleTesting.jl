__precompile__()

"""
*MultipleTesting* package

* Adjusting p-values for multiple testing
* Estimating the fraction π0 of tests under the null hypothesis
"""
module MultipleTesting

using Compat

using StatsBase
import StatsBase: fit, predict

using Distributions

export
    PValues,
    adjust,
    PValueAdjustment,
    Bonferroni,
    BenjaminiHochberg,
    BenjaminiHochbergAdaptive,
    BenjaminiYekutieli,
    BenjaminiLiu,
    Hochberg,
    Holm,
    Hommel,
    Sidak,
    ForwardStop,
    BarberCandes,
    qValues,
    estimate_pi0,
    Pi0Estimator,
    Storey,
    StoreyBootstrap,
    LeastSlope,
    Oracle,
    TwoStep,
    RightBoundary,
    CensoredBUM,
    CensoredBUMFit,
    BUM,
    BUMFit,
    FlatGrenander,
    isin,
    fit,
    BetaUniformMixtureModel,
    PValueCombination,
    combine,
    FisherCombination,
    LogitCombination,
    StoufferCombination,
    TippettCombination,
    SimesCombination,
    WilkinsonCombination,
    MinimumCombination,
    PriorityWeights,
    WeightedPValues,
    IsotonicRegression,
    Pi0RegressionEstimator,
    Pi0RegressionFit,
    BocaLeek,
    BocaLeekFit


include("types.jl")
include("shape-constrained.jl")
include("weights.jl")
include("utils.jl")
include("pval-adjustment.jl")
include("pi0-estimators.jl")
include("pi0-regression.jl")
include("qvalue.jl")
include("model.jl")
include("combinations.jl")

end
