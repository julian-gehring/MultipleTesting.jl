abstract MTPSimulator

immutable MTPSimulation{M <: MTPSimulator, T <: AbstractFloat, S}
  P::Vector{T}    # pvalues
  X::Array{S}    # covariate independent under the null
  H::Vector{Bool} # is hypothesis null or not?
  seed::Int # using MerseneTwister type would be cooler to do eventually
  simulator::M
end

immutable TwoSampleTTests <: MTPSimulator
  Δμ::Float64
  π0::Float64
  ngroup::Int
  m::Int
end

function Distributions.rand(sim::TwoSampleTTests, seed=1)
  srand(seed)
  m  = sim.m
  ngroup = sim.ngroup
  Δμ = sim.Δμ
  π1 = min(max( 1 - sim.π0, .0), 1.)
  Z = zeros(Float64, ngroup, 2)
  Z1 = sub(Z, :, 1) # Z scores in group 1
  Z2 = sub(Z, :, 2) # Z scores in group 2
  X  = zeros(Float64, m)
  P = zeros(Float64, m)
  H  = rand(Bernoulli(π1), m)
  normal = Normal(2)
  normal_shifted = Normal(Δμ+2)

  Ztmp = zero(Float64) # initialize outside of loop
  @inbounds for i=1:m
    rand!(normal, Z1)
    if H[i] == 0
      rand!(normal, Z2)
    else
      rand!(normal_shifted, Z2)
    end
    X[i] = var(Z)
    Ztmp = abs(sqrt(ngroup/2)*(mean(Z1)- mean(Z2))/(sqrt(0.5*var(Z1) + 0.5*var(Z2))))
    P[i] = 2*(1-cdf(TDist(2*ngroup-2), Ztmp))
  end
  MTPSimulation(P, X, map(Bool, H), seed, sim)
end
