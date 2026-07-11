# The dynamical-quantity inter-relationship graph: the operation edges,
# the chain from any derived quantity back to its source, and the tie
# from each edge to its exact pointwise relation (where one exists).

using AbstractQAtlas
using AbstractQAtlas: spectral_origin, spectral_chain, origin_relation

@testset "graph edges (which operation makes which quantity)" begin
    @test spectral_origin(RetardedGreensFunction()) == SpectralOrigin(SelfEnergy, :dyson)
    @test spectral_origin(SpectralFunction()) ==
        SpectralOrigin(RetardedGreensFunction, :neg_im_over_pi)
    @test spectral_origin(DensityOfStates()) ==
        SpectralOrigin(SpectralFunction, :bz_average)
    @test spectral_origin(DynamicalStructureFactor()) ==
        SpectralOrigin(DynamicalCorrelation, :spacetime_fourier)
    @test spectral_origin(NMRSpinRelaxationRate()) ==
        SpectralOrigin(DynamicalSusceptibility, :low_frequency_limit)
    # Kubo edge: the dynamical susceptibility (any response order) comes from
    # the correlation function; same source for linear and multi-time nonlinear
    @test spectral_origin(DynamicalSusceptibility(:x, :y)) ==
        SpectralOrigin(DynamicalCorrelation, :kubo)
    @test spectral_origin(DynamicalSusceptibility(:x, :y, :z)) ==
        SpectralOrigin(DynamicalCorrelation, :kubo)   # χ⁽²⁾(ω₁,ω₂), same origin
    # sources / off-graph quantities have no edge
    @test spectral_origin(SelfEnergy()) === nothing
    @test spectral_origin(DynamicalCorrelation()) === nothing
    @test spectral_origin(FreeEnergy()) === nothing
end

@testset "chains trace back to the source quantity" begin
    # the density of states is built from the self-energy: ρ ⟵ A ⟵ G^R ⟵ Σ
    @test spectral_chain(DensityOfStates()) ==
        [DensityOfStates, SpectralFunction, RetardedGreensFunction, SelfEnergy]
    @test spectral_chain(SpectralFunction()) ==
        [SpectralFunction, RetardedGreensFunction, SelfEnergy]
    @test spectral_chain(DynamicalStructureFactor()) ==
        [DynamicalStructureFactor, DynamicalCorrelation]
    # NMR ⟵ χ'' ⟵ (Kubo) correlation: 1/T₁ built from the correlation function
    @test spectral_chain(NMRSpinRelaxationRate()) ==
        [NMRSpinRelaxationRate, DynamicalSusceptibility, DynamicalCorrelation]
    @test spectral_chain(SelfEnergy()) == [SelfEnergy]           # source singleton
end

@testset "edges tie to their exact pointwise relation (or nothing)" begin
    # single-(q,ω)-point operations point at the @relation that realizes them
    @test origin_relation(:dyson) isa Dyson
    @test origin_relation(:neg_im_over_pi) isa SpectralFromGreens
    # transform / BZ-sum / limit operations have no pointwise relation here
    # (their evaluation is the functional sibling's job — issue #14)
    @test origin_relation(:bz_average) === nothing
    @test origin_relation(:spacetime_fourier) === nothing
    @test origin_relation(:low_frequency_limit) === nothing
    @test origin_relation(:kubo) === nothing            # transform of a multi-time correlation
    # every graph edge whose via is pointwise resolves to a registered relation
    for q in (RetardedGreensFunction(), SpectralFunction())
        rel = origin_relation(spectral_origin(q).via)
        @test rel !== nothing
        @test rel in all_relations(; domain=:spectral)
    end
end
