# The response-function derivative genealogy: the chain from any response
# back to the free energy, and the derivative order along it.
#
# Independent expectation: the textbook derivative structure —
# M = −∂F/∂h, χ = ∂M/∂h = ∂²F/∂h², U ⟵ F, C = ∂U/∂T — must be exactly
# what the genealogy reports, and every branch must root at FreeEnergy.

using AbstractQAtlas
using AbstractQAtlas:
    derivative_edge,
    differentiation_chain,
    potential_root,
    derivative_order,
    is_response,
    conjugate_field

@testset "genealogy edges" begin
    @test derivative_edge(Magnetization(:z)) == DerivativeEdge(FreeEnergy, MagneticField)
    @test derivative_edge(Susceptibility(:z, :z)) ==
        DerivativeEdge(Magnetization{:z}, MagneticField)
    @test derivative_edge(SpecificHeat()) == DerivativeEdge(Energy, Temperature)
    @test derivative_edge(ThermalEntropy()) == DerivativeEdge(FreeEnergy, Temperature)
    # Energy is dispatched through its {G} parameter
    @test derivative_edge(Energy()) == DerivativeEdge(FreeEnergy, InverseTemperature)
    # roots / non-genealogy quantities have no edge
    @test derivative_edge(FreeEnergy()) === nothing
    @test derivative_edge(PartitionFunction()) === nothing
    @test derivative_edge(CorrelationLength()) === nothing
end

@testset "is_response" begin
    @test is_response(Susceptibility(:z, :z))
    @test is_response(Magnetization(:x))
    @test is_response(SpecificHeat())
    @test !is_response(FreeEnergy())
    @test !is_response(PartitionFunction())
end

@testset "differentiation chain roots at the free energy" begin
    @test differentiation_chain(Susceptibility(:z, :z)) ==
        [typeof(Susceptibility(:z, :z)), typeof(Magnetization(:z)), FreeEnergy]
    @test differentiation_chain(Magnetization(:x)) ==
        [typeof(Magnetization(:x)), FreeEnergy]
    @test differentiation_chain(SpecificHeat()) == [SpecificHeat, Energy, FreeEnergy]
    @test differentiation_chain(FreeEnergy()) == [FreeEnergy]     # singleton root
    # every thermodynamic response roots at FreeEnergy
    for q in (
        Susceptibility(:x, :x),
        Magnetization(:z),
        SpecificHeat(),
        ThermalEntropy(),
        Energy(),
    )
        @test potential_root(q) === FreeEnergy
    end
end

@testset "derivative order along the chain" begin
    # χ is a SECOND field-derivative of F; M is a first; C none
    @test derivative_order(Susceptibility(:z, :z), MagneticField()) == 2
    @test derivative_order(Magnetization(:z), MagneticField()) == 1
    @test derivative_order(SpecificHeat(), MagneticField()) == 0
    # C is a first T-derivative of U (one T-edge on the U→…→F chain)
    @test derivative_order(SpecificHeat(), Temperature()) == 1
    # F itself: zero of anything
    @test derivative_order(FreeEnergy(), MagneticField()) == 0
end

@testset "nonlinear response: the genealogy extends recursively" begin
    # χ⁽²⁾_{x;yz} ⟵ χ⁽¹⁾_{x;y} ⟵ M_x ⟵ F  — each edge a field derivative
    χ2 = Susceptibility(:x, :y, :z)
    @test derivative_edge(χ2) ==
        DerivativeEdge(typeof(Susceptibility(:x, :y)), MagneticField)
    @test differentiation_chain(χ2) == [
        typeof(Susceptibility(:x, :y, :z)),
        typeof(Susceptibility(:x, :y)),
        typeof(Magnetization(:x)),
        FreeEnergy,
    ]
    @test potential_root(χ2) === FreeEnergy
    # χ⁽ⁿ⁾ is an (n+1)-th field-derivative of F
    @test derivative_order(Susceptibility(:x, :y), MagneticField()) == 2        # n=1
    @test derivative_order(Susceptibility(:x, :y, :z), MagneticField()) == 3    # n=2
    @test derivative_order(Susceptibility(:x, :x, :x, :x), MagneticField()) == 4  # n=3
    @test is_response(χ2)
end

@testset "conjugate fields" begin
    @test conjugate_field(Magnetization(:z)) == MagneticField()
    @test conjugate_field(Magnetization(:x)) == MagneticField()
    @test conjugate_field(ThermalEntropy()) == Temperature()
end

@testset "exact-formula companions (the relations the edges point at)" begin
    # χ = ∂M/∂h  (SusceptibilityResponse) and χ = β·Var(M) (SusceptibilityFDT)
    # are the same response two ways — both live in the registry now.
    @test check(SusceptibilityResponse(); χ=2.5, dM_dh=2.5)
    @test check(MagnetizationResponse(); M=1.0, dF_dh=-1.0)   # M = −∂F/∂h
    # the definitional and statistical routes agree on a shared value:
    dM_dh = 0.8               # ∂M/∂h at some point
    var_M = 1.6
    β = 0.5                   # β·Var(M) = 0.8 = ∂M/∂h  ⇒ same χ
    χ_def = solve(SusceptibilityResponse(), Val(:χ); dM_dh=dM_dh)
    χ_stat = solve(SusceptibilityFDT(), Val(:χ); var_M=var_M, β=β, N=1)
    @test χ_def ≈ χ_stat
end

@testset "grand-canonical branch: N roots at the GrandPotential (second root)" begin
    # N = −∂Ω/∂μ is the grand-canonical analogue of M = −∂F/∂h — a first
    # single-field derivative, but of the GRAND potential, not the free energy.
    @test derivative_edge(ParticleNumber()) ==
        DerivativeEdge(GrandPotential, ChemicalPotential)
    @test potential_root(ParticleNumber()) === GrandPotential
    @test differentiation_chain(ParticleNumber()) == [ParticleNumber, GrandPotential]
    @test is_response(ParticleNumber())
    @test derivative_order(ParticleNumber(), ChemicalPotential()) == 1
    @test conjugate_field(ParticleNumber()) === ChemicalPotential()
    # GrandPotential is itself a root (no edge); the canonical branch is untouched
    @test derivative_edge(GrandPotential()) === nothing
    @test differentiation_chain(GrandPotential()) == [GrandPotential]
    @test potential_root(Magnetization(:z)) === FreeEnergy
end
