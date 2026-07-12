# Spin-glass / disordered-magnet identities: the exact Nishimori-line
# results and the de Almeida–Thouless replica-symmetry-breaking boundary.

using AbstractQAtlas
using AbstractQAtlas:
    check, solve, slack, residual, domain, variables, tensor_rank, AbstractInequality

@testset "Edwards–Anderson order parameter" begin
    @test check(EdwardsAndersonOrderParameter(); q_EA=0.42, overlap=0.42, atol=1e-12)
    @test solve(EdwardsAndersonOrderParameter(), Val(:q_EA); overlap=0.7) ≈ 0.7
    @test domain(EdwardsAndersonOrderParameter()) == :spinglass
end

@testset "Nishimori-line exact identities" begin
    # U = −J tanh(βJ): the exact ±J energy per bond on the Nishimori line
    J, β = 1.3, 0.8
    U = -J * tanh(β * J)
    @test check(NishimoriEnergy(); U=U, J=J, β=β, atol=1e-12)
    @test solve(NishimoriEnergy(), Val(:U); J=J, T=1 / β) ≈ U          # β↔T convention
    # physical limits: T→∞ ⇒ U→0 (disordered), T→0 ⇒ U→−J (bonds satisfied)
    @test solve(NishimoriEnergy(), Val(:U); J=J, β=1e-6) ≈ 0.0 atol = 1e-5
    @test solve(NishimoriEnergy(), Val(:U); J=J, β=1e3) ≈ -J atol = 1e-9
    # gauge identity q = m: spin-glass order equals magnetization on the NL
    @test check(NishimoriMagnetizationOverlap(); q=0.55, m=0.55, atol=1e-12)
    @test !check(NishimoriMagnetizationOverlap(); q=0.55, m=0.30, atol=1e-9)
end

@testset "de Almeida–Thouless replica-symmetry stability" begin
    @test AlmeidaThoulessStability() isa AbstractInequality
    # SK model in zero field (h=0 ⇒ sech⁴(0)=1): the AT line is βJ = 1 (T_c = J).
    # replicon eigenvalue = 1 − (βJ)²
    @test slack(AlmeidaThoulessStability(); βJ=1.0, sech4_avg=1.0) == 0.0     # ON the AT line / SK transition
    @test check(AlmeidaThoulessStability(); βJ=0.8, sech4_avg=1.0)           # T > J: replica-symmetric (stable)
    @test !check(AlmeidaThoulessStability(); βJ=1.5, sech4_avg=1.0, atol=1e-9)  # T < J: RSB
    # in a field the local-field factor softens the transition (sech⁴ < 1)
    @test check(AlmeidaThoulessStability(); βJ=1.2, sech4_avg=0.5)
    @test slack(AlmeidaThoulessStability(); βJ=1.2, sech4_avg=0.5) ≈ 1 - 1.44 * 0.5 atol =
        1e-12
end

@testset "spin-glass quantities + wiring" begin
    @test tensor_rank(EdwardsAndersonParameter()) == 0
    @test tensor_rank(SpinGlassSusceptibility()) == 0
    @test variables(NishimoriEnergy()) == (:U, :J, :β)
    @test domain(AlmeidaThoulessStability()) == :spinglass
end
