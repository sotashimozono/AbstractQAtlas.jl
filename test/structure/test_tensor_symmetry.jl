# Intrinsic permutation symmetry of the nonlinear response tensors —
# the structural fact that makes the higher-order susceptibility a proper
# tensor: its field indices (with their frequencies) are interchangeable.

using AbstractQAtlas
using AbstractQAtlas:
    intrinsic_permutation_symmetric, canonical_component, permutation_equivalent

@testset "which quantities carry intrinsic permutation symmetry" begin
    @test intrinsic_permutation_symmetric(Susceptibility(:x, :y, :z))
    @test intrinsic_permutation_symmetric(DynamicalSusceptibility(:x, :y, :z))
    @test intrinsic_permutation_symmetric(Conductivity(:x, :y, :z))
    @test !intrinsic_permutation_symmetric(Energy())
    @test !intrinsic_permutation_symmetric(SpontaneousMagnetization())
    @test !intrinsic_permutation_symmetric(Magnetization(:x))   # rank-1, nothing to permute
end

@testset "canonical component sorts the FIELD indices (response index fixed)" begin
    # the response index (first) is fixed; field indices sort
    @test canonical_component(Susceptibility(:x, :z, :y)) === Susceptibility(:x, :y, :z)
    @test canonical_component(Susceptibility(:z, :z, :x)) === Susceptibility(:z, :x, :z)
    # linear (one field index): already canonical
    @test canonical_component(Susceptibility(:x, :y)) === Susceptibility(:x, :y)
    # dynamical and conductivity families too
    @test canonical_component(DynamicalSusceptibility(:x, :z, :y)) ===
        DynamicalSusceptibility(:x, :y, :z)
    @test canonical_component(Conductivity(:z, :y, :x)) === Conductivity(:z, :x, :y)
end

@testset "permutation equivalence" begin
    # χ⁽²⁾_{x;yz} = χ⁽²⁾_{x;zy}: same response index, permuted fields
    @test permutation_equivalent(Susceptibility(:x, :y, :z), Susceptibility(:x, :z, :y))
    # different RESPONSE index ⇒ not equivalent
    @test !permutation_equivalent(Susceptibility(:x, :y, :z), Susceptibility(:y, :x, :z))
    # different field content ⇒ not equivalent
    @test !permutation_equivalent(Susceptibility(:x, :y, :y), Susceptibility(:x, :y, :z))
    # reflexive on the linear component
    @test permutation_equivalent(Susceptibility(:x, :y), Susceptibility(:x, :y))
    # different families never equivalent
    @test !permutation_equivalent(Susceptibility(:x, :y, :z), Conductivity(:x, :y, :z))
    # non-symmetric quantities: not equivalent (no permutation freedom)
    @test !permutation_equivalent(Energy(), Energy())
end
