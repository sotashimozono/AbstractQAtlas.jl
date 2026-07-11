# Wick's theorem vs INDEPENDENT many-body Fock-space ED.
#
# A 3-site free-fermion chain is solved exactly in the full 8-dimensional
# Fock space (with explicit Jordan–Wigner sign bookkeeping).  Both the
# 2-point matrix G AND the 4-point functions are computed from the
# many-body thermal state directly; Wick's determinant factorization of
# the 4-points from G is then a genuinely independent identity check
# (ED trace vs determinant formula).

using AbstractQAtlas
using LinearAlgebra
using Test

# ── Fock-space machinery (3 sites, 8 states, JW convention) ─────────────
# |n⟩ = (c†₁)^{n₁}(c†₂)^{n₂}(c†₃)^{n₃}|0⟩, bit k of the integer label is
# the occupation of site k.  c_i |n⟩ carries sign (−1)^{Σ_{k<i} n_k}.
const NS = 3
const DIM = 2^NS

function annihilator(i::Int)
    c = zeros(Float64, DIM, DIM)
    for n in 0:(DIM - 1)
        if (n >> (i - 1)) & 1 == 1
            sign = iseven(count_ones(n & ((1 << (i - 1)) - 1))) ? 1.0 : -1.0
            c[n - (1 << (i - 1)) + 1, n + 1] = sign
        end
    end
    return c
end
const CS = [annihilator(i) for i in 1:NS]
const CDS = [Matrix(c') for c in CS]

@testset "Wick vs Fock-space ED (3-site tight binding, thermal)" begin
    t = 1.0
    H = zeros(DIM, DIM)
    for i in 1:(NS - 1)
        H .+= -t .* (CDS[i] * CS[i + 1] .+ CDS[i + 1] * CS[i])
    end

    for β in (0.7, 2.5)
        ρ = exp(-β .* H)
        ρ ./= tr(ρ)

        # 2-point matrix from the many-body state
        G = [tr(ρ * CDS[i] * CS[j]) for i in 1:NS, j in 1:NS]

        @testset "β = $β" begin
            # 2-point round trip through wick_contraction (n = 1 determinant)
            for i in 1:NS, j in 1:NS
                @test wick_contraction(G, [i], [j]) ≈ G[i, j] atol = 1e-12
            end

            # densities and density-density correlations: ED vs Wick
            for i in 1:NS, j in 1:NS
                nn_ed = tr(ρ * (CDS[i] * CS[i]) * (CDS[j] * CS[j]))
                @test wick_density_correlation(G, i, j) ≈ nn_ed atol = 1e-12
            end

            # a genuinely off-diagonal 4-point: ⟨c†₁ c†₂ c₂ c₃⟩
            #   operator string  c†_{cr₁} c†_{cr₂} c_{an₂} c_{an₁}
            #   with cr = [1,2], an = [3,2]  (an₂ = 2, an₁ = 3)
            val_ed = tr(ρ * CDS[1] * CDS[2] * CS[2] * CS[3])
            @test wick_contraction(G, [1, 2], [3, 2]) ≈ val_ed atol = 1e-12

            # and a 6-point: ⟨c†₁c†₂c†₃ c₃c₂c₁⟩ = ⟨n₁n₂n₃⟩
            val6_ed = tr(ρ * (CDS[1] * CS[1]) * (CDS[2] * CS[2]) * (CDS[3] * CS[3]))
            @test wick_contraction(G, [1, 2, 3], [1, 2, 3]) ≈ val6_ed atol = 1e-12
        end
    end
end

@testset "input validation" begin
    G = [0.5 0.0; 0.0 0.5]
    @test_throws ErrorException wick_contraction(G, [1, 2], [1])
end
