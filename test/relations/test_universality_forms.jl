# FSS forms: synthetic exponent round trips.
#
# Data generated FROM a form with a known exponent must return exactly
# that exponent under the standard log-log extraction — an arithmetic
# round trip that pins the exponent-combination conventions (γ/ν
# vs −γ/ν etc.) the docstrings promise.

using AbstractQAtlas
using Test

# least-squares slope of y vs x (hand-rolled: no Statistics dependency)
function _slope(x, y)
    n = length(x)
    mx = sum(x) / n
    my = sum(y) / n
    return sum((x .- mx) .* (y .- my)) / sum((x .- mx) .^ 2)
end

@testset "fss_peak_scaling round trip (2D Ising γ/ν = 7/4)" begin
    ratio = 7 / 4
    Ls = [8, 16, 32, 64]
    χmax = [2.31 * fss_peak_scaling(L; ratio=ratio) for L in Ls]
    slope = _slope(log.(Ls), log.(χmax))
    @test isapprox(slope, ratio; atol=1e-12)
end

@testset "order-parameter / correlation-length / susceptibility forms" begin
    ts = [-0.2, -0.1, -0.05, -0.025]
    β = 1 / 8
    Ms = [order_parameter_form(t; β=β) for t in ts]
    @test isapprox(_slope(log.(abs.(ts)), log.(Ms)), β; atol=1e-12)

    ν = 1.0
    ξs = [correlation_length_form(t; ν=ν) for t in ts]
    @test isapprox(_slope(log.(abs.(ts)), log.(ξs)), -ν; atol=1e-12)

    γ = 7 / 4
    χs = [susceptibility_form(t; γ=γ) for t in ts]
    @test isapprox(_slope(log.(abs.(ts)), log.(χs)), -γ; atol=1e-12)
end

@testset "collapse_coordinates" begin
    Tc = 2.269185314213022
    # the collapse pivot: at T = Tc the scaling variable is exactly 0, ∀L
    for L in (8, 16, 32)
        c = collapse_coordinates(Tc, L, Tc; ν=1.0, ratio=1 / 8)
        @test c.x == 0.0
        @test c.scale == L^(1 / 8)
    end
    # scaling-variable convention: x = (T − Tc)·L^{1/ν}
    c = collapse_coordinates(Tc + 0.01, 16, Tc; ν=1.0, ratio=-7 / 4)
    @test c.x ≈ 0.01 * 16.0
    @test c.scale ≈ 16.0^(-7 / 4)
    # perfect-collapse round trip: synthetic O(T, L) = L^{-ratio} f(x)
    f(x) = exp(-x^2)                    # any smooth universal curve
    ratio = 1 / 8
    ν = 1.0
    for L1 in (8, 64), L2 in (16, 32)
        x = 0.37                        # same scaling variable for both sizes
        T1 = Tc + x * L1^(-1 / ν)
        T2 = Tc + x * L2^(-1 / ν)
        O1 = L1^(-ratio) * f(collapse_coordinates(T1, L1, Tc; ν=ν, ratio=ratio).x)
        O2 = L2^(-ratio) * f(collapse_coordinates(T2, L2, Tc; ν=ν, ratio=ratio).x)
        s1 = collapse_coordinates(T1, L1, Tc; ν=ν, ratio=ratio).scale
        s2 = collapse_coordinates(T2, L2, Tc; ν=ν, ratio=ratio).scale
        @test isapprox(O1 * s1, O2 * s2; rtol=1e-12)   # sizes collapse
    end
end
