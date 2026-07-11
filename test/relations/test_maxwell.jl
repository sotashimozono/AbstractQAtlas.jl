# Maxwell relations, Clausius–Clapeyron, and Gibbs–Duhem.
#
# The two most-used Maxwell relations are checked on the ideal gas by
# finite-differencing TWO independent equations of state (S and p, or S
# and V) — a genuine cross-check, not an algebraic tautology.  The
# remaining relations are checked on constructed consistent values.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve

# ideal gas, Nk = 1:  p = T/V,  S = ln V + (3/2) ln T,  V = T/p
p_TV(T, V) = T / V
S_TV(T, V) = log(V) + 1.5 * log(T)
S_Tp(T, p) = log(T / p) + 1.5 * log(T)
V_Tp(T, p) = T / p

@testset "Maxwell (Helmholtz): (∂S/∂V)_T = (∂p/∂T)_V — ideal gas" begin
    h = 1e-6
    for (T, V) in ((2.0, 3.0), (1.3, 5.0), (4.0, 1.7))
        dS_dV = (S_TV(T, V + h) - S_TV(T, V - h)) / (2h)   # = 1/V
        dp_dT = (p_TV(T + h, V) - p_TV(T - h, V)) / (2h)   # = 1/V
        @test check(MaxwellHelmholtz(); dS_dV=dS_dV, dp_dT=dp_dT, atol=1e-6)
    end
end

@testset "Maxwell (Gibbs): (∂S/∂p)_T = −(∂V/∂T)_p — ideal gas" begin
    h = 1e-6
    for (T, p) in ((2.0, 0.8), (1.5, 1.2), (3.3, 0.5))
        dS_dp = (S_Tp(T, p + h) - S_Tp(T, p - h)) / (2h)   # = −1/p
        dV_dT = (V_Tp(T + h, p) - V_Tp(T - h, p)) / (2h)   # = 1/p
        @test check(MaxwellGibbs(); dS_dp=dS_dp, dV_dT=dV_dT, atol=1e-6)
    end
end

@testset "Maxwell (Internal / Enthalpy): sign structure" begin
    # (∂T/∂V)_S = −(∂p/∂S)_V  and  (∂T/∂p)_S = (∂V/∂S)_p
    @test check(MaxwellInternal(); dT_dV=-0.4, dp_dS=0.4)
    @test !check(MaxwellInternal(); dT_dV=0.4, dp_dS=0.4)
    @test check(MaxwellEnthalpy(); dT_dp=0.9, dV_dS=0.9)
    @test !check(MaxwellEnthalpy(); dT_dp=0.9, dV_dS=-0.9)
    @test solve(MaxwellInternal(), Val(:dT_dV); dp_dS=0.7) == -0.7
end

@testset "Clausius–Clapeyron: dp/dT = L/(T ΔV)" begin
    L, T, ΔV = 1.5, 2.0, 0.5
    @test check(ClausiusClapeyron(); dp_dT=L / (T * ΔV), L=L, T=T, ΔV=ΔV, atol=1e-13)
    @test solve(ClausiusClapeyron(), Val(:L); dp_dT=L / (T * ΔV), T=T, ΔV=ΔV) ≈ L
    @test !check(ClausiusClapeyron(); dp_dT=0.0, L=L, T=T, ΔV=ΔV)
end

@testset "Gibbs–Duhem: S dT − V dp + N dμ = 0" begin
    S, V, N, dT, dp = 2.0, 3.0, 1.5, 0.1, 0.05
    dμ = (V * dp - S * dT) / N                     # constructed to satisfy the constraint
    @test check(GibbsDuhem(); S=S, dT=dT, V=V, dp=dp, N=N, dμ=dμ, atol=1e-14)
    @test solve(GibbsDuhem(), Val(:dμ); S=S, dT=dT, V=V, dp=dp, N=N) ≈ dμ
    # exact arithmetic when the inputs are rational
    r = residual(
        GibbsDuhem(); S=2 // 1, dT=1 // 1, V=1 // 1, dp=2 // 1, N=1 // 1, dμ=0 // 1
    )
    @test r isa Rational && r == 0 // 1
end
