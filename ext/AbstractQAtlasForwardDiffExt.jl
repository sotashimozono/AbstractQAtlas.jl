# ext/AbstractQAtlasForwardDiffExt.jl — AD realization of the response
# genealogy via ForwardDiff.
#
# Each method evaluates the derivative a supplied-derivative relation
# needs directly from the potential function, so a downstream calculation
# can go from "a free-energy function F(h)" to "the magnetization" (and
# on to `check(MagnetizationResponse(); …)`) without hand-coding the
# derivative.

module AbstractQAtlasForwardDiffExt

using AbstractQAtlas
import AbstractQAtlas: thermal_derivative        # extended below → must import
using AbstractQAtlas:
    response_order, Magnetization, Susceptibility, ThermalEntropy, SpecificHeat, Energy
using ForwardDiff: derivative

# n-th derivative of a scalar function by nested ForwardDiff (n small —
# response orders are 1–3).
_nth(f, x, n::Integer) = n == 0 ? f(x) : _nth(y -> derivative(f, y), x, n - 1)

# M_α = −∂F/∂h  (first field-derivative of the free energy)
thermal_derivative(::Magnetization, F, h) = -derivative(F, h)

# χ⁽ⁿ⁾ = −∂ⁿ⁺¹F/∂hⁿ⁺¹ : an order-n susceptibility is the (n+1)-th field-
# derivative of F (single-field / diagonal; response_order = n).
thermal_derivative(χ::Susceptibility, F, h) = -_nth(F, h, response_order(χ) + 1)

# S = −∂F/∂T
thermal_derivative(::ThermalEntropy, F, T) = -derivative(F, T)

# C = ∂U/∂T
thermal_derivative(::SpecificHeat, U, T) = derivative(U, T)

# U = ∂(βF)/∂β  (Gibbs–Helmholtz; pass the βF function of β)
thermal_derivative(::Energy, βF, β) = derivative(βF, β)

end # module AbstractQAtlasForwardDiffExt
