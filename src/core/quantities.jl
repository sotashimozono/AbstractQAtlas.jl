# core/quantities.jl — the model-independent quantity vocabulary.
#
# Ported from QAtlas `src/core/quantities.jl`: the abstract quantity
# hierarchy plus a CURATED set of stable concrete tags (the ones the
# relations layer and its first consumers reference).  The full tag set
# (~60 structs: NMR, Loschmidt/DQPT, entanglement family, structure
# factors, correlation families, velocities, CFT data, …) migrates here
# incrementally as QAtlas adopts this package — tracked upstream; do not
# bulk-copy it ahead of need.
#
# Three tags that were previously stranded inside a *model* file
# (`PartitionFunction`, `CriticalTemperature`, `SpontaneousMagnetization`
# in QAtlas' IsingSquare.jl) are given their generic home here — they are
# meaningful for any statistical-mechanics model.

# ─── Abstract quantity hierarchy ────────────────────────────────────────

"""
    AbstractThermalPotential <: AbstractQuantity

Scalar thermodynamic potentials and their densities (energy, free
energy, entropy, specific heat, …).
"""
abstract type AbstractThermalPotential <: AbstractQuantity end

"""
    AbstractMagnetization <: AbstractQuantity

Order-parameter expectation values ⟨M_α⟩ and their site-resolved
variants.
"""
abstract type AbstractMagnetization <: AbstractQuantity end

"""
    AbstractSusceptibility <: AbstractQuantity

Linear and nonlinear response susceptibilities `χ⁽ⁿ⁾_{α;β₁…βₙ}` — see
[`Susceptibility`](@ref) for the arbitrary-order tensor.
"""
abstract type AbstractSusceptibility <: AbstractQuantity end

"""
    AbstractTwoPointCorrelation <: AbstractQuantity

Two-point correlation functions ⟨O_i O_j⟩ (connected or full).
"""
abstract type AbstractTwoPointCorrelation <: AbstractQuantity end

"""
    AbstractStructureFactor <: AbstractQuantity

Momentum-space structure factors S(q).
"""
abstract type AbstractStructureFactor <: AbstractQuantity end

"""
    AbstractGap <: AbstractQuantity

Spectral gaps (mass, charge, spin channels).
"""
abstract type AbstractGap <: AbstractQuantity end

"""
    AbstractVelocity <: AbstractQuantity

Characteristic velocities (Fermi, Luttinger, Lieb-Robinson, …).
"""
abstract type AbstractVelocity <: AbstractQuantity end

"""
    AbstractEntanglementMeasure <: AbstractQuantity

Entanglement measures (von Neumann / Rényi entropies, negativity,
mutual information, …).
"""
abstract type AbstractEntanglementMeasure <: AbstractQuantity end

export AbstractThermalPotential, AbstractMagnetization, AbstractSusceptibility
export AbstractTwoPointCorrelation, AbstractStructureFactor, AbstractGap
export AbstractVelocity, AbstractEntanglementMeasure

# The tensor traits (`tensor_rank`, `index_spaces`, `indices`) — the
# honest successors of the old fused `component` label — live in
# `core/indices.jl`; concrete tensor quantities add their methods below.

# ─── Scalar thermodynamics ──────────────────────────────────────────────

"""
    Energy{G}() <: AbstractThermalPotential
    Energy()                 # G = :natural — model-and-BC-natural granularity
    Energy(:total)           # explicit ⟨H⟩
    Energy(:per_site)        # explicit ⟨H⟩ / N

Ground-state / thermal energy expectation.  The type parameter `G` makes
the granularity (total vs per-site) a dispatch axis instead of a hidden
docstring contract.

`Energy()` resolves to the model's native granularity via the
[`native_energy_granularity`](@ref) trait.  Use the explicit
constructors when the caller needs a specific granularity.
"""
struct Energy{G} <: AbstractThermalPotential
    function Energy{G}() where {G}
        G isa Symbol || error("Energy granularity must be a Symbol, got $(typeof(G))")
        G in (:natural, :total, :per_site) ||
            error("unknown Energy granularity :$G; expected :natural, :total, or :per_site")
        return new{G}()
    end
end
Energy() = Energy{:natural}()
Energy(g::Symbol) = Energy{g}()
export Energy

"""
    native_energy_granularity(model, bc) -> :total | :per_site

Trait declaring which granularity the given `model` returns natively for
[`Energy`](@ref) at boundary condition `bc`.  Every model that supports
`Energy` must add a method per supported BC.  A missing method is caught
at the call site as a `MethodError`, which is intentional: it forces new
models to declare the convention rather than silently inheriting an
unrelated default.
"""
function native_energy_granularity end
export native_energy_granularity

"""
    FreeEnergy() <: AbstractThermalPotential

Helmholtz free energy per site, `f = -β⁻¹ log Z / N`.
"""
struct FreeEnergy <: AbstractThermalPotential end
export FreeEnergy

"""
    SpecificHeat() <: AbstractThermalPotential

Specific heat per site, `c_v(β) = β² (⟨H²⟩ − ⟨H⟩²) / N`.

The defining identity is available as a first-class relation:
[`SpecificHeatFDT`](@ref).
"""
struct SpecificHeat <: AbstractThermalPotential end
export SpecificHeat

"""
    ThermalEntropy() <: AbstractThermalPotential

Thermal entropy per site, `s = β(ε − f)`.
"""
struct ThermalEntropy <: AbstractThermalPotential end
export ThermalEntropy

"""
    PartitionFunction() <: AbstractThermalPotential

The partition function `Z(β) = Σ exp(-βE)` itself (finite systems).

Generic home for a tag that previously lived inside a model file:
any statistical-mechanics model with a finite configuration space can
register it.
"""
struct PartitionFunction <: AbstractThermalPotential end
export PartitionFunction

# ─── Order parameters and responses (tensors in spin space) ─────────────
#
# These were the quantities blurred into scalars by baking a component
# into the type name (`MagnetizationX`, `SusceptibilityZZ`).  The honest
# form carries the spin-axis index/indices as type parameters, so a
# component is a *selection* of the tensor and off-diagonal components
# (χ_xy) are expressible — not only the diagonal.  `_axis` validates that
# a parameter is a bare spin-axis symbol.

function _axis(a)
    return if a isa Symbol
        a
    else
        error("spin-axis index must be a Symbol (e.g. :x), got $(repr(a))")
    end
end

"""
    Magnetization{A}() <: AbstractMagnetization
    Magnetization(a::Symbol)

Uniform magnetization component `⟨M_A⟩` per site — a rank-1 tensor in
[`SpinAxis`](@ref) space, `A ∈ {:x, :y, :z, …}`.  `Magnetization(:z)`
replaces the old `MagnetizationZ`.
"""
struct Magnetization{A} <: AbstractMagnetization
    Magnetization{A}() where {A} = (_axis(A); new{A}())
end
Magnetization(a::Symbol) = Magnetization{a}()
tensor_rank(::Type{<:Magnetization}) = 1
index_spaces(::Type{<:Magnetization}) = (SpinAxis(),)
indices(::Type{Magnetization{A}}) where {A} = (A,)
export Magnetization

"""
    SpontaneousMagnetization() <: AbstractMagnetization

The spontaneous (symmetry-broken) order-parameter magnitude `M(T)` in
the ordered phase; identically zero above `T_c`.  A scalar — the
magnitude, not a spin component — so `tensor_rank == 0`; its critical
exponent is β (see `critical_scaling`).
"""
struct SpontaneousMagnetization <: AbstractMagnetization end
export SpontaneousMagnetization

function _axistuple(I)
    return if (I isa Tuple && !isempty(I) && all(a -> a isa Symbol, I))
        I
    else
        error(
            "index parameter must be a non-empty Tuple of spin-axis Symbols, got $(repr(I))"
        )
    end
end

"""
    Susceptibility{I}() <: AbstractSusceptibility
    Susceptibility(α, β₁, …, βₙ)          # each a Symbol

Susceptibility of **arbitrary response order** — the `n`-th order term
of the order-parameter response to its conjugate field,

`χ⁽ⁿ⁾_{α; β₁…βₙ} = ∂ⁿ⟨M_α⟩ / ∂h_{β₁}…∂h_{βₙ}`,

a rank-`(n+1)` tensor in [`SpinAxis`](@ref) space whose index parameter
`I = (α, β₁, …, βₙ)` carries one response direction `α` and `n` field
directions.  The response order is `n = length(I) − 1`
([`response_order`](@ref)):

- `Susceptibility(:x, :y)` — **linear** `χ_xy = ∂M_x/∂h_y` (order 1),
  the off-diagonal component the fused `SusceptibilityXX/ZZ` names could
  not express;
- `Susceptibility(:x, :y, :z)` — **second-order nonlinear**
  `χ⁽²⁾_{x;yz} = ∂²M_x/∂h_y∂h_z` (order 2);
- `Susceptibility(:x, :x, :x, :x)` — third-order `χ⁽³⁾`, and so on.

The genealogy is recursive: `χ⁽ⁿ⁾ ⟵ χ⁽ⁿ⁻¹⁾ ⟵ … ⟵ M ⟵ F`
(`derivative_edge`), so `derivative_order(χ⁽ⁿ⁾, MagneticField) == n + 1`.
The linear component's defining identity is [`SusceptibilityFDT`](@ref)
`χ_AB = β·Cov(M_A, M_B)`.
"""
struct Susceptibility{I} <: AbstractSusceptibility
    function Susceptibility{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "Susceptibility needs ≥2 indices (1 response + ≥1 field), got $(repr(I))",
            );
            new{I}()
        )
    end
end
Susceptibility(idx::Symbol...) = Susceptibility{idx}()
tensor_rank(::Type{Susceptibility{I}}) where {I} = length(I)
index_spaces(::Type{Susceptibility{I}}) where {I} = ntuple(_ -> SpinAxis(), length(I))
indices(::Type{Susceptibility{I}}) where {I} = I
response_order(::Type{Susceptibility{I}}) where {I} = length(I) - 1
export Susceptibility

# ─── Two-point correlations (tensors in spin space) ─────────────────────

"""
    SpinCorrelation{A,B}() <: AbstractTwoPointCorrelation
    SpinCorrelation(a::Symbol, b::Symbol)

Two-point spin correlation `⟨S^A_i S^B_j⟩` — a rank-2 tensor in
[`SpinAxis`](@ref) space (`SpinCorrelation(:z, :z)` replaces the old
`ZZCorrelation`).  At criticality its decay is governed by the anomalous
dimension η — see the `correlation_decay` correspondence.
"""
struct SpinCorrelation{A,B} <: AbstractTwoPointCorrelation
    SpinCorrelation{A,B}() where {A,B} = (_axis(A); _axis(B); new{A,B}())
end
SpinCorrelation(a::Symbol, b::Symbol) = SpinCorrelation{a,b}()
tensor_rank(::Type{<:SpinCorrelation}) = 2
index_spaces(::Type{<:SpinCorrelation}) = (SpinAxis(), SpinAxis())
indices(::Type{SpinCorrelation{A,B}}) where {A,B} = (A, B)
export SpinCorrelation

"""
    Conductivity{I}() <: AbstractQuantity
    Conductivity(μ, ν₁, …, νₙ)            # each a Symbol

Electrical conductivity of **arbitrary response order** — the `n`-th
order current response `j_μ = Σ σ⁽ⁿ⁾_{μ; ν₁…νₙ} E_{ν₁}…E_{νₙ}`, a
rank-`(n+1)` tensor in [`SpatialDirection`](@ref) space with one current
direction `μ` and `n` field directions.  `response_order = length(I) − 1`:

- `Conductivity(:x, :y)` — **linear** `σ_xy` (order 1); its Hall
  component is quantized by [`TKNN`](@ref);
- `Conductivity(:x, :y, :z)` — **second-order** `σ⁽²⁾` (nonlinear /
  photogalvanic response), and so on.
"""
struct Conductivity{I} <: AbstractQuantity
    function Conductivity{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "Conductivity needs ≥2 indices (1 current + ≥1 field), got $(repr(I))"
            );
            new{I}()
        )
    end
end
Conductivity(idx::Symbol...) = Conductivity{idx}()
tensor_rank(::Type{Conductivity{I}}) where {I} = length(I)
index_spaces(::Type{Conductivity{I}}) where {I} = ntuple(_ -> SpatialDirection(), length(I))
indices(::Type{Conductivity{I}}) where {I} = I
response_order(::Type{Conductivity{I}}) where {I} = length(I) - 1
export Conductivity

# ─── Dynamical & spectral quantities ────────────────────────────────────
#
# The frequency-resolved family.  These tags name the quantities; the
# identities relating them (`A = −Im G^R/π`, Dyson, DOS = BZ-average of A,
# detailed balance, the NMR relaxation relations) live in
# `relations/spectral.jl` and the transform/sum genealogy in
# `structure/spectral.jl`.  Actually *evaluating* the ω-dependence
# (Kramers–Kronig, analytic continuation) is out of scope here — see the
# functional sibling (issue #14).

"""
    AbstractPropagator <: AbstractQuantity

Single-particle propagators — retarded/advanced/Matsubara Green's
functions and the self-energy — the `(q, ω)`-resolved objects the Dyson
equation relates.
"""
abstract type AbstractPropagator <: AbstractQuantity end
export AbstractPropagator

"""
    RetardedGreensFunction() <: AbstractPropagator

The retarded single-particle Green's function `G^R(q, ω)`.  Its spectral
representation `A = −Im G^R/π` and the Dyson equation
`G^{-1} = G₀^{-1} − Σ` are in `relations/spectral.jl`.
"""
struct RetardedGreensFunction <: AbstractPropagator end
export RetardedGreensFunction
# G_ab(q,ω): rank-2 in orbital space
tensor_rank(::Type{RetardedGreensFunction}) = 2
index_spaces(::Type{RetardedGreensFunction}) = (OrbitalIndex(), OrbitalIndex())
frequency_arguments(::Type{RetardedGreensFunction}) = 1

"""
    SelfEnergy() <: AbstractPropagator

The single-particle self-energy `Σ(q, ω)` — the Dyson correction
`G^{-1} = G₀^{-1} − Σ` between the bare and full propagators.
"""
struct SelfEnergy <: AbstractPropagator end
export SelfEnergy
tensor_rank(::Type{SelfEnergy}) = 2
index_spaces(::Type{SelfEnergy}) = (OrbitalIndex(), OrbitalIndex())
frequency_arguments(::Type{SelfEnergy}) = 1

"""
    SpectralFunction() <: AbstractQuantity

The single-particle spectral function `A(q, ω) = −(1/π) Im G^R(q, ω)`,
normalized by `∫ A(q, ω) dω = 1`.
"""
struct SpectralFunction <: AbstractQuantity end
export SpectralFunction
tensor_rank(::Type{SpectralFunction}) = 2
index_spaces(::Type{SpectralFunction}) = (OrbitalIndex(), OrbitalIndex())
frequency_arguments(::Type{SpectralFunction}) = 1

"""
    DensityOfStates() <: AbstractQuantity

The density of states `ρ(ω) = (1/N) Σ_q A(q, ω)` — the Brillouin-zone
average of the [`SpectralFunction`](@ref).
"""
struct DensityOfStates <: AbstractQuantity end
export DensityOfStates
frequency_arguments(::Type{DensityOfStates}) = 1

"""
    DynamicalCorrelation() <: AbstractTwoPointCorrelation

The space-and-time-resolved correlation `⟨A(r, t) A(0, 0)⟩` whose
space-time Fourier transform (any spatial dimension) is the
[`DynamicalStructureFactor`](@ref).
"""
struct DynamicalCorrelation <: AbstractTwoPointCorrelation end
export DynamicalCorrelation
frequency_arguments(::Type{DynamicalCorrelation}) = 1

"""
    DynamicalStructureFactor() <: AbstractStructureFactor

The dynamical structure factor `S(q, ω)` — the space-time Fourier
transform of the [`DynamicalCorrelation`](@ref); obeys detailed balance
`S(q, −ω) = e^{−βω} S(q, ω)` and the fluctuation–dissipation link to the
[`DynamicalSusceptibility`](@ref).
"""
struct DynamicalStructureFactor <: AbstractStructureFactor end
export DynamicalStructureFactor
# S_αβ(q,ω): rank-2 in spin space
tensor_rank(::Type{DynamicalStructureFactor}) = 2
index_spaces(::Type{DynamicalStructureFactor}) = (SpinAxis(), SpinAxis())
frequency_arguments(::Type{DynamicalStructureFactor}) = 1

"""
    DynamicalSusceptibility{I}() <: AbstractSusceptibility
    DynamicalSusceptibility(α, β₁, …, βₙ)          # each a Symbol

The dynamical susceptibility of **arbitrary response order** — the
frequency-domain (multi-time) counterpart of the static
[`Susceptibility`](@ref).  The linear `DynamicalSusceptibility(:x, :y)`
is `χ_xy(ω)`, one frequency argument, and its imaginary part `χ''(q, ω)`
is the dissipative response of the fluctuation–dissipation theorem and
the NMR relaxation rate.

The `n`-th order term `DynamicalSusceptibility(α, β₁, …, βₙ)` is
`χ⁽ⁿ⁾_{α;β₁…βₙ}(ω₁, …, ωₙ)`: the field is applied at `n` distinct times,
so the response is intrinsically **multi-time** — `frequency_arguments
== n` (`response_order`).  `DynamicalSusceptibility(:x, :y, :z)` is the
second-order `χ⁽²⁾(ω₁, ω₂)` of two-dimensional coherent spectroscopy
(Wan & Armitage, Phys. Rev. Lett. 122, 257401 (2019)).  Its microscopic
Kubo expression is the `n`-fold nested-commutator response function
(Kubo, J. Phys. Soc. Jpn. 12, 570 (1957)); see `structure/spectral.jl`.

The static `Susceptibility{I}` of the same order is the zero-frequency
limit, `χ⁽ⁿ⁾(0, …, 0)`.
"""
struct DynamicalSusceptibility{I} <: AbstractSusceptibility
    function DynamicalSusceptibility{I}() where {I}
        return (
            length(_axistuple(I)) >= 2 || error(
                "DynamicalSusceptibility needs ≥2 indices (1 response + ≥1 field), got $(repr(I))",
            );
            new{I}()
        )
    end
end
DynamicalSusceptibility(idx::Symbol...) = DynamicalSusceptibility{idx}()
tensor_rank(::Type{DynamicalSusceptibility{I}}) where {I} = length(I)
function index_spaces(::Type{DynamicalSusceptibility{I}}) where {I}
    return ntuple(_ -> SpinAxis(), length(I))
end
indices(::Type{DynamicalSusceptibility{I}}) where {I} = I
response_order(::Type{DynamicalSusceptibility{I}}) where {I} = length(I) - 1
# multi-time: an n-th order dynamical response depends on n frequencies
frequency_arguments(::Type{DynamicalSusceptibility{I}}) where {I} = length(I) - 1
export DynamicalSusceptibility

"""
    NMRSpinRelaxationRate() <: AbstractQuantity

The NMR spin–lattice relaxation rate `1/T₁` — set by the low-frequency
limit of the dissipative dynamical susceptibility (Moriya),
`1/T₁ ∝ T · lim_{ω→0} Σ_q |A_hf(q)|² χ''(q, ω)/ω`.
"""
struct NMRSpinRelaxationRate <: AbstractQuantity end
export NMRSpinRelaxationRate

"""
    NMRRelaxationExponent() <: AbstractQuantity

The low-temperature scaling exponent `θ_NMR` of `1/T₁ ∝ T^{θ_NMR}`,
fixed by the operator scaling dimension via `θ_NMR = 2Δ_op − 1`.
"""
struct NMRRelaxationExponent <: AbstractQuantity end
export NMRRelaxationExponent

"""
    ScalingDimension() <: AbstractQuantity

The scaling dimension `Δ_op` of a local operator at a quantum critical
point — the input to dynamical scaling relations such as the NMR
exponent `θ_NMR = 2Δ_op − 1`.
"""
struct ScalingDimension <: AbstractQuantity end
export ScalingDimension

# ─── Criticality ────────────────────────────────────────────────────────

"""
    CriticalTemperature() <: AbstractQuantity

Critical temperature `T_c` of a finite-temperature phase transition.

Generic home for a tag that previously lived inside a model file.
"""
struct CriticalTemperature <: AbstractQuantity end
export CriticalTemperature

"""
    CorrelationLength() <: AbstractQuantity

Correlation length ξ (units of lattice spacing).
"""
struct CorrelationLength <: AbstractQuantity end
export CorrelationLength

"""
    UniversalityClass() <: AbstractQuantity

The universality class a model's transition belongs to (returned as a
[`Universality`](@ref) tag by implementing atlases).
"""
struct UniversalityClass <: AbstractQuantity end
export UniversalityClass

"""
    CentralCharge() <: AbstractQuantity

Central charge `c` of the critical theory's CFT.
"""
struct CentralCharge <: AbstractQuantity end
export CentralCharge

# ─── Topology ───────────────────────────────────────────────────────────

"""
    TopologicalInvariant() <: AbstractQuantity

The model's topological invariant (winding number, Chern number, ℤ₂
index, Pfaffian sign, … — the concrete meaning is declared by the
implementing model).  Generic *computations* of standard invariants on
Bloch maps live in this package's relations layer
([`winding_number`](@ref), [`chern_number`](@ref)).
"""
struct TopologicalInvariant <: AbstractQuantity end
export TopologicalInvariant
