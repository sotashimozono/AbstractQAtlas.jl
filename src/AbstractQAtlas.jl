"""
    AbstractQAtlas

The model-independent layer of the QAtlas ecosystem, in the spirit of
`AbstractFFTs`: concrete atlases (QAtlas) *implement* this package,
never the reverse.

Two responsibilities:

1. **core/** — the abstract type vocabulary for physical quantities:
   `AbstractQAtlasModel`, `BoundaryCondition` (`Infinite`/`OBC`/`PBC`),
   `AbstractQuantity` and its hierarchy, the generic
   [`fetch`](@ref AbstractQAtlas.fetch) verb, and the
   [`Universality`](@ref) machinery — so atlases and third packages
   share dispatch types without depending on a full atlas.

2. **relations/** — generic, model-independent physics relations as
   first-class tested objects with a uniform three-verb interface
   ([`residual`](@ref) / [`check`](@ref) / [`solve`](@ref)):
   critical-exponent scaling laws, fluctuation–dissipation identities,
   Wick's theorem, standard topological invariants, and finite-size
   scaling forms.

**Values do not live here, and neither do model-specific laws.**  This
package owns only what holds **universally within a domain** — independent
of the system's symmetry, Hamiltonian, or any individual detail (Wick's
theorem, the fluctuation–dissipation and Maxwell relations, the entropy
inequalities, Kramers–Kronig, the scaling laws, …).  Reference *numbers*
(critical temperatures, exact magnetizations, exponent tables) AND
*model-specific* relations (the Drude mobility `μ=eτ/m`, the ±J Nishimori-
line energy, the SK de Almeida–Thouless line, single-band `R_H=1/ne`, …)
belong to the implementing atlas (QAtlas), not here.  The library is a
**universal yardstick**: apply its relations to measured quantities to
check whether a system obeys the laws that must hold regardless of its
details.
"""
module AbstractQAtlas

using LinearAlgebra: det, eigen, Hermitian

# core — type vocabulary
include("core/types.jl")
include("core/indices.jl")
include("core/representations.jl")
include("core/quantities.jl")
include("core/universality.jl")
include("core/distributions.jl")
include("core/fields.jl")

# structure — model-independent definitional correspondences between the
# core quantities: transition classification, the critical
# quantity⇄exponent map from which the scaling forms are derived, and the
# response-function derivative genealogy rooted at the free energy.
include("structure/transitions.jl")
include("structure/criticality.jl")
include("structure/response.jl")
include("structure/tensor_symmetry.jl")
include("structure/spectral.jl")
include("structure/fourier.jl")

# relations — model-independent identities and forms
include("relations/interface.jl")
include("relations/scaling.jl")
include("relations/thermodynamic.jl")
include("relations/fundamental.jl")
include("relations/statistics.jl")
include("relations/wick.jl")
include("relations/topology.jl")
include("relations/spectral.jl")
include("relations/transport.jl")
include("relations/quantum.jl")
include("relations/ensembles.jl")
include("relations/entanglement.jl")
include("relations/cft.jl")

# automatic-differentiation entry point (methods live in ext/, ForwardDiff)
include("autodiff.jl")

end # module AbstractQAtlas
