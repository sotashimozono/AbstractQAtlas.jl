# ext/AbstractQAtlasZygoteExt.jl — reverse-mode AD realization of the first-order
# response genealogy via Zygote.
#
# `thermal_gradient(F, x⃗)` returns the whole first-order response vector −∇F in a
# single reverse pass — the reverse-mode companion of the per-component,
# forward-mode `thermal_derivative` (the ForwardDiff extension).  For a free
# energy of a high-dimensional field vector this gets every component
# (`M_α = −∂F/∂h_α`) at once, where forward mode would take one pass per
# component.  Same generic verb + erroring-fallback seam as `thermal_derivative`
# / `fetch` (core-functions.md, pillars 2 & 4): the interface is in the parent,
# the backend at the leaf.

module AbstractQAtlasZygoteExt

using AbstractQAtlas
import AbstractQAtlas: thermal_gradient        # extended below → must import
using Zygote: gradient

# −∇F: the full first-order response (a vector for a field vector, a scalar for a
# scalar field) in one reverse-mode pass.  `Zygote.gradient` returns a 1-tuple.
# Typed on the point argument so these ADD methods to the generic `thermal_gradient`
# rather than overwriting its erroring fallback (which would break precompilation).
thermal_gradient(F, x::AbstractVector) = -gradient(F, x)[1]
thermal_gradient(F, x::Number) = -gradient(F, x)[1]

end # module AbstractQAtlasZygoteExt
