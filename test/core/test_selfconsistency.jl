# Self-consistency of the relation ↔ quantity web: every declared link
# resolves, the reverse index is consistent, the accessibility guarantee
# holds, and the quantity-graph edges have no dangling nodes.

using AbstractQAtlas
using AbstractQAtlas:
    quantities,
    relations_constraining,
    all_relations,
    AbstractQuantity,
    AbstractRelation,
    derivative_edge,
    spectral_origin,
    variables

@testset "every relation → quantity link resolves to a real quantity leaf" begin
    for r in all_relations()
        for T in quantities(r)
            @test T isa Type
            @test T <: AbstractQuantity            # no dangling link to a non-quantity
        end
    end
    # at least a substantial fraction is annotated (the web is populated)
    @test count(r -> !isempty(quantities(r)), all_relations()) >= 40
end

@testset "reverse index `relations_constraining` matches `quantities`" begin
    # the reverse index is exactly the relations whose `quantities` contain the type
    for q in (
        Susceptibility(:z, :z),
        SpecificHeat(),
        RetardedGreensFunction(),
        VonNeumannEntropy(),
        Conductivity(:x, :y),
    )
        got = Set(typeof.(relations_constraining(q)))
        want = Set(
            typeof(r) for r in all_relations() if any(T -> typeof(q) <: T, quantities(r))
        )
        @test got == want
        # instance and type queries agree
        @test Set(typeof.(relations_constraining(typeof(q)))) == got
    end
end

@testset "accessibility guarantee: key quantities reach their laws" begin
    # from a quantity you can find the universal laws it must obey
    susc = Set(typeof.(relations_constraining(Susceptibility(:z, :z))))
    @test SusceptibilityFDT in susc            # χ = β·Var(M)
    @test SusceptibilityResponse in susc        # χ = ∂M/∂h
    @test Dyson in Set(typeof.(relations_constraining(RetardedGreensFunction())))
    @test WiedemannFranz in
        Set(typeof.(relations_constraining(ThermalConductivity(:x, :x))))
    @test KitaevPreskillTEE in
        Set(typeof.(relations_constraining(TopologicalEntanglementEntropy())))
    # round-trip: each annotated relation is found by each of its quantities
    for r in all_relations(), T in quantities(r)
        q = try
            T()
        catch
            try
                T(:x)
            catch
                try
                    T(:x, :y)
                catch
                    nothing
                end
            end
        end
        q === nothing && continue
        @test r in relations_constraining(q)
    end
end

@testset "cross-relation coherence: Susceptibility ↔ Magnetization two ways" begin
    # χ appears in the FDT (χ=β·Var(M)) AND the response (χ=∂M/∂h); the
    # genealogy independently says χ is a field-derivative of M — the
    # relation web and the derivative graph must name the SAME neighbour.
    @test Magnetization in quantities(SusceptibilityFDT())      # FDT links χ to M
    @test derivative_edge(Susceptibility(:z, :z)).parent == Magnetization{:z}  # graph: χ ⟵ M
    # SpecificHeat's fluctuation (β²Var(E)) and its ≥ 0 stability agree on the quantity
    @test SpecificHeat in quantities(SpecificHeatFDT())
    @test SpecificHeat in quantities(SpecificHeatPositivity())
end

@testset "quantity-graph edges have no dangling nodes" begin
    using InteractiveUtils: subtypes
    _leaves(T, acc=Type[]) = (
        foreach(S -> isabstracttype(S) ? _leaves(S, acc) : push!(acc, S), subtypes(T));
        acc
    )
    for T in _leaves(AbstractQuantity)
        q = try
            if T === ThermalAverage
                ThermalAverage(Susceptibility(:x, :y), Canonical(1.0))
            else
                T()
            end
        catch
            try
                T(:x)
            catch
                try
                    T(:x, :y)
                catch
                    nothing
                end
            end
        end
        q === nothing && continue
        # every response-genealogy edge points at a real quantity type
        e = derivative_edge(q)
        e === nothing || @test e.parent isa Type
        # every spectral-graph edge's source is a real quantity type
        o = spectral_origin(q)
        o === nothing || @test o.from isa Type
    end
end
