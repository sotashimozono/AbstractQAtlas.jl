# The unified quantity-relationship graph: related_quantities / quantity_graph
# / quantity_path / quantity_graph_jsonl fold the four separate typed-edge
# accessors (derivative_edge / spectral_origin / fourier_conjugate_quantity +
# the relation↔quantity links) into ONE graph.  These tests check that the
# unification is FAITHFUL (no edge invented, none dropped), that paths are
# well-formed, and that the export is consistent — against the independent
# underlying accessors, not the graph's own output.

using AbstractQAtlas
using AbstractQAtlas:
    QuantityEdge,
    related_quantities,
    quantity_graph,
    quantity_neighbors,
    quantity_path,
    quantity_graph_jsonl,
    derivative_edge,
    spectral_origin,
    fourier_conjugate_quantity,
    relations_constraining,
    quantities,
    _family,
    _rep_quantity,
    AbstractQuantity

@testset "related_quantities faithfully mirrors the underlying accessors" begin
    # every :derivative / :spectral / :fourier edge from related_quantities is
    # exactly what the standalone accessor says, and vice versa — no invention.
    for T in (
        Susceptibility{(:z, :z)},
        SpecificHeat,
        Energy,
        SpectralFunction,
        DensityOfStates,
        StaticStructureFactor,
        Magnetization{:z},
    )
        q = _rep_quantity(T)
        q === nothing && continue
        edges = related_quantities(q)

        de = derivative_edge(q)
        if de === nothing
            @test !any(e -> e.kind === :derivative, edges)
        else
            d = only(filter(e -> e.kind === :derivative, edges))
            @test d.to === _family(de.parent)
        end

        so = spectral_origin(q)
        if so === nothing
            @test !any(e -> e.kind === :spectral, edges)
        else
            s = only(filter(e -> e.kind === :spectral, edges))
            @test s.to === _family(so.from)
            @test s.detail == string(so.via)
        end

        fc = fourier_conjugate_quantity(typeof(q))
        if fc === nothing
            @test !any(e -> e.kind === :fourier, edges)
        else
            f = only(filter(e -> e.kind === :fourier, edges))
            @test f.to === _family(fc)
        end

        # every :law edge corresponds to a real shared relation naming both
        for e in filter(e -> e.kind === :law, edges)
            rels = relations_constraining(q)
            @test any(
                r ->
                    string(nameof(typeof(r))) == e.detail &&
                    any(U -> _family(U) === e.to, quantities(r)),
                rels,
            )
        end
    end
end

@testset "known physics edges are present (independent expectations)" begin
    # χ_zz is tied to M_z BOTH as a field-derivative AND by the FDT — two
    # independent edge kinds between the same pair.
    chi = related_quantities(Susceptibility(:z, :z))
    @test any(e -> e.kind === :derivative && e.to === Magnetization, chi)
    @test any(
        e -> e.kind === :law && e.to === Magnetization && e.detail == "SusceptibilityFDT",
        chi,
    )
    # A(ω) = −Im G/π  is a :spectral edge A → G
    aw = related_quantities(SpectralFunction())
    @test any(
        e ->
            e.kind === :spectral &&
            e.to === RetardedGreensFunction &&
            e.detail == "neg_im_over_pi",
        aw,
    )
    # S(q) ↔ ⟨SS⟩(r) is a Fourier pair
    sq = related_quantities(StaticStructureFactor())
    @test any(e -> e.kind === :fourier && e.to === SpinCorrelation, sq)
end

@testset "quantity_graph is deduplicated and well-typed" begin
    g = quantity_graph()
    @test !isempty(g)
    # all four edge kinds are represented (nothing silently dropped)
    for k in (:derivative, :spectral, :fourier, :law)
        @test any(e -> e.kind === k, g)
    end
    # every endpoint is a genuine quantity family
    for e in g
        @test e.from <: AbstractQuantity
        @test e.to <: AbstractQuantity
        @test e.from !== e.to                         # no self-loops
    end
    # no duplicate unordered {kind, {from,to}, detail} edge
    keys = Set{Any}()
    for e in g
        k1 = (e.kind, e.from, e.to, e.detail)
        k2 = (e.kind, e.to, e.from, e.detail)
        @test !(k1 in keys) && !(k2 in keys)
        push!(keys, k1)
    end
    # caching: a second call returns the identical object
    @test quantity_graph() === g
end

@testset "quantity_path finds well-formed shortest paths" begin
    # same family → empty path
    @test quantity_path(SpecificHeat(), SpecificHeat()) == QuantityEdge[]
    # a connected pair: SpecificHeat and Magnetization share the FreeEnergy root
    p = quantity_path(SpecificHeat(), Magnetization(:z))
    @test p !== nothing
    @test !isempty(p)
    # every edge on the path is a real graph edge
    g = quantity_graph()
    for e in p
        @test any(
            x ->
                x.kind === e.kind &&
                x.detail == e.detail &&
                (
                    (x.from === e.from && x.to === e.to) ||
                    (x.from === e.to && x.to === e.from)
                ),
            g,
        )
    end
    # consecutive edges share a node (the chain is actually connected)
    for i in 1:(length(p) - 1)
        shared = Set((p[i].from, p[i].to)) ∩ Set((p[i + 1].from, p[i + 1].to))
        @test !isempty(shared)
    end
    # node-incidence over a simple path: the two terminal nodes (odd incidence)
    # are exactly the requested endpoints; every internal node has even
    # incidence — orientation-independent, so it doesn't depend on how each
    # accessor happened to orient its edge.
    incidence = Dict{Type,Int}()
    for e in p
        incidence[e.from] = get(incidence, e.from, 0) + 1
        incidence[e.to] = get(incidence, e.to, 0) + 1
    end
    terminals = Set(n for (n, c) in incidence if isodd(c))
    @test terminals == Set((SpecificHeat, Magnetization))
end

@testset "quantity_neighbors surfaces incident edges (both directions)" begin
    # FreeEnergy is a derivative ROOT: it never appears as a `from` of a
    # :derivative edge, but several quantities point AT it — quantity_neighbors
    # must still surface those (related_quantities alone would not).
    nb = quantity_neighbors(FreeEnergy)
    @test !isempty(nb)
    @test all(e -> e.from === FreeEnergy || e.to === FreeEnergy, nb)
    @test any(e -> e.kind === :derivative && e.to === FreeEnergy, nb)
end

@testset "quantity_graph_jsonl export is consistent" begin
    io = IOBuffer()
    quantity_graph_jsonl(io)
    lines = split(strip(String(take!(io))), '\n')
    g = quantity_graph()
    nodes = Set{Type}()
    for e in g
        push!(nodes, e.from)
        push!(nodes, e.to)
    end
    # summary line first
    @test occursin("\"nodes\":$(length(nodes))", lines[1])
    @test occursin("\"edges\":$(length(g))", lines[1])
    # one edge object per graph edge, each carrying the four fields
    @test length(lines) == length(g) + 1
    for ln in lines[2:end]
        @test occursin("\"kind\":", ln)
        @test occursin("\"from\":", ln)
        @test occursin("\"to\":", ln)
        @test occursin("\"detail\":", ln)
    end
end
