# test/report/test_card.jl — the report/card contract (schema-v2).
#
# `report(model, quantity, bc; value, …)` packages an oracle's computed value into a
# schema-valid `Card` (the reporter-facing sibling of `fetch`); `card_jsonl` streams
# them NaN-safely.  AbstractQAtlas is model-independent, so the tests define their
# own trivial model to exercise the hub construction.

using AbstractQAtlas
using Test
const AQ = AbstractQAtlas

struct _DemoModel <: AbstractQAtlasModel end

@testset "report: hub, status, provenance, error_bar" begin
    c = report(
        _DemoModel(), VonNeumannEntropy(), PBC(64);
        value=0.87, err=0.01, route=:monte_carlo,
        provenance="ClassicalMonteCarlo@metropolis", mechanism="L=64",
        atol=1e-3, refs=["Calabrese2004"],
    )
    @test c isa Card
    @test c.schema_version == 2
    @test c.hub == "_DemoModel/VonNeumannEntropy/PBC"      # TypeName(model)/TypeName(quantity)/TypeName(bc)
    @test c.status == :ok && c.subject == 0.87 && c.error_bar == 0.01
    @test c.independence == :measured                      # a bare measurement, no cross-check
    @test c.provenance == "ClassicalMonteCarlo@metropolis"
    @test c.mechanism == "L=64" && c.atol == 1e-3 && c.refs == ["Calabrese2004"]

    # instances and types give the SAME hub (report accepts either)
    @test report(_DemoModel, VonNeumannEntropy, Infinite; value=1.0, route=:dmrg,
        provenance="x").hub == "_DemoModel/VonNeumannEntropy/Infinite"
end

@testset "report: non-finite → :divergent, never a raw token" begin
    for bad in (NaN, Inf, -Inf)
        c = report(_DemoModel(), VonNeumannEntropy(), Infinite();
            value=bad, route=:monte_carlo, provenance="p")
        @test c.status == :divergent && c.subject === nothing
    end
    # Hermitian-round-off imaginary part → real; a genuinely complex value → divergent
    @test report(_DemoModel(), VonNeumannEntropy(), Infinite(); value=1.0 + 1e-12im,
        route=:ed_finite_size, provenance="p").subject == 1.0
    @test report(_DemoModel(), VonNeumannEntropy(), Infinite(); value=1.0 + 1.0im,
        route=:ed_finite_size, provenance="p").status == :divergent
end

@testset "report: independence class from route + cross-check" begin
    mk(route, ind) = report(_DemoModel(), VonNeumannEntropy(), Infinite();
        value=1.0, route=route, provenance="p", independent=ind)
    @test mk(:monte_carlo, ()).independence == :measured               # no cross-check
    @test mk(:ed_finite_size, (1.0, 1.0)).independence == :structural  # mechanically independent
    @test mk(:sum_rule, (1.0,)).independence == :asserted              # corroborated but asserted
    # a non-finite cross-check corroborates nothing (dropped from independent[])
    @test isempty(mk(:ed_finite_size, (NaN,)).independent)
    # unknown route rejected loudly
    @test_throws ErrorException report(_DemoModel(), VonNeumannEntropy(), Infinite();
        value=1.0, route=:made_up, provenance="p")
end

@testset "card_jsonl: schema-valid, NaN-safe JSONL" begin
    cards = [
        report(_DemoModel(), VonNeumannEntropy(), PBC(64); value=0.5, err=0.01,
            route=:monte_carlo, provenance="mc", refs=["A2020"]),
        report(_DemoModel(), VonNeumannEntropy(), Infinite(); value=NaN,
            route=:dmrg, provenance="d"),
    ]
    io = IOBuffer()
    card_jsonl(io, cards)
    lines = split(strip(String(take!(io))), '\n')
    @test length(lines) == 2
    # first line: the ok card carries its fields
    @test occursin("\"schema_version\":2", lines[1])
    @test occursin("\"hub\":\"_DemoModel/VonNeumannEntropy/PBC\"", lines[1])
    @test occursin("\"provenance\":\"mc\"", lines[1])
    @test occursin("\"refs\":[\"A2020\"]", lines[1])
    @test occursin("\"route\":\"monte_carlo\"", lines[1])
    # second line: the divergent card emits null, never a raw NaN token
    @test !occursin("NaN", lines[2])
    @test occursin("\"subject\":null", lines[2])
    @test occursin("\"status\":\"divergent\"", lines[2])
end
