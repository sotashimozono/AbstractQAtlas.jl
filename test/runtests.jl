ENV["GKSwstype"] = "100"

using AbstractQAtlas
using Test, Aqua
const dirs = ["core", "structure", "relations", "report", "ext"]

const FIG_BASE = joinpath(pkgdir(AbstractQAtlas), "docs", "src", "assets")
const PATHS = Dict()
mkpath.(values(PATHS))

# The runnable CASES: "aqua" (the module's own Aqua suite) + every dir/test_*.jl file.
# `test/ci/plan_shards.jl` splits this SAME list across CI shards and hands each shard
# its subset in the `CI_CASES` env var (space-separated).  An empty `CI_CASES` — a
# local run, or the reusable's unsharded fallback — runs everything.
function all_cases()
    cases = String["aqua"]
    for dir in dirs
        dpath = joinpath(@__DIR__, dir)
        isdir(dpath) || continue
        for f in
            sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dpath)))
            push!(cases, "$dir/$f")
        end
    end
    return cases
end

const _CI_CASES = String.(split(strip(get(ENV, "CI_CASES", ""))))
const CASES = isempty(_CI_CASES) ? all_cases() : _CI_CASES

@testset "tests" begin
    println("Passed arguments ARGS = $(copy(ARGS)) to tests.")
    println("Running $(length(CASES)) case(s): $(join(CASES, ", "))")
    @test !isempty(CASES)                       # never a silent empty run (bad CI_CASES)
    @time for case in CASES
        if case == "aqua"
            @testset "Aqua tests" begin
                Aqua.test_all(AbstractQAtlas)
            end
        else
            @testset "$case" begin
                @time begin
                    println("  Including $(case)")
                    include(joinpath(@__DIR__, case))
                end
            end
        end
    end
end
