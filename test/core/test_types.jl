# Core type vocabulary: subtyping, BC size semantics, fetch fallback.

using AbstractQAtlas
using AbstractQAtlas: _bc_size, fetch

@testset "type tree" begin
    @test Universality(:Ising) isa AbstractQAtlasModel
    @test Universality(:Ising) === Universality{:Ising}()
    @test Infinite() isa BoundaryCondition
    @test OBC(4) isa BoundaryCondition
    @test PBC(4) isa BoundaryCondition
    @test CriticalExponents() isa AbstractQuantity
    @test SpecificHeat() isa AbstractThermalPotential
    @test SpecificHeat() isa AbstractQuantity
    @test MagnetizationZ() isa AbstractMagnetization
    @test SusceptibilityZZ() isa AbstractSusceptibility
    @test PartitionFunction() isa AbstractThermalPotential
    @test SpontaneousMagnetization() isa AbstractMagnetization
    @test CriticalTemperature() isa AbstractQuantity
    @test TopologicalInvariant() isa AbstractQuantity
end

@testset "component trait" begin
    @test component(MagnetizationX()) == :x
    @test component(MagnetizationY()) == :y
    @test component(MagnetizationZ()) == :z
    @test component(SusceptibilityXX()) == :xx
    @test component(SusceptibilityZZ()) == :zz
    # default: no component
    @test component(SpecificHeat()) === nothing
    @test component(CriticalExponents()) === nothing
end

@testset "Energy granularity" begin
    @test Energy() === Energy{:natural}()
    @test Energy(:total) === Energy{:total}()
    @test Energy(:per_site) === Energy{:per_site}()
    @test_throws ErrorException Energy(:bogus)
end

@testset "BC size semantics (N = 0 sentinel)" begin
    @test _bc_size(OBC(12), NamedTuple()) == 12
    @test _bc_size(PBC(8), NamedTuple()) == 8
    # sentinel: N from kwargs
    @test _bc_size(OBC(), (N=24,)) == 24
    @test _bc_size(PBC(; N=0), (N=6,)) == 6
    # bc.N wins over kwargs
    @test _bc_size(OBC(12), (N=99,)) == 12
    # unresolvable
    @test_throws ErrorException _bc_size(OBC(), NamedTuple())
    @test_throws ErrorException _bc_size(Infinite(), (N=4,))
end

@testset "fetch fallback throws informatively" begin
    struct _TTModel <: AbstractQAtlasModel end
    err = try
        fetch(_TTModel(), SpecificHeat(), Infinite())
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("_TTModel", err.msg)
    @test occursin("SpecificHeat", err.msg)
    @test occursin("Infinite", err.msg)
end
