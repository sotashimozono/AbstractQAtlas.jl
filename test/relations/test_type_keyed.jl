# test/relations/test_type_keyed.jl — the type-keyed relation interface.
#
# The prototype migrates the undecorated Green's-function relations (Dyson,
# SpectralFromGreens, the six Keldysh relations) to TYPE keys.  These tests keep
# the PHYSICS assertions of the symbol-keyed suites and add: (a) that the
# auto-derived `quantities` reproduces the deleted hand-links exactly, (b) that
# the type-keyed bag / residual / check / solve / report agree with the kernels,
# (c) the typed β-or-T conversion, and (d) the collision-freedom that motivated
# the redesign — distinct quantity types can never share a key.
# See docs/design/type-keyed-interface.md.

using AbstractQAtlas
using Test
const AQ = AbstractQAtlas

@testset "VariableKey + bag basics" begin
    @test VariableKey(RetardedGreensFunction) ==
        VariableKey(RetardedGreensFunction, Global())
    @test VariableKey(RetardedGreensFunction) != VariableKey(AdvancedGreensFunction)
    @test hash(VariableKey(SpectralFunction)) == hash(VariableKey(SpectralFunction))
    b = bag(SpectralFunction => 1.0, RetardedGreensFunction => 2.0 + im)
    @test b isa AQ.Bag
    @test b[VariableKey(SpectralFunction)] == 1.0
    @test length(b) == 2
    # a bare Type key and an explicit VariableKey key are the same slot
    @test bag(SelfEnergy => 3)[VariableKey(SelfEnergy)] == 3
end

@testset "auto-derived quantities reproduce the deleted hand-links" begin
    # the exact tuples that were in quantity_links.jl before migration — the
    # proof that deleting them is behavior-preserving (order-insensitive).
    expected = Dict(
        Dyson() => (RetardedGreensFunction, SelfEnergy),
        SpectralFromGreens() => (SpectralFunction, RetardedGreensFunction),
        KeldyshComponent() =>
            (KeldyshGreensFunction, GreaterGreensFunction, LesserGreensFunction),
        KeldyshCausality() => (
            RetardedGreensFunction,
            AdvancedGreensFunction,
            GreaterGreensFunction,
            LesserGreensFunction,
        ),
        AdvancedRetardedConjugate() => (AdvancedGreensFunction, RetardedGreensFunction),
        KeldyshFDT() =>
            (KeldyshGreensFunction, RetardedGreensFunction, AdvancedGreensFunction),
        KMSGreaterLesser() => (LesserGreensFunction, GreaterGreensFunction),
        SpectralFromKeldysh() =>
            (SpectralFunction, RetardedGreensFunction, AdvancedGreensFunction),
    )
    for (rel, qs) in expected
        @test Set(quantities(rel)) == Set(qs)
        # every listed subject is a real AbstractQuantity leaf/family
        @test all(T -> T <: AbstractQuantity, quantities(rel))
        # the reverse index resolves (registry queryable without hand-links)
        for T in quantities(rel)
            @test rel in relations_constraining(T())
        end
    end
    # fields are NOT quantities: KMS's β::InverseTemperature is excluded
    @test !(InverseTemperature in quantities(KMSGreaterLesser()))
    # supplied (untyped) slots are excluded: Dyson's bare G₀ is not a quantity
    @test length(quantities(Dyson())) == 2                # G and Σ, not G₀
end

@testset "type-keyed residual / check on the Green's identities" begin
    # KeldyshComponent  G^K = G^> + G^<   (independent numbers)
    @test residual(
        KeldyshComponent(),
        bag(
            KeldyshGreensFunction => 5.0,
            GreaterGreensFunction => 2.0,
            LesserGreensFunction => 3.0,
        ),
    ) == 0.0
    @test !check(
        KeldyshComponent(),
        bag(
            KeldyshGreensFunction => 5.0,
            GreaterGreensFunction => 2.0,
            LesserGreensFunction => 3.5,
        ),
    )

    # a single-pole retarded propagator; the whole RAK web is consistent
    ε, η = 0.4, 0.05
    for ω in (-0.3, 0.0, 0.9)
        GR = 1 / (ω - ε + im * η)
        GA = conj(GR)                      # AdvancedRetardedConjugate
        Aw = -imag(GR) / π                 # A = −Im G^R/π
        Gles = 0.1 - 0.2im                 # arbitrary; only G^> − G^< is fixed…
        Ggtr = Gles + (GR - GA)            # …so that G^> − G^< = G^R − G^A (KeldyshCausality)
        @test check(
            AdvancedRetardedConjugate(),
            bag(AdvancedGreensFunction => GA, RetardedGreensFunction => GR);
            atol=1e-12,
        )
        @test check(
            SpectralFromGreens(),
            bag(SpectralFunction => Aw, RetardedGreensFunction => GR);
            atol=1e-12,
        )
        # SpectralFromKeldysh reduces to the same A once G^A = conj(G^R)
        @test check(
            SpectralFromKeldysh(),
            bag(
                SpectralFunction => Aw,
                RetardedGreensFunction => GR,
                AdvancedGreensFunction => GA,
            );
            atol=1e-12,
        )
        # KeldyshCausality  G^R − G^A = G^> − G^<
        @test check(
            KeldyshCausality(),
            bag(
                RetardedGreensFunction => GR,
                AdvancedGreensFunction => GA,
                GreaterGreensFunction => Ggtr,
                LesserGreensFunction => Gles,
            );
            atol=1e-12,
        )
    end
end

@testset "supplied (untyped) slots via extras" begin
    # Dyson: G₀ is a supplied slot (no free-propagator quantity type)
    G0, Σ = 1 / (0.5 + im * 0.1), 0.2 + 0.05im
    G = inv(inv(G0) - Σ)               # construct the full G so the identity holds
    @test abs(residual(Dyson(), bag(RetardedGreensFunction => G, SelfEnergy => Σ); G0=G0)) <
        1e-12
    # KeldyshFDT: the distribution h is a supplied slot
    GR = 1 / (0.3 + im * 0.1)
    GA = conj(GR)
    h = 2.3
    GK = h * (GR - GA)
    @test abs(
        residual(
            KeldyshFDT(),
            bag(
                KeldyshGreensFunction => GK,
                RetardedGreensFunction => GR,
                AdvancedGreensFunction => GA,
            );
            h=h,
        ),
    ) < 1e-12
    # a missing supplied slot is a loud, type/name-named error, not a silent skip
    @test_throws ErrorException residual(
        Dyson(), bag(RetardedGreensFunction => G, SelfEnergy => Σ)
    )
end

@testset "type-keyed solve preserves exact arithmetic" begin
    @test solve(
        KeldyshComponent(),
        KeldyshGreensFunction,
        bag(GreaterGreensFunction => 2 // 1, LesserGreensFunction => 3 // 1),
    ) == 5 // 1
    # solve for a summand
    @test solve(
        KeldyshComponent(),
        GreaterGreensFunction,
        bag(KeldyshGreensFunction => 5 // 1, LesserGreensFunction => 3 // 1),
    ) == 2 // 1
    # asking for a type the relation does not carry is a loud error
    @test_throws ErrorException solve(
        KeldyshComponent(), SelfEnergy, bag(GreaterGreensFunction => 1)
    )
end

@testset "typed β-or-T (InverseTemperature ⇄ Temperature)" begin
    # KMSGreaterLesser  G^< = ζ e^{−βω} G^>   (fermion ζ = −1)
    Ggtr, ζ, ω, β = 2.0, -1, 1.0, 0.5
    Gles = ζ * exp(-β * ω) * Ggtr
    rβ = residual(
        KMSGreaterLesser(),
        bag(
            LesserGreensFunction => Gles,
            GreaterGreensFunction => Ggtr,
            InverseTemperature => β,
        );
        ζ=ζ,
        ω=ω,
    )
    rT = residual(
        KMSGreaterLesser(),
        bag(
            LesserGreensFunction => Gles,
            GreaterGreensFunction => Ggtr,
            Temperature => 1 / β,
        );
        ζ=ζ,
        ω=ω,
    )
    @test rβ == rT
    @test abs(rβ) < 1e-12
end

@testset "type-keyed bag report / applicability" begin
    GR = 1 / (0.2 + im * 0.1)
    GA = conj(GR)
    b = bag(
        RetardedGreensFunction => GR,
        AdvancedGreensFunction => GA,
        SpectralFunction => -imag(GR) / π,
    )
    matched = Set(nameof(typeof(r)) for r in applicable_relations(b))
    # exactly the type-keyed relations whose identity types are all present
    @test :AdvancedRetardedConjugate in matched
    @test :SpectralFromKeldysh in matched
    @test :SpectralFromGreens in matched
    # a relation needing a type absent from the bag is NOT matched
    @test !(:KeldyshComponent in matched)     # needs G^K/G^>/G^<
    # every matched relation passes on this consistent bag
    @test check_all(b; atol=1e-12)
    rep = relation_report(b; atol=1e-12)
    @test !isempty(rep) && all(row -> row.pass, rep)
end

@testset "collision-freedom: the redesign's whole point" begin
    # The old symbol keying let `:G`/`:GR` drift and `:S` mean four things.
    # With type keys, G^R and G^A are distinct slots that cannot be confused:
    GR = 1 / (0.2 + im * 0.1)
    GA = conj(GR)
    Aw = -imag(GR) / π
    good = bag(
        SpectralFunction => Aw, RetardedGreensFunction => GR, AdvancedGreensFunction => GA
    )
    @test check(SpectralFromKeldysh(), good; atol=1e-12)
    # swap the two propagators' values → the identity now fails, proving the two
    # TYPE keys are honored as distinct roles (a symbol fusion could not tell them apart)
    swapped = bag(
        SpectralFunction => Aw, RetardedGreensFunction => GA, AdvancedGreensFunction => GR
    )
    @test !check(SpectralFromKeldysh(), swapped; atol=1e-12)

    # a legacy symbol-only relation is NEVER pulled into a type-keyed bag report,
    # so a Green's bag can't accidentally trigger, say, a thermoelectric relation
    @test isempty(variable_types(KelvinRelation()))              # not migrated
    @test !(KelvinRelation() in applicable_relations(good))
end
