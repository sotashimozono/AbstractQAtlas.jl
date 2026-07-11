# Topological invariants vs KNOWN phase structures.
#
# SSH winding and QWZ Chern regimes are unambiguous, independently known
# results; FHS returns exact integers already on coarse grids (that is
# the method's selling point), so the assertions are exact.

using AbstractQAtlas
using AbstractQAtlas: residual, check, solve
using Test

@testset "SSH winding number" begin
    ssh(v, w) = k -> (v + w * cos(k), w * sin(k))
    # topological: v < w winds once (counterclockwise ⇒ +1 in our sign)
    @test winding_number(ssh(0.0, 1.0)) == 1
    @test winding_number(ssh(0.3, 1.0)) == 1
    # trivial: v > w does not enclose the origin
    @test winding_number(ssh(1.0, 0.3)) == 0
    @test winding_number(ssh(2.0, 1.0)) == 0
    # near-critical but resolved
    @test winding_number(ssh(0.9, 1.0)) == 1
    @test winding_number(ssh(1.1, 1.0)) == 0
    # gap closing v = w is rejected
    @test_throws ErrorException winding_number(ssh(1.0, 1.0))
end

@testset "QWZ Chern number (FHS on a coarse grid)" begin
    σx = [0.0 1.0; 1.0 0.0]
    σy = [0.0 -im; im 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    qwz(m) = (kx, ky) -> sin(kx) .* σx .+ sin(ky) .* σy .+ (m + cos(kx) + cos(ky)) .* σz

    c_p1 = chern_number(qwz(1.0), 1)    # 0 <  m < 2
    c_m1 = chern_number(qwz(-1.0), 1)   # −2 < m < 0
    # regime structure of the QWZ phase diagram (sign-convention-free):
    @test abs(c_p1) == 1
    @test abs(c_m1) == 1
    @test c_p1 == -c_m1                  # the two regimes carry opposite Chern
    @test chern_number(qwz(3.0), 1) == 0   # trivial |m| > 2
    @test chern_number(qwz(-3.0), 1) == 0
    # FHS exactness: a different coarse grid gives the identical integer
    @test chern_number(qwz(1.0), 1; nk=16) == c_p1
    @test chern_number(qwz(1.0), 1; nk=37) == c_p1
    # gap closing (m = 2: Dirac point) is rejected
    @test_throws ErrorException chern_number(qwz(2.0), 1)
end

@testset "TKNN quantization relation" begin
    @test residual(TKNN(); σxy=1, C=1) == 0
    @test check(TKNN(); σxy=-1, C=-1)
    @test !check(TKNN(); σxy=0, C=1)
    @test solve(TKNN(), Val(:σxy); C=2) == 2
    @test solve(TKNN(), Val(:C); σxy=-1) == -1
end
