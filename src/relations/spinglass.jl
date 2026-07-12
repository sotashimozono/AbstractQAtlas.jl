# relations/spinglass.jl — disordered / spin-glass identities.
#
# The exact, model-independent facts of a quenched-disordered magnet: the
# Edwards–Anderson order parameter, the exact Nishimori-line identities
# (gauge symmetry of the ±J model), and the de Almeida–Thouless
# replica-symmetry-breaking boundary.  Domain tag :spinglass.
#
# References (doiget-verified): Edwards & Anderson, J. Phys. F 5, 965
# (1975); Nishimori, Prog. Theor. Phys. 66, 1169 (1981); de Almeida &
# Thouless, J. Phys. A 11, 983 (1978).

"""
    EdwardsAndersonOrderParameter <: AbstractRelation

The Edwards–Anderson order parameter as the replica self-overlap (Edwards
& Anderson, J. Phys. F 5, 965 (1975)),

`q_EA = [⟨s_i⟩²]`

(the disorder average of the squared thermal magnetization =
`(1/N) Σ_i [⟨s_i⟩²]`, the caller-supplied `overlap`).

Variables: `q_EA`, `overlap` = `[⟨s_i⟩²]`.
"""
@relation :spinglass EdwardsAndersonOrderParameter(q_EA, overlap) = q_EA - overlap

"""
    NishimoriEnergy <: AbstractRelation

The exact internal energy per bond of the ±J Ising model on the Nishimori
line (Nishimori, Prog. Theor. Phys. 66, 1169 (1981)),

`U = −J tanh(βJ)`,

a gauge-symmetry consequence that holds for ANY lattice/dimension — a
rare closed-form result in a disordered system, and a hard check for a
disorder-averaged simulation on the Nishimori line.

Variables: `U`, `J`, `β` (or `T`).
"""
@relation :spinglass NishimoriEnergy(U, J, β) = U + J * tanh(β * J)

"""
    NishimoriMagnetizationOverlap <: AbstractRelation

The Nishimori gauge identity: on the Nishimori line the spin-glass order
parameter equals the ferromagnetic magnetization (Nishimori, Prog. Theor.
Phys. 66, 1169 (1981)),

`q = m`,

so the spin-glass and ferromagnetic order coincide — there is no
replica-symmetry-broken spin-glass phase below the Nishimori line.

Variables: `q`, `m`.
"""
@relation :spinglass NishimoriMagnetizationOverlap(q, m) = q - m

"""
    AlmeidaThoulessStability <: AbstractInequality

The de Almeida–Thouless replica-symmetric stability criterion (de Almeida
& Thouless, J. Phys. A 11, 983 (1978)): the replicon eigenvalue is
non-negative in the replica-symmetric phase,

`1 − (βJ)² [sech⁴(βh)] ≥ 0`

(slack = the replicon eigenvalue; `= 0` on the AT line, `< 0` in the
replica-symmetry-broken phase).  `sech4_avg = [sech⁴(βh)]` is the
disorder-averaged local-field factor.

Variables: `βJ`, `sech4_avg` = `[sech⁴(βh)]`.
"""
@inequality :spinglass AlmeidaThoulessStability(βJ, sech4_avg) = 1 - βJ^2 * sech4_avg
