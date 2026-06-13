# Off-Line Zeros of the Riemann ξ-Function — v16 verification suite

Numerical verification code for the manuscript

> **Off-Line Zeros of the Riemann ξ-Function: A Constraint Network, an Exactly Solvable Collision Model, Explicit Frozen-Field Bounds, a Lifetime–Deficit Dictionary for Weil Positivity, and an Interference-Channel Negative Control**
> Mesut Ismail, Department of Programming and Computer Technologies, Technical University of Sofia.
> ORCID: [0009-0001-0496-964X](https://orcid.org/0009-0001-0496-964X)

This single package reproduces every quantitative claim of the five-part manuscript. Each script prints its computed values next to the paper values.

## Quick start

```bash
pip install mpmath numpy
python run_all.py
```

The band-wide deficit landscape (Part III) and the DH on-line zero regeneration are the slow steps (a few minutes total at 25-digit precision); the cached zeros in `dh_online_true.pkl` make repeat runs fast for the other scripts.

## Map: manuscript part → script → reproduced numbers

| Part | Script | Key reproduced numbers |
|---|---|---|
| I (constraints) | `verify_constraints_partI.py` | self-consistent depth `h_thr=√(2/S_on)`; Blaschke modulus `\|Φ_off\|=1` on the critical line (phase-blind); Apollonius radius `R=h²/(2AN₀)→0` |
| I–II (flow) | `verify_flow_bounds.py` | invariant `I=195.37` conserved to `6e-14`; `tau_dyn=0.4305`; `tau_2=0.4233`; `eta_2=0.4867`; `x1*(T0)=3.1323`; budget `13.32`; frozen field `tau*(T0)=0.229`; dividend crossing `tau*=0.2` at `T=1.25e15`; scaling `tau*·logT → 9.4` |
| II (tensor) | `verify_gram_budget_partII.py` | Gram entries `G_kk` (unconditional prime sums), Speiser class IV dominant; `E/B→∞` structural divergence (with the currency-mismatch note) |
| III (deficit) | `verify_lifetime_deficit.py` | aggregate descent `Σ y_k² ↓ ≤ -2m` (shielding-robust); both DH witnesses inside `[δ/3e^{1/4w²}, δ]` windows; band-wide triangle closure; margin minima at widest-gap centres; clipping `+0.077` (free) vs `+0.0018` (measured), skewness `-0.50` |
| V (interference) | `verify_core_lemma23.py` | `κ=0.284079`; `c6=(1+κ²)log6=1.9364`; envelope `30.81 @ n=2856`; FE & reality |
| V | `verify_residue_krein.py` | `Res H'/H=1`, `Res L1'/L1=0`; Krein `\|R(ρ1)-(-c1/c2)\|=7.7e-11`, `\|R(ρ1)\|=1.000000000` |
| V | `verify_weil_triangle.py` | witnesses `-0.3451/-0.2786/-0.2295/-0.1159` (ρ1), `-0.0548/-0.0509/-0.0480` (ρ2); triangle `ζ +0.0226/+0.0461/+0.0230`, `DH -0.2786/+0.0744/+0.3457`; support split `-0.0094/+0.3550`, composite fraction `1.03` |
| V | `verify_flow_invariant.py` | invariant `A²-12B` conserved to `3e-14`; `tau=0.4305`; DH lifetimes `tau2=0.0451`, `tau1=0.1839`; staircase `4→2→0` |
| V | `verify_qw_index.py` | witness depth ratio `5.59` vs CCM eigenvalue ratio `5.37` (paper 5.5/5.4); staircase reference table |

## Methodological notes (read before trusting a number)

**DH on-line zeros and the wide-gap claim.** Both off-line witnesses sit in on-line gaps of length **4.54** and **4.35** against a mean spacing **≈ 1.7**. This is reproduced only with correct zero detection: a naive sign-change scan of `Re H(½+it)` overcounts (`Re H` has near-zero touches that are not zeros), inflating the count to ≈ 44 in `[80,120]` (mean spacing 0.84) and shrinking the apparent gaps. Confirming each candidate with `|H| < 1e-3` rejects these and gives **23** true zeros, mean spacing **1.692**, gaps **4.539 / 4.345** — the paper values to three digits. `dh_online_zeros()` does the confirmed detection; the cache covers `[72,127]`.

**The capacity-triangle archimedean factor of 2.** The archimedean side carries a factor of 2 (both sides of the critical line, matching the `2 ∂_t θ / 2π` of the Weil form). With it the triangle closes to `0.0004` (ζ) and `0.0074` (DH).

**Scope of the QW index verification.** The full QW matrix in the Connes–Consani–Moscovici eigenbasis (the `0→2→4` staircase table) was produced by the author's CCM-basis spectral code, using the exact inner product `<V_a|T(n)|V_b>` and the `sin²(tL/2)/L` normalization. `verify_qw_index.py` does **not** reconstruct that matrix from quadrature; it reproduces the scalar invariant (the depth ratio of the two negative planes, `5.59`) that cross-checks the published eigenvalue ratio (`5.37`) independently of the matrix assembly.

**A recorded correction.** An early implementation of the archimedean partial-fraction block carried a sign error in the denominator (`ω_m-ω_k` instead of `ω_k-ω_m`), shifting the low-frequency spectrum by `O(1)`; the corrected identity yields the ζ control at `-1e-3`.

## Precision

All spectral quantities use `mpmath` at 25 digits; frequency sweeps use `numpy`. The completed DH function's functional equation and reality hold to `~1e-25` in this implementation. Matrix entries carry quadrature errors of order `3e-2` against signals of order `0.4–2.3`; the two-zero lifetimes are frozen-field upper bounds.

## Files

**Part I/II machinery** (`verify_01`–`verify_10`, updated to v16 numbering) — these reproduce the quantities named in the manuscript's Code-availability statement: the collision-time ratio eta ~ 0.46 (CV 2.7%), the structural constants C_n in [1.002, 2.306] and S_on, the Gram metric tensor G_ij, and the energy budget E(T)/B(T).

```
verify_01_fundamental_identity.py   V'(h,t) = A N0/h^2 + Re[zeta'/zeta] > 0
verify_02_curvature_test.py         curvature detector (residual after subtracting on-line zeros)
verify_03_statistical_balance.py    sigma-equivalents (heuristic, Layer B)
verify_04_inner_function.py         inner/Blaschke |Phi|=1 on the critical line
verify_05_concavity.py              concavity of log|zeta| in 1999 gaps; turning points
verify_06_tunneling.py              tunneling invariant M h0 = pi s/sqrt(C_n)
verify_07_self_consistency.py       self-consistent-depth invariant h0^2 S_on = 2
verify_08_speiser.py                Speiser depth scaling; detector switch-off
verify_09_collision.py              full-ODE collision-time ratio eta ~ 0.46 (Layer B)
verify_10_energy_budget.py          Gram tensor G_ij and E(T)/B(T); two-zero ladder 0.4233
```

**New v16 verifiers** — exact closed forms, the lifetime-deficit dictionary (Part III), and the interference channel / DH negative control (Part V).

```
ic_core.py                    core: kappa, coefficients c_n, L(s,chi), H, Xi_DH, the two zeros
config.py                     zeros loader (2000 zeros, cached), S_on, constants
verify_constraints_partI.py   Part I:  constraint-network identities and detectors
verify_flow_bounds.py         Part I-II: collision model, frozen field tau*=0.229, dividend
verify_gram_budget_partII.py  Part II: supply envelope (cap 0.977, util 2.4%, sliver 0.453)
verify_lifetime_deficit.py    Part III: equivalence, aggregate descent, deficit landscape
verify_core_lemma23.py        Part V: Theorem A, Lemma 2.3, c6=1.9364, envelope, FE/reality
verify_residue_krein.py       Part V: Theorem 3.1 (residue split, Krein level set)
verify_weil_triangle.py       Part V: Weil witnesses, capacity triangle, support split
verify_flow_invariant.py      Part V: closed-form lifetimes, index path 4->2->0
verify_qw_index.py            Part V: witness-ratio cross-check, staircase reference
run_all.py                    runs both families, group by group
dh_online_true.pkl            cached true DH on-line zeros in [72,127]
Ismail_v16.tex / .pdf         the manuscript
```

All Section/Theorem numbers in the legacy scripts have been updated to the v16 numbering; the dead value 0.081 (a truncated-Trudgian artefact) and the Caster effect (retired in the audit) have been removed, and eta is described as an approximately stable ratio (CV 2.7%), not a universal constant.

## AI assistance

The author declares that the research programme, the mathematical ideas, the innovations, and the overall approach of this work are his own. During the preparation of the manuscript the author used Claude (Anthropic) as an AI assistant: for surveying existing literature and related concepts; for drafting, structuring, and editing the exposition; and in the development of this Python verification code. All mathematical content was verified by the author, both directly and through this numerical verification suite (2000 nontrivial zeros of ζ at 25-digit precision, mpmath), with quantitative claims cross-checked by independent methods. The AI tool is not an author; the author reviewed, edited, and verified all content and bears full and sole responsibility for the correctness and integrity of this work.

## License

MIT.
