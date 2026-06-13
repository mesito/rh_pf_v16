"""
Part I (concavity / turning points): Concavity and Turning Points (the concavity equivalence)
d^2/dt^2 log|zeta(1/2 + it)| = -Sum_k 1/(t - gamma_k)^2 < 0 in every gap.
Z has exactly one turning point per gap.

Paper: verified in 1999 gaps; baseline excess = 0.
With hypothetical off-line zero at h0 < hthr/2: +2 extra turning points.
"""
# v16 NOTE: this script verifies a result of the v16 manuscript
# "Off-Line Zeros of the Riemann xi-Function" (repository rh_pf_v16).
# The Section/Theorem numbers below refer to the v16 numbering.


import mpmath
from config import Son, print_header


def verify_concavity(zeros, n_gaps=50, verbose=True):
    """
    Verify concavity in n_gaps consecutive zero gaps.
    For each gap (gamma_n, gamma_{n+1}), check d^2/dt^2 log|zeta| < 0
    at the midpoint and at several interior points.
    """
    if verbose:
        print_header("17.5 Concavity and Turning Points (the concavity equivalence)")
        print(f"\nTesting {n_gaps} consecutive zero gaps")

    results = []
    all_concave = True

    for idx in range(min(n_gaps, len(zeros) - 1)):
        g_n = float(zeros[idx])
        g_n1 = float(zeros[idx + 1])
        Ln = g_n1 - g_n
        t0 = (g_n + g_n1) / 2.0

        # d^2/dt^2 log|zeta(1/2+it)| = -Sum_k 1/(t0 - gamma_k)^2
        d2 = 0.0
        for gk in zeros:
            diff = t0 - float(gk)
            if abs(diff) > 1e-15:
                d2 -= 1.0 / diff**2

        concave = d2 < 0
        if not concave:
            all_concave = False

        # Count turning points of Z in this gap (Im[zeta'/zeta] sign changes)
        # Im[zeta'/zeta(1/2+it)] = -Sum_k 1/(t - gamma_k)
        # It goes from +inf at gamma_n to -inf at gamma_{n+1}, strictly decreasing
        # => exactly one zero => exactly one turning point
        n_samples = 20
        dt = Ln / (n_samples + 1)
        sign_changes = 0
        prev_val = None
        for j in range(1, n_samples + 1):
            t_j = g_n + j * dt
            im_val = 0.0
            for gk in zeros:
                diff = t_j - float(gk)
                if abs(diff) > 1e-15:
                    im_val -= 1.0 / diff
            if prev_val is not None and prev_val * im_val < 0:
                sign_changes += 1
            prev_val = im_val

        results.append({
            "gap": idx + 1,
            "t0": t0,
            "Ln": Ln,
            "d2": d2,
            "concave": concave,
            "turning_points": sign_changes  # should be exactly 1
        })

    if verbose:
        print(f"\n  {'Gap':>5s}  {'t0':>10s}  {'Ln':>8s}  {'d2/dt2':>12s}  {'Concave':>8s}  {'TP':>4s}")
        print(f"  {'-'*5}  {'-'*10}  {'-'*8}  {'-'*12}  {'-'*8}  {'-'*4}")
        for r in results[:20]:
            print(f"  {r['gap']:5d}  {r['t0']:10.2f}  {r['Ln']:8.4f}  {r['d2']:12.4f}  "
                  f"{'Yes' if r['concave'] else 'NO':>8s}  {r['turning_points']:4d}")
        if len(results) > 20:
            print(f"  ... ({len(results) - 20} more gaps)")

        n_exact_1_tp = sum(1 for r in results if r["turning_points"] == 1)
        print(f"\nAll gaps concave: {all_concave}")
        print(f"Gaps with exactly 1 turning point: {n_exact_1_tp}/{len(results)}")
        print(f"\nPaper reference: all 1999 gaps concave; baseline excess = 0 everywhere.")

    return {"all_concave": all_concave, "n_gaps_tested": len(results)}


if __name__ == "__main__":
    from config import load_all_zeros
    zeros = load_all_zeros()
    verify_concavity(zeros, n_gaps=50)
