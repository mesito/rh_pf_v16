"""
Part I (fundamental identity): Fundamental Identity (the fundamental identity)
V'(h,t) = A*N0/h^2 + Re[zeta'/zeta(1/2 + h + it)] > 0

Verified at grid of (h, t) points.
Paper reports V' > 0 at all 70 test points (35 original + 35 extended).
"""
# v16 NOTE: this script verifies a result of the v16 manuscript
# "Off-Line Zeros of the Riemann xi-Function" (repository rh_pf_v16).
# The Section/Theorem numbers below refer to the v16 numbering.


import mpmath
from config import A_CONST, print_header


def verify_fundamental_identity(zeros, verbose=True):
    if verbose:
        print_header("17.1 Fundamental Identity (the fundamental identity)")

    N0 = len(zeros)
    A = A_CONST

    h_values = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    # Original + extended t values
    t_values_orig = [17.5, 25, 35, 50, 70]
    t_values_ext = [500, 1000, 1500, 2000, 2500]
    t_values = t_values_orig + t_values_ext

    results = []
    all_positive = True
    smallest_V = None
    smallest_info = None

    for t in t_values:
        for h in h_values:
            s = mpmath.mpc(0.5 + h, t)
            zeta_prime_over_zeta = mpmath.re(mpmath.diff(mpmath.zeta, s) / mpmath.zeta(s))
            V_prime = A * N0 / h**2 + zeta_prime_over_zeta

            positive = float(V_prime) > 0
            if not positive:
                all_positive = False
            if smallest_V is None or float(V_prime) < float(smallest_V):
                smallest_V = V_prime
                smallest_info = (h, t)

            results.append((h, t, float(V_prime), positive))

    if verbose:
        print(f"\nN0 = {N0} zeros, A = {A}")
        print(f"Tested {len(results)} points ({len(h_values)} h-values x {len(t_values)} t-values)")
        print(f"\nSample results:")
        print(f"  {'h':>8s}  {'t':>8s}  {'V_prime':>15s}  {'V>0?':>5s}")
        print(f"  {'-'*8}  {'-'*8}  {'-'*15}  {'-'*5}")
        for h, t, vp, pos in results[:15]:
            print(f"  {h:8.3f}  {t:8.1f}  {vp:15.1f}  {'Yes' if pos else 'NO'}")
        if len(results) > 15:
            print(f"  ... ({len(results)-15} more rows)")

        print(f"\nSmallest V' = {float(smallest_V):.1f} at h={smallest_info[0]}, t={smallest_info[1]}")
        print(f"All positive: {all_positive}")
        print(f"\nPaper reference: smallest V'(0.5, 17.5) ~ 35.7; all 70 points positive.")

    return {"all_positive": all_positive, "smallest_V": float(smallest_V),
            "n_points": len(results)}


if __name__ == "__main__":
    from config import load_all_zeros
    zeros = load_all_zeros()
    verify_fundamental_identity(zeros)
