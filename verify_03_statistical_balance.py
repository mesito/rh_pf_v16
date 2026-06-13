"""
Part I (statistical balance, heuristic): Statistical Balance (the statistical-balance heuristic)
Sigma-equivalents of the needed fluctuation at h0 = 1/log(T).

Paper Table 2:
  T=50:    174 sigma
  T=100:   536 sigma
  T=1000:  16067 sigma
  T=10^4:  355320 sigma
  T=10^6:  1.10e8 sigma
"""
# v16 NOTE: this script verifies a result of the v16 manuscript
# "Off-Line Zeros of the Riemann xi-Function" (repository rh_pf_v16).
# The Section/Theorem numbers below refer to the v16 numbering.


import mpmath
from config import A_CONST, print_header


def verify_statistical_balance(zeros, verbose=True):
    if verbose:
        print_header("17.3 Statistical Balance (the statistical-balance heuristic)")

    N0 = len(zeros)
    A = A_CONST

    # T values from Table 2
    test_cases = [50, 100, 1000, 10000, 1000000]

    if verbose:
        print(f"\nN0 = {N0}, A = {A}")
        print(f"h0 = 1/log(T) (sign barrier depth)")
        print(f"\n  {'T':>10s}  {'Needed |f|':>15s}  {'Typical sigma':>15s}  {'Sigma-equiv':>15s}")
        print(f"  {'-'*10}  {'-'*15}  {'-'*15}  {'-'*15}")

    results = []
    for T in test_cases:
        logT = mpmath.log(T)
        h0 = 1.0 / float(logT)

        # Needed |f| = A*N0 / h0^2
        # For the purpose of this table, N0 scales as c*T*logT
        # At given T: use N0(T) ~ (c * T * logT) / (2*pi), with c ~ 1
        # But the paper uses actual N0 from loaded zeros when T is small
        # For large T: extrapolate
        if T <= 2515:
            # Count zeros up to T
            n0_T = sum(1 for g in zeros if float(g) <= T)
        else:
            # Extrapolate: N0(T) ~ T/(2*pi) * log(T/(2*pi*e))
            n0_T = int(float(T / (2 * mpmath.pi) * mpmath.log(T / (2 * mpmath.pi * mpmath.e))))

        needed = A * n0_T / h0**2

        # Typical sigma from Selberg CLT: variance ~ (1/2) log log T for log|zeta|.
        # (v14 Table 2 convention: typical_sigma = sqrt(0.5 * log log T), no 1/h0 factor.)
        typical_sigma = mpmath.sqrt(0.5 * mpmath.log(logT))
        sigma_equiv = needed / float(typical_sigma)

        results.append((T, float(needed), float(typical_sigma), sigma_equiv))
        if verbose:
            if needed > 1e6:
                print(f"  {T:10d}  {float(needed):15.2e}  {float(typical_sigma):15.2f}  {sigma_equiv:15.2e}")
            else:
                print(f"  {T:10d}  {float(needed):15.0f}  {float(typical_sigma):15.2f}  {sigma_equiv:15.0f}")

    if verbose:
        print(f"\nPaper reference (Table 2):")
        print(f"  T=50:    needed=143, sigma=0.83, equiv=174 sigma")
        print(f"  T=1000:  needed=15794, sigma=0.98, equiv=16067 sigma")
        print(f"  T=10^6:  needed=1.26e8, sigma=1.15, equiv=1.10e8 sigma")
        print(f"\nNote: exact values depend on A and N0 counting convention.")

    return results


if __name__ == "__main__":
    from config import load_all_zeros
    zeros = load_all_zeros()
    verify_statistical_balance(zeros)
