"""
Part II: Gram Tensor and Energy-Budget Ratio (Part II)

Two E/B bounds from the Gram/energy-budget bound:
  (i)  Speiser-free: E/B >= G_{II}/B = pi*log^3(T)/2  [valid at ALL depths]
  (ii) Full bound:   E/B >= c0*log^5(T), c0 = 120/(2pi)^5  [Theorem 17.2]

The log^5 bound is the headline result; it is dominated by G_{IV,IV} (Speiser class).
"""
# v16 NOTE: this script verifies a result of the v16 manuscript
# "Off-Line Zeros of the Riemann xi-Function" (repository rh_pf_v16).
# The Section/Theorem numbers below refer to the v16 numbering.


import math
from config import print_header


def gram_diagonal(logT):
    """Gram tensor diagonal elements at height T (using logT directly)."""
    two_pi = 2 * math.pi
    return {
        "I":   logT**2 / (4 * math.pi**2),
        "II":  logT**4 / 4,                       # class II: b=1/2, h0-independent
        "III": 6 * logT**4 / two_pi**4,
        "IV":  120 * logT**6 / two_pi**6,          # dominant term -> log^5 in E/B
        "V":   math.log(logT / two_pi) if logT > two_pi else 0.01
    }


def verify_energy_budget(verbose=True):
    if verbose:
        print_header("Part II: Energy-Budget Ratio (Part II)")

    log10_values = [12, 23, 100, 1000]

    # --- Both bounds ---
    if verbose:
        print("\nthe Gram/energy-budget bound: Two E/B bounds")
        print(f"\n  {'T':>10s}  {'G_II/B':>12s}  {'pi*log^3/2':>12s}  "
              f"{'E/B (full)':>12s}  {'c0*log^5':>12s}")
        print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")

    c0 = 120 / (2 * math.pi)**5
    ratios = []

    for log10T in log10_values:
        logT = log10T * math.log(10)
        G = gram_diagonal(logT)
        E = sum(G.values())
        B = logT / (2 * math.pi)

        ratio_II = G["II"] / B          # = pi*log^3/2
        sf_pred  = math.pi * logT**3 / 2

        ratio_full = E / B
        c0_pred    = c0 * logT**5

        ratios.append((log10T, ratio_full))

        if verbose:
            print(f"  10^{log10T:>4d}  {ratio_II:12.2e}  {sf_pred:12.2e}  "
                  f"{ratio_full:12.2e}  {c0_pred:12.2e}")

    if verbose:
        print(f"\n  c0 = 120/(2*pi)^5 = {c0:.6f}")
        print(f"\n  (i)  G_II/B  = pi*log^3(T)/2   -- valid at ALL depths incl. self-consistency")
        print(f"  (ii) E/B     >= c0*log^5(T)     -- Theorem 17.2 (outside SC regime)")
        print(f"  Both -> infinity as T -> infinity.")

    # --- Individual condition costs ---
    if verbose:
        print(f"\nIndividual condition costs at T = 10^12 (Section 25):")
        logT = 12 * math.log(10)
        G = gram_diagonal(logT)
        h0 = math.pi / logT
        log_h0 = abs(math.log(h0))
        B = logT / (2 * math.pi)

        for name, key, scaling in [
            ("(a) Apollonius", "I",   "log^2 T"),
            ("(b) Sign barrier", "II",  "log^4 T"),
            ("(c) SC / (e) Curv", "III", "log^4 T"),
            ("(h) Speiser",     "IV",  "log^6 T"),
            ("(g) Inner/Euler", "V",   "log log T"),
        ]:
            eps_k = G[key] / log_h0
            print(f"    {name:>20s} (class {key:>3s}):  eps = {eps_k:12.1f}   [{scaling}]")
        print(f"    {'Budget':>20s}:               B = {B:12.2f}   [log T]")
        print(f"\n    Speiser (h) alone exceeds budget by factor ~{G['IV']/log_h0/B:.0f}")

    # --- Lambda bound (the frozen-field bound / Corollary, v14) ---
    if verbose:
        print(f"\nFrozen-field bound (v16, full Trudgian):")
        import math
        # FULL Trudgian bound on S(T): 0.112 log t + 0.278 loglog t + 2.510
        # (truncating to the leading 0.111 log t understates the gap ~2x and is NOT legitimate)
        T0 = 3e12; Sbar = 0.112*math.log(T0)+0.278*math.log(math.log(T0))+2.510
        # extremal nearest-neighbour distance in z-coordinates
        Lmax = 2*math.pi/math.log(T0)*(1+2*Sbar)   # = 3.1323
        half_sq = (Lmax / 2.0)**2

        # The UNCONDITIONAL bound uses the critical-strip depth ceiling h0 < 1/2
        # (NOT the self-consistency depth 0.807), together with the rigorous
        # two-zero collision integral.  The full backward-flow collision time is
        # bounded above by this integral via the comparison principle.
        #   tau_2(h0=1/2, L_n=Lmax) = int_0^{1/2} h(h^2 + (L/2)^2)/((L/2)^2 + 5h^2) dh
        # Closed form (v16): tau_2(y0,a) = y0^2/10 + 2 a^2/25 ln(1 + 5 y0^2/a^2)
        # at y0=1 (strip ceiling, z-coords), a = x1* = Lmax.  Validated == 0.4233.
        y0 = 1.0
        tau2 = y0**2/10.0 + 2.0*Lmax**2/25.0 * math.log(1.0 + 5.0*y0**2/Lmax**2)
        eta2_at_Cn1 = 1.0/5.0 + 4.0/25.0*math.log(6)  # = 0.4867, ratio at C_n=1

        print(f"  Critical-strip depth ceiling (standard time): y0 = 2 h0 < 1")
        print(f"  Full Trudgian: Sbar(T0) = {Sbar:.4f}, budget 2*Sbar = {2*Sbar:.2f} (paper 13.32)")
        print(f"  Extremal nearest-neighbour distance: x1* = {Lmax:.4f} (paper 3.1323)")
        print(f"  Static two-zero ladder value: tau_2 = {tau2:.4f} (paper 0.4233 at y0=1)")
        print(f"  Note: the budget-extremal frozen-field bound is tau*(T0)=0.229")
        print(f"  (see verify_flow_bounds.py); the dynamic two-zero value is 0.4305.")
        print(f"  The dead value 0.081 of v14 came from a truncated Trudgian (L=1.614)")
        print(f"  and a 4x clock error; both are corrected in v16.")

    return {"ratios": ratios, "c0": c0}


if __name__ == "__main__":
    verify_energy_budget()
