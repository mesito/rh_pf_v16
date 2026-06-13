"""
Part I-II (collision-time ratio, Layer B): Collision Time -- FULL-ODE approximate ratio (Layer B).

LAYER B (numerical observation). This computes the FULL backward-flow
collision-time ratio eta = tau_ODE / (h0^2/2) using nearby on-line zeros.
The ratio is approximately stable, eta ~ 0.46 with per-gap CV ~ 2.7%
across the tested gaps -- an approximately stable ratio (NOT a universal constant).

IMPORTANT: This Layer-B full-ODE eta is DISTINCT from, and NOT used by, the
frozen-field bound tau*(T0)=0.229 (v16).  That bound (Layer A) uses the rigorous
TWO-ZERO collision integral (an absolutely convergent quadrature, see
verify_10_energy_budget.py) together with the critical-strip ceiling h0 < 1/2.
Do not conflate the two.

The collision time is computed via numerical quadrature:
  tau = integral_0^{h0} h / (1 + 2*h^2 * P_on(h, t0)) dh
where P_on(h, t0) = Sum_j 1/((t0 - gamma_j)^2 + h^2).

eta depends weakly on gap structure: corr(eta,s) ~ -0.96, corr(eta,Son) ~ +0.42, corr(eta,T) ~ 0.00.
"""
# v16 NOTE: this script verifies a result of the v16 manuscript
# "Off-Line Zeros of the Riemann xi-Function" (repository rh_pf_v16).
# The Section/Theorem numbers below refer to the v16 numbering.


import math
import statistics
from config import print_header


def collision_time_quadrature(t0, h0, nearby_gammas, n_points=2000):
    """
    Compute collision time via composite Simpson's rule:
      tau = integral from 0 to h0 of h / (1 + 2*h^2 * P_on(h)) dh
    """
    a, b = 1e-12, h0
    n = n_points if n_points % 2 == 0 else n_points + 1
    dx = (b - a) / n

    def integrand(h):
        P_on = sum(1.0 / ((t0 - gj)**2 + h**2) for gj in nearby_gammas)
        return h / (1.0 + 2.0 * h**2 * P_on)

    total = integrand(a) + integrand(b)
    for i in range(1, n):
        x = a + i * dx
        total += (4 if i % 2 else 2) * integrand(x)

    return total * dx / 3.0


def verify_collision_time(zeros, n_gaps=50, verbose=True):
    if verbose:
        print_header("17.9 Collision Time (Theorem 16.1)")
        print(f"\nComputing collision time via Simpson quadrature")

    n_gaps = min(n_gaps, len(zeros) - 1)
    n_nearby = min(200, len(zeros))

    eta_values = []
    son_values = []
    s_values = []
    results = []

    # Select gaps spread across range
    step = max(1, n_gaps // 25)
    gap_indices = sorted(set(
        list(range(0, min(10, n_gaps))) +
        list(range(0, n_gaps, step)) +
        [n_gaps - 1]
    ))
    gap_indices = [i for i in gap_indices if i < n_gaps]

    for idx in gap_indices:
        g_n = float(zeros[idx])
        g_n1 = float(zeros[idx + 1])
        Ln = g_n1 - g_n
        t0 = (g_n + g_n1) / 2.0

        son_val = sum(1.0 / (t0 - float(gk))**2 for gk in zeros if abs(t0 - float(gk)) > 1e-12)
        h0 = math.sqrt(2.0 / son_val)

        # Select nearby zeros
        nearby = sorted(zeros, key=lambda g: abs(float(g) - t0))
        nearby_f = [float(g) for g in nearby[:n_nearby]]

        tau = collision_time_quadrature(t0, h0, nearby_f, n_points=2000)

        h0_sq_half = h0**2 / 2.0
        eta = tau / h0_sq_half

        logT = math.log(max(t0, 2.0))
        L_bar = 2 * math.pi / logT
        s_n = Ln / L_bar

        eta_values.append(eta)
        son_values.append(son_val)
        s_values.append(s_n)

        results.append({
            "gap": idx + 1, "t0": t0, "Son": son_val,
            "h0": h0, "s": s_n, "tau": tau,
            "h0sq2": h0_sq_half, "eta": eta,
        })

    if verbose:
        print(f"\n  {'Gap':>5s}  {'t0':>8s}  {'Son':>8s}  {'h0':>8s}  {'s':>6s}  "
              f"{'tau':>10s}  {'h0^2/2':>10s}  {'eta':>8s}")
        print(f"  {'-'*5}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*6}  "
              f"{'-'*10}  {'-'*10}  {'-'*8}")

        for r in results:
            print(f"  {r['gap']:5d}  {r['t0']:8.1f}  {r['Son']:8.4f}  {r['h0']:8.3f}  {r['s']:6.2f}  "
                  f"{r['tau']:10.4f}  {r['h0sq2']:10.4f}  {r['eta']:8.4f}")

        if len(eta_values) > 1:
            eta_mean = statistics.mean(eta_values)
            eta_std = statistics.stdev(eta_values)
            eta_cv = eta_std / eta_mean * 100

            print(f"\neta = {eta_mean:.4f} +/- {eta_std:.4f} (CV = {eta_cv:.2f}%)")

            # Correlations
            if len(eta_values) >= 5:
                def corr(x, y):
                    n = len(x)
                    mx, my = statistics.mean(x), statistics.mean(y)
                    cov = sum((a-mx)*(b-my) for a,b in zip(x,y))/(n-1)
                    sx, sy = statistics.stdev(x), statistics.stdev(y)
                    return cov/(sx*sy) if sx>0 and sy>0 else 0

                t_vals = [r['t0'] for r in results]
                print(f"Correlation eta vs Son: r = {corr(eta_values, son_values):+.3f}")
                print(f"Correlation eta vs s:   r = {corr(eta_values, s_values):+.3f}")
                print(f"Correlation eta vs T:   r = {corr(eta_values, t_vals):+.3f}")

        print(f"\nPaper reference: eta ~ 0.46, CV ~ 2.7% (approximately stable, Layer B)")
        print(f"  (weak gap-structure dependence: corr(eta,s) ~ -0.96, corr(eta,Son) ~ +0.42)")

    return {
        "eta_mean": statistics.mean(eta_values) if eta_values else 0,
        "eta_std": statistics.stdev(eta_values) if len(eta_values) > 1 else 0,
        "n_gaps": len(results)
    }


if __name__ == "__main__":
    from config import load_all_zeros
    zeros = load_all_zeros()
    verify_collision_time(zeros, n_gaps=50)
