"""
verify_flow_bounds.py -- Part I/II of the v16 manuscript: the exactly solvable
collision model, the frozen-field budget-extremal bound, the verification
dividend, and the two-zero ladder. All in standard de Bruijn time.

Reproduces:
  invariant I = y^4 + 10 Q^2 y^2 + Q^4 conserved; tau_dyn = 0.4305;
  static two-zero tau_2 = 0.4233; eta_2 = 1/5 + 4/25 ln6 = 0.4867;
  x1*(T0) = 3.1323; budget B_S(T0) = 13.32;
  frozen-field tau*(T0) = 0.229; dividend crossing tau*=0.2 at T ~ 1.25e15;
  scaling law tau* -> 9.4/log T.
"""
import math
import numpy as np
from ic_core import banner

def Sbar(t):
    """Full Trudgian bound on S(T): 0.112 log t + 0.278 loglog t + 2.510."""
    return 0.112*math.log(t) + 0.278*math.log(math.log(t)) + 2.510

def tau_static(y0, a):
    """Static two-zero closed form: y0^2/10 + 2a^2/25 ln(1 + 5 y0^2/a^2)."""
    return y0**2/10 + 2*a**2/25*math.log(1 + 5*y0**2/a**2)

def tau_dyn(y0, Q0):
    I0 = y0**4 + 10*Q0**2*y0**2 + Q0**4
    return (math.sqrt(I0) - (Q0**2 - y0**2))/12

def tau_star(logT, Nz=20000, Ny=1500):
    """Frozen-field budget-extremal collision time at verification height e^logT."""
    Sb = Sbar(math.exp(logT))
    a = 4*math.pi/logT
    x1 = (1 + 2*Sb)*a/2
    xj = x1 + np.arange(Nz)*a
    xN = xj[-1] + a
    ys = np.linspace(1e-6, 1.0, Ny)
    S = np.array([np.sum(4*y/(xj*xj + y*y)) + (4/a)*(math.pi/2 - math.atan(xN/y)) for y in ys])
    return float(np.trapezoid(1.0/(1.0/ys + S), ys))

def run():
    banner("Part I: exactly solvable collision model")
    y0, Q0 = 1.0, 3.1323
    A0, B0 = y0**2 - Q0**2, -(y0**2)*(Q0**2)
    I0 = A0**2 - 12*B0
    # invariant conserved along the coefficient flow A(t)=A0-12t, B(t)=B0-2A0 t+12 t^2
    drift = max(abs(((A0-12*t)**2 - 12*(B0 - 2*A0*t + 12*t*t)) - I0)
                for t in np.linspace(0, tau_dyn(y0, Q0), 300))
    print(f"  invariant I = y^4+10Q^2y^2+Q^4 = {y0**4+10*Q0**2*y0**2+Q0**4:.4f} = A0^2-12B0 = {I0:.4f}")
    print(f"  conservation drift along flow: {drift:.2e}")
    print(f"  tau_dyn (dynamic two-zero, exact) = {tau_dyn(1.0,3.1323):.4f}  (paper 0.4305)")
    print(f"  tau_2  (static two-zero)          = {tau_static(1.0,3.1323):.4f}  (paper 0.4233)")
    print(f"  eta_2 = 1/5 + 4/25 ln6            = {1/5 + 4/25*math.log(6):.4f}  (paper 0.4867)")
    print(f"  pure-conjugate limit  Q->inf: tau -> y0^2/2 = 0.5 (de Bruijn anchor)")

    banner("Part II: frozen-field bound and the verification dividend")
    T0 = 3e12; logT0 = math.log(T0)
    x1 = 2*math.pi/logT0*(1 + 2*Sbar(T0))
    print(f"  x1*(T0) = 2pi/logT0 (1+2 Sbar) = {x1:.4f}  (paper 3.1323)")
    print(f"  budget B_S(T0) = 2 Sbar(T0)    = {2*Sbar(T0):.2f}   (paper 13.32)")
    print(f"  tau*(T0) frozen field          = {tau_star(logT0):.4f}  (paper 0.229)")
    print("  verification dividend curve:")
    print(f"    {'logT':>6} {'T':>11} {'Sbar':>7} {'tau*':>8}")
    for lT, Ts in [(28.73,'3.0e12'),(34.76,'1.25e15'),(40,'2.4e17'),(50,'5.2e21'),(120,'1.3e52')]:
        print(f"    {lT:>6.2f} {Ts:>11} {Sbar(math.exp(lT)):>7.2f} {tau_star(lT):>8.4f}")
    # crossing tau*=0.2
    lo, hi = 28.73, 80
    for _ in range(24):
        m = (lo+hi)/2
        if tau_star(m, 15000, 1000) > 0.2: lo = m
        else: hi = m
    print(f"  crossing tau*=0.2 at logT={(lo+hi)/2:.2f}, T={math.exp((lo+hi)/2):.2e}  (paper 1.25e15)")
    print(f"  scaling law tau* * logT: " + ", ".join(f"{tau_star(lT)*lT:.2f}" for lT in [80,120,200]) + " -> 9.4")

if __name__ == "__main__":
    run()
