"""
verify_flow_invariant.py -- Lemma 8.1 (closed-form lifetimes) and the
Morse index path of DH under the de Bruijn flow.

The backward heat flow on the quartic z^4 + A z^2 + B is linear on coefficients:
    Adot = -12,  Bdot = -2A,  so  A^2 - 12B = y^4 + 10 Q^2 y^2 + Q^4  is conserved,
    collision time  tau = (sqrt(y^4+10Q^2 y^2+Q^4) - (Q^2-y^2))/12.
Validated against direct ODE integration of the four-body system.
"""
import math
import numpy as np
from ic_core import banner

def tau_closed(y0, Q0):
    I = y0**4 + 10*Q0**2*y0**2 + Q0**4
    return (math.sqrt(I) - (Q0**2 - y0**2))/12

def tau_ode(y0, Q0):
    """Direct integration of the symmetric four-body flow dot z = 2 sum 1/(z-z')."""
    y, Q = y0, Q0
    t, dt = 0.0, 1e-6
    while y > 1e-7:
        dy = -1.0/y - 4*y/(y*y+Q*Q)
        dQ =  1.0/Q + 4*Q/(Q*Q+y*y)
        y += dy*dt; Q += dQ*dt; t += dt
        if t > 10: break
    return t

def run():
    banner("Lemma 8.1: closed-form lifetimes and the DH index path")
    # invariant conservation: the flow is LINEAR on quartic coefficients,
    # A(t)=A0-12t, B(t)=B0-2 int A = B0-2A0 t+6t^2, so A^2-12B is exactly conserved.
    y0, Q0 = 1.0, 3.1323
    A0, B0 = y0**2 - Q0**2, -(y0**2)*(Q0**2)
    I0 = A0**2 - 12*B0
    drift = 0.0
    for t in np.linspace(0, tau_closed(y0, Q0), 200):
        A = A0 - 12*t; B = B0 - 2*A0*t + 12*t*t   # B' = -2A = -2(A0-12t)
        drift = max(drift, abs((A*A - 12*B) - I0))
    print(f"  invariant A^2-12B along coefficient flow (exact): drift = {drift:.2e}")
    print(f"  I0 = y^4+10Q^2y^2+Q^4 = {y0**4+10*Q0**2*y0**2+Q0**4:.4f} = A0^2-12B0 = {I0:.4f}")
    print(f"  closed-form tau = {tau_closed(1.0,3.1323):.4f}, ODE tau = {tau_ode(1.0,3.1323):.4f} (dynamic two-zero; v16 0.4305)")
    print(f"  linear law d(Q^2-y^2)/dt = 12: check via I(t=0)={I0:.2f}, q_T^2=sqrt(I0)={math.sqrt(I0):.4f}")

    # DH lifetimes of the two negative planes -> index path
    print("  DH negative-plane lifetimes (two-zero upper bounds):")
    # depths y=2h; on-line neighbor Q from the wide-gap half-width in z-coords
    cases = [("rho2 (h2=0.1508)", 2*0.150830, 114.16),
             ("rho1 (h1=0.3085)", 2*0.308517, 85.70)]
    # neighbour distance Q: use the measured wide-gap half-width *2 (z=2gamma)
    # gap widths 4.345 (rho2) and 4.539 (rho1); half-width in gamma -> z: *2
    Qs = {"rho2 (h2=0.1508)": 2*(4.345/2), "rho1 (h1=0.3085)": 2*(4.539/2)}
    taus = {}
    for nm, y, t0 in cases:
        Q = Qs[nm]
        taus[nm] = tau_closed(y, Q)
        print(f"    {nm}: y={y:.4f}, Q={Q:.3f} -> tau = {taus[nm]:.4f}")
    print(f"  paper lifetimes: tau2=0.0449 (rho2), tau1=0.1819 (rho1)")
    print(f"  index path: ind_-(Q_t^DH) : 4 --[t~{taus['rho2 (h2=0.1508)']:.3f}]--> 2 --[t~{taus['rho1 (h1=0.3085)']:.3f}]--> 0")
    print(f"  (paper: 4 ->[0.045]-> 2 ->[0.182]-> 0)")

if __name__ == "__main__":
    run()
