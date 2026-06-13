"""
verify_residue_krein.py -- Theorem 3.1 (residue split) and the Krein level set.

Checks:
  * Res_{rho} H'/H = 1, Res L1'/L1 = 0, Res I_D = 1, at rho1 and an on-line zero.
  * R(rho1) = L2/L1 = -c1/c2 to input precision; |R(rho1)| = 1.
  * generic interior |R| values are O(1) and != 1.
"""
import mpmath as mp
from ic_core import KAPPA, L_chi, H, RHO1, banner

mp.mp.dps = 25
I = mp.mpc(0, 1)
# H = c1 L1 + c2 L2 with L1 = L(s,chi), L2 = L(s,chibar), c1=(1-i k)/2, c2=(1+i k)/2
c1 = (1 - I*KAPPA)/2
c2 = (1 + I*KAPPA)/2
L1 = lambda s: L_chi(s, conj=False)
L2 = lambda s: L_chi(s, conj=True)
R  = lambda s: L2(s)/L1(s)

def residue(f, z0, r=mp.mpf('0.008'), n=64):
    """Contour residue (1/2pi i) oint f dz by n-node trapezoid on circle radius r.
       dz = i r e^{i th} d(th), d(th)=2pi/n; the 2pi i cancels, leaving sum/n."""
    tot = mp.mpc(0)
    for k in range(n):
        th = 2*mp.pi*k/n
        z = z0 + r*mp.e**(I*th)
        tot += f(z)*(I*r*mp.e**(I*th))          # f * dz/dth
    return tot/n

def run():
    banner("Theorem 3.1: residue split and Krein level set")
    dHH = lambda s: mp.diff(H, s)/H(s)
    dL1 = lambda s: mp.diff(L1, s)/L1(s)
    ID  = lambda s: dHH(s) - dL1(s)
    for name, z0 in [("rho1", RHO1), ("on-line gamma=87.64747", mp.mpf('0.5')+mp.mpf('87.64747')*I)]:
        rH  = residue(dHH, z0)
        rL1 = residue(dL1, z0)
        rID = residue(ID,  z0)
        print(f"  {name}:")
        print(f"    Res H'/H   = {mp.nstr(rH,6)}   (paper: 1)")
        print(f"    Res L1'/L1 = {mp.nstr(rL1,3)}   (paper: 0)")
        print(f"    Res I_D    = {mp.nstr(rID,6)}   (paper: 1)")
    # Krein level set at rho1
    target = -c1/c2
    val = R(RHO1)
    print(f"  Krein: R(rho1) = {mp.nstr(val,10)}")
    print(f"         -c1/c2  = {mp.nstr(target,10)}")
    print(f"         |R(rho1)-(-c1/c2)| = {mp.nstr(abs(val-target),3)}  (paper: 8e-11)")
    print(f"         |R(rho1)| = {mp.nstr(abs(val),10)}  (paper: 1.000000000)")
    # generic interior values
    pts = [mp.mpf('0.6')+mp.mpf('80')*I, mp.mpf('0.9')+mp.mpf('90')*I, mp.mpf('0.7')+mp.mpf('100')*I]
    gv = [mp.nstr(abs(R(p)),4) for p in pts]
    print(f"         generic |R| at non-zeros: {gv}  (paper: 0.187, 4.062, 1.266)")

if __name__ == "__main__":
    run()
