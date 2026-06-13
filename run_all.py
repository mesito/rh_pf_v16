"""
run_all.py -- run every verification for the v16 manuscript.

Two families of scripts:
  * verify_01..10  -- the Part I/II machinery (fundamental identity, curvature
    detector, statistical balance, inner function, concavity, tunneling,
    self-consistency, Speiser, collision-time ratio eta, Gram/energy-budget),
    updated to v16 numbering. These reproduce the quantities referenced in the
    Code-availability statement: eta ~ 0.46 (CV 2.7%), C_n in [1.002,2.306],
    S_on, the Gram tensor G_ij, and E(T)/B(T).
  * verify_<name>  -- the new v16 verifiers for the exact closed-form bounds,
    the lifetime-deficit dictionary (Part III), and the interference channel
    / DH negative control (Part V).

Usage:  python run_all.py
Heavy scripts (residue, Weil triangle, deficit landscape) take a few minutes
each at 25-digit precision; run them individually if needed.
"""
import time, importlib, traceback

LEGACY = [
    ("verify_01_fundamental_identity", "Part I: fundamental identity V'(h,t)>0"),
    ("verify_02_curvature_test",       "Part I: curvature detector"),
    ("verify_03_statistical_balance",  "Part I: statistical balance (heuristic)"),
    ("verify_04_inner_function",       "Part I: inner/Blaschke |Phi|=1"),
    ("verify_05_concavity",            "Part I: concavity / turning points"),
    ("verify_06_tunneling",            "Part I: tunneling invariant"),
    ("verify_07_self_consistency",     "Part I: self-consistent-depth invariant"),
    ("verify_08_speiser",              "Part I: Speiser depth scaling"),
    ("verify_09_collision",            "Part I-II: collision-time ratio eta (Layer B)"),
    ("verify_10_energy_budget",        "Part II: Gram tensor and E(T)/B(T)"),
]
NEW = [
    ("verify_constraints_partI",  "Part I: constraint-network identities"),
    ("verify_flow_bounds",        "Part I-II: collision model, frozen field, dividend"),
    ("verify_gram_budget_partII", "Part II: supply envelope, positivity sliver"),
    ("verify_lifetime_deficit",   "Part III: lifetime-deficit dictionary"),
    ("verify_core_lemma23",       "Part V: Euler product / interference witness"),
    ("verify_residue_krein",      "Part V: residue split, Krein level set"),
    ("verify_weil_triangle",      "Part V: Weil witnesses, capacity triangle"),
    ("verify_flow_invariant",     "Part V: closed-form lifetimes, index path"),
    ("verify_qw_index",           "Part V: index-staircase cross-checks"),
]

# entry-point function names for the legacy scripts
LEGACY_FN = {
    "verify_01_fundamental_identity": "verify_fundamental_identity",
    "verify_02_curvature_test": "verify_curvature_test",
    "verify_03_statistical_balance": "verify_statistical_balance",
    "verify_04_inner_function": "verify_inner_function",
    "verify_05_concavity": "verify_concavity",
    "verify_06_tunneling": "verify_tunneling",
    "verify_07_self_consistency": "verify_self_consistency",
    "verify_08_speiser": "verify_speiser",
    "verify_09_collision": "verify_collision_time",
    "verify_10_energy_budget": "verify_energy_budget",
}

def run_group(group, header, zeros=None):
    print("\n" + "#"*72 + f"\n# {header}\n" + "#"*72)
    for mod, desc in group:
        print(f"\n----- {mod}  ({desc}) -----")
        try:
            m = importlib.import_module(mod)
            if mod in LEGACY_FN and hasattr(m, LEGACY_FN[mod]):
                getattr(m, LEGACY_FN[mod])(zeros)
            elif hasattr(m, "run"):
                m.run()
            elif hasattr(m, "main"):
                m.main()
            else:
                print("  (no entry point; import-only)")
        except Exception:
            print("  ERROR:\n" + traceback.format_exc())

def main():
    t0 = time.time()
    print("="*72)
    print("v16 manuscript -- full numerical verification suite")
    print("Off-Line Zeros of the Riemann xi-Function")
    print("repository: https://github.com/mesito/rh_pf_v16")
    print("="*72)
    from config import load_all_zeros
    print("\nLoading 2000 zeros for the Part I/II machinery (cached after first run)...")
    zeros = load_all_zeros(verbose=False)
    run_group(LEGACY, "PART I/II MACHINERY (verify_01..10, v16 numbering)", zeros)
    run_group(NEW, "NEW v16 VERIFIERS (closed forms, Parts III & V)")
    print(f"\n{'='*72}\n# complete in {time.time()-t0:.0f}s\n{'='*72}")

if __name__ == "__main__":
    main()
