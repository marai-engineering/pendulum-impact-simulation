import numpy as np
import matplotlib.pyplot as plt

from src.impact_solver import solve_case

base = {
    "g": 9.81,
    "m1": 1.0, "L1": 0.7,
    "m2": 0.7, "L2": 0.35,
    "m0": 0.2,
    "phi0_deg": 40,
    "phiB_deg": 40,
    "k_t": 300,
    "w2_before": 0.0,
}

# print one reference case
case0 = dict(base)
case0["e"] = 0.6
res0 = solve_case(case0)

print("\n--- Example case results (e = 0.6) ---")
for k, v in res0.items():
    print(f"{k} = {v:.4f}")

# parameter sweep
e_vals = np.linspace(0.0, 1.0, 51)
phiA_vals = []
w1_after_vals = []

for e in e_vals:
    case = dict(base)
    case["e"] = float(e)
    res = solve_case(case)
    phiA_vals.append(res["phiA_deg"])
    w1_after_vals.append(res["w1_after"])

plt.figure()
plt.plot(e_vals, phiA_vals)
plt.xlabel("Restitution coefficient e [-]")
plt.ylabel("Rebound angle phiA [deg]")
plt.title("Effect of restitution on rebound angle")
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/e_sweep_phiA.png", dpi=200)

plt.figure()
plt.plot(e_vals, w1_after_vals)
plt.xlabel("Restitution coefficient e [-]")
plt.ylabel("omega1_after [rad/s]")
plt.title("Effect of restitution on omega1_after")
plt.grid(True)
plt.tight_layout()
plt.savefig("figures/e_sweep_omega1_after.png", dpi=200)

plt.show()