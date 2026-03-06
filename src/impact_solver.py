import numpy as np

def solve_case(data):
    # ---------- inputs ----------
    g = data["g"]
    m1 = data["m1"]; L1 = data["L1"]
    m2 = data["m2"]; L2 = data["L2"]
    m0 = data["m0"]
    phi0_deg = data["phi0_deg"]
    phiB_deg = data["phiB_deg"]
    e = data["e"]
    k_t = data["k_t"]
    w2_before = data["w2_before"]

    # ---------- angles ----------
    phi0 = np.deg2rad(phi0_deg)
    phiB = np.deg2rad(phiB_deg)

    # ---------- inertias ----------
    theta_O = (1/3) * m1 * L1**2 + m0 * L1**2
    theta_Q = (1/3) * m2 * L2**2 + m0 * L2**2

    # ---------- energy before impact (your convention) ----------
    C1 = g * (m1*(L1/2) + m0*L1)
    U0 = -C1 * np.sin(phi0)
    UB = -C1 * np.cos(phiB)

    T_B = U0 - UB
    if T_B <= 0:
        raise ValueError(f"T_B must be > 0. Got {T_B:.6f}. Check U0/UB convention.")

    w1_before = np.sqrt((2 * T_B) / theta_O)

    # ---------- impact (cos projection, consistent) ----------
    r1 = L1 * np.cos(phiB)
    r2 = L2

    u1 = w1_before * r1
    u2 = w2_before * r2

    mred1 = theta_O / (r1**2)
    mred2 = theta_Q / (r2**2)

    v1_after = (mred1*u1 + mred2*u2 - mred2*e*(u1 - u2)) / (mred1 + mred2)
    v2_after = (mred1*u1 + mred2*u2 + mred1*e*(u1 - u2)) / (mred1 + mred2)

    w1_after = v1_after / r1
    w2_after = v2_after / r2

    # ---------- swing-back angle ----------
    T_after = 0.5 * theta_O * w1_after**2
    U2 = -C1 * np.cos(phiB)
    U3 = T_after + U2

    cos_phiA = -U3 / C1
    cos_phiA = np.clip(cos_phiA, -1.0, 1.0)
    phiA_deg = np.rad2deg(np.arccos(cos_phiA))

    # ---------- bar 2 natural frequency ----------
    m_tot_2 = m2 + m0
    COM_2 = (m2*(L2/2) + m0*L2) / m_tot_2
    k_eff = k_t - m_tot_2 * g * COM_2

    if k_eff <= 0:
        raise ValueError(f"k_eff must be > 0. Got {k_eff:.6f}.")

    wn = np.sqrt(k_eff / theta_Q)

    results = {
        "theta_O": theta_O,
        "theta_Q": theta_Q,
        "U0": U0,
        "UB": UB,
        "T_B": T_B,
        "w1_before": w1_before,
        "w1_after": w1_after,
        "w2_after": w2_after,
        "phiA_deg": phiA_deg,
        "COM_2": COM_2,
        "k_eff": k_eff,
        "wn": wn,
    }
    return results

if __name__ == "__main__":

    case = {
        "g": 9.81,
        "m1": 1.0,
        "L1": 0.7,
        "m2": 0.7,
        "L2": 0.35,
        "m0": 0.2,
        "phi0_deg": 40,
        "phiB_deg": 40,
        "e": 0.6,
        "k_t": 300,
        "w2_before": 0.0,
    }

    res = solve_case(case)

    print("\n--- Results ---")
    for k, v in res.items():
        print(f"{k} = {v:.4f}")