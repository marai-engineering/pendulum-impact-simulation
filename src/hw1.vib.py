import numpy as np
import matplotlib.pyplot as plt


def inertia_bar_plus_tip(m_bar: float, L: float, m_tip: float) -> float:
    """Moment of inertia about the pin for a slender rod about one end + tip mass at the end."""
    return (1/3) * m_bar * L**2 + m_tip * L**2


def com_bar_plus_tip(m_bar: float, L: float, m_tip: float) -> float:
    """Center of mass distance from the pin for rod (COM at L/2) + tip mass at L."""
    m_tot = m_bar + m_tip
    return (m_bar * (L/2) + m_tip * L) / m_tot


def impact_1d_restitution(u1: float, u2: float, mred1: float, mred2: float, e: float) -> tuple[float, float]:
    """
    1D impact along normal direction with restitution e using reduced masses.
    Returns (v1_after, v2_after).
    """
    v1 = (mred1*u1 + mred2*u2 - mred2*e*(u1 - u2)) / (mred1 + mred2)
    v2 = (mred1*u1 + mred2*u2 + mred1*e*(u1 - u2)) / (mred1 + mred2)
    return v1, v2


def simulate_bar2(phi_dot0: float, wn: float, t_end: float = 1.0, n: int = 2001):
    """Analytical free vibration response: phi(t) = (phi_dot0/wn) sin(wn t), with phi(0)=0."""
    t = np.linspace(0.0, t_end, n)
    phi = (phi_dot0 / wn) * np.sin(wn * t)
    phi_dot = phi_dot0 * np.cos(wn * t)
    return t, phi, phi_dot


def main():
    # --- constants / inputs ---
    g = 9.81

    m1 = 1.0
    L1 = 0.7
    m2 = 0.7
    L2 = 0.35
    m0 = 0.2

    k_t = 300.0
    e = 0.6

    phi0_deg = 40.0
    phiB_deg = 40.0
    w2_before = 0.0

    # --- angles ---
    phi0 = np.deg2rad(phi0_deg)
    phiB = np.deg2rad(phiB_deg)

    # --- inertias ---
    theta_O = inertia_bar_plus_tip(m1, L1, m0)
    theta_Q = inertia_bar_plus_tip(m2, L2, m0)

    # --- energies for bar 1 (your convention) ---
    C1 = g * (m1*(L1/2) + m0*L1)
    U0 = -C1 * np.sin(phi0)
    UB = -C1 * np.cos(phiB)
    T_B = U0 - UB

    if T_B <= 0:
        raise ValueError(f"T_B must be positive. Got T_B={T_B:.6f}. Check U0/UB convention.")

    w1_before = np.sqrt((2 * T_B) / theta_O)

    # --- impact geometry: use consistent projection (cos) ---
    r1 = L1 * np.cos(phiB)
    r2 = L2

    u1 = w1_before * r1
    u2 = w2_before * r2

    # --- reduced masses ---
    mred1 = theta_O / (r1**2)
    mred2 = theta_Q / (r2**2)

    v1_after, v2_after = impact_1d_restitution(u1, u2, mred1, mred2, e)

    w1_after = v1_after / r1
    w2_after = v2_after / r2

    # --- swing-back angle for bar 1 using energy after impact ---
    T_after = 0.5 * theta_O * w1_after**2
    U2 = -C1 * np.cos(phiB)
    U3 = T_after + U2

    cos_phiA = -U3 / C1
    cos_phiA = np.clip(cos_phiA, -1.0, 1.0)
    phiA = np.arccos(cos_phiA)
    phiA_deg = np.rad2deg(phiA)

    # --- natural frequency of bar 2 (small vibrations) ---
    m_tot2 = m2 + m0
    s2 = com_bar_plus_tip(m2, L2, m0)              # <-- FIXED
    k_grav = m_tot2 * g * s2
    k_eff = k_t - k_grav

    if k_eff <= 0:
        raise ValueError(f"k_eff must be positive for oscillations. Got k_eff={k_eff:.6f}.")

    wn = np.sqrt(k_eff / theta_Q)

    # --- simulate bar 2 response after impact ---
    t, phi2, phi2_dot = simulate_bar2(phi_dot0=w2_after, wn=wn, t_end=0.5, n=2001)

    # --- prints (clean portfolio style) ---
    print("\n=== RESULTS ===")
    print(f"theta_O = {theta_O:.6f} kg*m^2")
    print(f"theta_Q = {theta_Q:.6f} kg*m^2")
    print(f"U0 = {U0:.6f} J, UB = {UB:.6f} J, T_B = {T_B:.6f} J")
    print(f"omega1_before = {w1_before:.6f} rad/s")
    print(f"omega1_after  = {w1_after:.6f} rad/s")
    print(f"omega2_after  = {w2_after:.6f} rad/s")
    print(f"phiA = {phiA_deg:.6f} deg")
    print(f"s2 = {s2:.6f} m, k_grav = {k_grav:.6f} Nm, k_eff = {k_eff:.6f} Nm")
    print(f"wn = {wn:.6f} rad/s")

    # --- plots (save for GitHub) ---
    plt.figure()
    plt.plot(t, np.rad2deg(phi2))
    plt.xlabel("t [s]")
    plt.ylabel("phi2 [deg]")
    plt.title("Bar 2 free vibration after impact (small-angle)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("bar2_response_phi_deg.png", dpi=200)

    plt.figure()
    plt.plot(t, phi2_dot)
    plt.xlabel("t [s]")
    plt.ylabel("phi2_dot [rad/s]")
    plt.title("Bar 2 angular velocity after impact")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("bar2_response_phi_dot.png", dpi=200)

    plt.show()


if __name__ == "__main__":
    main()