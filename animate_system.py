import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

from src.impact_solver import solve_case


def smooth_angle(theta0, w0, theta1, w1, T, t):
    s = np.clip(t / T, 0.0, 1.0)

    h00 = 2*s**3 - 3*s**2 + 1
    h10 = s**3 - 2*s**2 + s
    h01 = -2*s**3 + 3*s**2
    h11 = s**3 - s**2

    return h00*theta0 + h10*T*w0 + h01*theta1 + h11*T*w1


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

L1 = case["L1"]
L2 = case["L2"]

phi0 = np.deg2rad(case["phi0_deg"])
phiB = np.deg2rad(case["phiB_deg"])
phiA = np.deg2rad(res["phiA_deg"])

w1_before = res["w1_before"]
w1_after = res["w1_after"]
w2_after = res["w2_after"]
wn = res["wn"]

theta1_init = -(np.pi/2 - phi0)
theta1_impact = phiB
theta1_return = theta1_impact - phiA

T_pre = 2 * abs(theta1_impact - theta1_init) / w1_before
T_post = 2 * abs(theta1_impact - theta1_return) / w1_after
T_extra = 1.2
T_total = T_pre + T_post + T_extra

O = np.array([0.0, 0.0])

B = np.array([
    L1 * np.sin(theta1_impact),
    -L1 * np.cos(theta1_impact)
])

Q = B + np.array([0.0, -L2])


def theta1_t(t):
    if t <= T_pre:
        return smooth_angle(theta1_init, 0.0, theta1_impact, w1_before, T_pre, t)
    elif t <= T_pre + T_post:
        tau = t - T_pre
        return smooth_angle(theta1_impact, -w1_after, theta1_return, 0.0, T_post, tau)
    else:
        return theta1_return


def theta2_t(t):
    if t <= T_pre:
        return 0.0
    tau = t - T_pre
    return (w2_after / wn) * np.sin(wn * tau)


def tip1(theta):
    return np.array([
        L1 * np.sin(theta),
        -L1 * np.cos(theta)
    ])


def tip2(theta):
    return Q + np.array([
        L2 * np.sin(theta),
        L2 * np.cos(theta)
    ])


fig, ax = plt.subplots(figsize=(7, 7))
ax.set_aspect("equal")
ax.grid(True)

ax.set_title("Two-bar impact dynamics animation")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")

ax.set_xlim(-L1 - 0.35, B[0] + 0.35)
ax.set_ylim(min(-L1 - 0.2, Q[1] - 0.2), 0.2)

ax.plot(
    [B[0] - 0.15, B[0] + 0.25],
    [B[1], B[1]],
    "--",
    linewidth=1.2,
    label="impact line n"
)

ax.plot(O[0], O[1], "ko", markersize=6)
ax.text(O[0] + 0.02, O[1] + 0.02, "O")

ax.plot(Q[0], Q[1], "ko", markersize=6)
ax.text(Q[0] + 0.02, Q[1] - 0.04, "Q")

line1, = ax.plot([], [], linewidth=3, label="bar 1")
line2, = ax.plot([], [], linewidth=3, label="bar 2")

mass1, = ax.plot([], [], "o", markersize=12)
mass2, = ax.plot([], [], "o", markersize=12)

# dashed rotational paths
path1, = ax.plot([], [], "--", linewidth=1.5)
path2, = ax.plot([], [], "--", linewidth=1.5)

time_text = ax.text(0.02, 0.96, "", transform=ax.transAxes, va="top")
phase_text = ax.text(0.02, 0.90, "", transform=ax.transAxes, va="top")
angle_text = ax.text(0.02, 0.84, "", transform=ax.transAxes, va="top")

ax.legend(loc="lower left")

fps = 40
n_frames = int(T_total * fps)
t_vals = np.linspace(0.0, T_total, n_frames)

# full rotational paths
ang1 = np.linspace(min(theta1_init, theta1_return) - 0.1, theta1_impact + 0.1, 200)
x_path1 = L1 * np.sin(ang1)
y_path1 = -L1 * np.cos(ang1)

ang2 = np.linspace(-0.35, 0.35, 200)
x_path2 = Q[0] + L2 * np.sin(ang2)
y_path2 = Q[1] + L2 * np.cos(ang2)


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    mass1.set_data([], [])
    mass2.set_data([], [])
    path1.set_data(x_path1, y_path1)
    path2.set_data(x_path2, y_path2)
    time_text.set_text("")
    phase_text.set_text("")
    angle_text.set_text("")
    return line1, line2, mass1, mass2, path1, path2, time_text, phase_text, angle_text


def update(frame):
    t = t_vals[frame]

    th1 = theta1_t(t)
    th2 = theta2_t(t)

    P1 = tip1(th1)
    P2 = tip2(th2)

    line1.set_data([O[0], P1[0]], [O[1], P1[1]])
    line2.set_data([Q[0], P2[0]], [Q[1], P2[1]])

    mass1.set_data([P1[0]], [P1[1]])
    mass2.set_data([P2[0]], [P2[1]])

    time_text.set_text(f"t = {t:.2f} s")

    if t <= T_pre:
        phase_text.set_text("Phase: before impact")
        angle_text.set_text(f"bar 1 angle = {np.rad2deg(th1):.1f} deg")
    elif t <= T_pre + T_post:
        phase_text.set_text("Phase: after impact")
        angle_text.set_text(f"bar 1 angle = {np.rad2deg(th1):.1f} deg")
    else:
        phase_text.set_text("Phase: bar 2 vibration continues")
        angle_text.set_text(f"bar 2 angle = {np.rad2deg(th2):.2f} deg")

    return line1, line2, mass1, mass2, path1, path2, time_text, phase_text, angle_text


anim = FuncAnimation(
    fig,
    update,
    frames=n_frames,
    init_func=init,
    interval=1000 / fps,
    blit=True
)

Path("figures").mkdir(exist_ok=True)
gif_path = "figures/two_bar_impact_animation.gif"
anim.save(gif_path, writer=PillowWriter(fps=fps))

print(f"\nAnimation saved to: {gif_path}")

plt.show()