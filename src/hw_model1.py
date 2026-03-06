

import numpy as np

g = 9.81   #m/s**2
m1 = 1.0   #kg   
L1 = 0.7   #m  
m2 = 0.7    #kg   
L2 = 0.35  #m
k_t = 300 #Nm


m0 = 0.2       
phi0_deg = 40  
phiB_deg = 40  
e = 0.6
w2_before = 0.0 



phiB = np.deg2rad(phiB_deg)
phi0 = np.deg2rad(phi0_deg)

theta_O = (1/3) * m1 * L1**2 + m0 * L1**2   
theta_Q = (1/3) * m2 * L2**2 + m0 * L2**2  

C1 = g * (m1*(L1/2) + m0*L1) 
U0 = -C1 * np.sin(phi0)
UB = -C1 * np.cos(phiB)

print("U0 =", U0)
print("UB =", UB)

T_B = U0 - UB

w1_before = np.sqrt((2*T_B) / theta_O)
u1 = w1_before * L1 * np.cos(phiB)
u2 = w2_before * L2

r1 = L1 * np.cos(phiB)
r2 = L2

mred1 = theta_O / (r1**2)
mred2 = theta_Q / (r2**2)
v1_after = (mred1*u1 + mred2*u2 - mred2*e*(u1 - u2)) / (mred1 + mred2)
v2_after = (mred1*u1 + mred2*u2 + mred1*e*(u1 - u2)) / (mred1 + mred2)

w1_after = v1_after / r1
w2_after = v2_after / r2



print("\n--- Inertia (theta) ---")
print(f"theta_O = {theta_O:.4f} kg*m^2")
print(f"theta_Q = {theta_Q:.4f} kg*m^2")

print("\n--- Energy ---")
print(f"U0  = {U0:.4f} J")
print(f"UB  = {UB:.4f} J")
print(f"T_B = {T_B:.4f} J")
print(f"omega1_before = {w1_before:.4f} rad/s")

print("\n--- Reduced masses ---")
print(f"r1 = {r1:.4f} m,  mred1 = {mred1:.4f} kg")
print(f"r2 = {r2:.4f} m,  mred2 = {mred2:.4f} kg")

print("\n--- Impact velocities along n ---")
print(f"u1 = {u1:.4f} m/s, u2 = {u2:.4f} m/s")
print(f"v1_after = {v1_after:.4f} m/s, v2_after = {v2_after:.4f} m/s")

print("\n--- Angular velocities after impact ---")
print(f"omega1_after = {w1_after:.4f} rad/s")
print(f"omega2_after = {w2_after:.4f} rad/s")

T_after=0.5*theta_O*w1_after**2

U2 = -C1 * np.cos(phiB)
U3 = T_after + U2    

cos_phiA = -U3 / C1
cos_phiA = np.clip(cos_phiA, -1.0, 1.0)

phiA = np.arccos(cos_phiA)
phiA_deg = np.rad2deg(phiA)

print(f"phiA = {phiA_deg :.4f} deg")

x=phiA_deg/42.503
print(f"x=, {x:.4f}")

m_tot_2 = m2+m0

COM_2 = (m2*(L2/2)+m0*L2)/m_tot_2
k_eff=k_t-m_tot_2*g*COM_2

print(f"keff=, {k_eff:.4f}")

wn= np.sqrt(k_eff / theta_Q)
print(f"wn= {wn:.4f}")

