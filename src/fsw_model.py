import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# === PARAMETERS ===
omega = 100.0                 # Angular speed (rad/s)
Rs = 5.0                     # mm, shoulder outer radius
Rp = 1                     # mm, shoulder inner radius
tau_Cu = 15e6                  # Pa, shear stress for material A
tau_Al = 10e6                  # Pa, shear stress for material B
N_r = 150                     # Radial grid points
N_theta = 150                 # Angular grid points
k = 0.1                   # Parabolic constant (mm^-1) for z = -k*(r-Rp)*(r-Rs)

# Convert mm to meters for area integration
Rs_m = Rs * 1e-3
Rp_m = Rp * 1e-3

# === SHOULDER (inverted parabolic surface) ===
r_sh = np.linspace(Rp, Rs, N_r)          # mm
dr = (Rs - Rp) / (N_r - 1)               # mm
dr_m = dr * 1e-3                         # m

# === Material A: theta in [pi/2, 3pi/2] ===
theta_A = np.linspace(np.pi/2, 3*np.pi/2, N_theta)
dtheta_A = (theta_A[-1] - theta_A[0]) / (N_theta - 1)

R_sh_A, Theta_sh_A = np.meshgrid(r_sh, theta_A)
X_sh_A = R_sh_A * np.cos(Theta_sh_A)
Y_sh_A = R_sh_A * np.sin(Theta_sh_A)
Z_sh_A = -k * (R_sh_A - Rp) * (R_sh_A - Rs)

# Heat flux density (W/m^2)
q_A = omega * tau_Cu* (R_sh_A * 1e-3)

# Incremental heat generation (W)
dQ_A = omega * tau_Cu* (R_sh_A * 1e-3)**2 * dr_m * dtheta_A
Q_total_A = np.sum(dQ_A)

# === Material B: theta in [0, pi/2] and [3pi/2, 2pi] ===
theta_B1 = np.linspace(0, np.pi/2, N_theta//2)
theta_B2 = np.linspace(3*np.pi/2, 2*np.pi, N_theta//2)
dtheta_B = (np.pi/2) / (N_theta//2 - 1)

# Segment 1
R_sh_B1, Theta_sh_B1 = np.meshgrid(r_sh, theta_B1)
X_sh_B1 = R_sh_B1 * np.cos(Theta_sh_B1)
Y_sh_B1 = R_sh_B1 * np.sin(Theta_sh_B1)
Z_sh_B1 = -k * (R_sh_B1 - Rp) * (R_sh_B1 - Rs)
q_B1 = omega * tau_Al * (R_sh_B1 * 1e-3)
dQ_B1 = omega * tau_Al * (R_sh_B1 * 1e-3)**2 * dr_m * dtheta_B
Q_total_B1 = np.sum(dQ_B1)

# Segment 2
R_sh_B2, Theta_sh_B2 = np.meshgrid(r_sh, theta_B2)
X_sh_B2 = R_sh_B2 * np.cos(Theta_sh_B2)
Y_sh_B2 = R_sh_B2 * np.sin(Theta_sh_B2)
Z_sh_B2 = -k * (R_sh_B2 - Rp) * (R_sh_B2 - Rs)
q_B2 = omega * tau_Al * (R_sh_B2 * 1e-3)
dQ_B2 = omega * tau_Al * (R_sh_B2 * 1e-3)**2 * dr_m * dtheta_B
Q_total_B2 = np.sum(dQ_B2)

Q_total_B = Q_total_B1 + Q_total_B2

# === GLOBAL NORMALIZATION FOR COLORMAP ===
q_max = max(q_A.max(), q_B1.max(), q_B2.max())

# === PLOT FOR MATERIAL A ===
fig_A = plt.figure(figsize=(10, 8))
ax_A = fig_A.add_subplot(111, projection='3d')
surf_A = ax_A.plot_surface(X_sh_A, Y_sh_A, Z_sh_A, facecolors=plt.cm.hot(q_A / q_max),
                          rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

mappable_A = plt.cm.ScalarMappable(cmap='hot')
mappable_A.set_array(q_A)
fig_A.colorbar(mappable_A, ax=ax_A, shrink=0.6, label="Heat Flux (W/m²)")

ax_A.set_title(f"Heat Map on Inverted Parabolic Shoulder (Cu)")
ax_A.set_xlabel("X (mm)")
ax_A.set_ylabel("Y (mm)")
ax_A.set_zlabel("Z (mm)")
ax_A.view_init(elev=30, azim=60)

# Equal aspect ratio
max_range = np.array([X_sh_A.max() - X_sh_A.min(),
                      Y_sh_A.max() - Y_sh_A.min(),
                      abs(Z_sh_A.min())]).max() / 2.0
mid_x = (X_sh_A.max() + X_sh_A.min()) * 0.5
mid_y = (Y_sh_A.max() + Y_sh_A.min()) * 0.5
mid_z = Z_sh_A.min() / 2.0
ax_A.set_xlim(mid_x - max_range, mid_x + max_range)
ax_A.set_ylim(mid_y - max_range, mid_y + max_range)
ax_A.set_zlim(Z_sh_A.min() - max_range, max_range)

plt.tight_layout()

# === PLOT FOR MATERIAL B ===
fig_B = plt.figure(figsize=(10, 8))
ax_B = fig_B.add_subplot(111, projection='3d')

surf_B1 = ax_B.plot_surface(X_sh_B1, Y_sh_B1, Z_sh_B1, facecolors=plt.cm.hot(q_B1 / q_max),
                           rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
surf_B2 = ax_B.plot_surface(X_sh_B2, Y_sh_B2, Z_sh_B2, facecolors=plt.cm.hot(q_B2 / q_max),
                           rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

mappable_B = plt.cm.ScalarMappable(cmap='hot')
mappable_B.set_array(np.concatenate([q_B1.ravel(), q_B2.ravel()]))
fig_B.colorbar(mappable_B, ax=ax_B, shrink=0.6, label="Heat Flux (W/m²)")

ax_B.set_title(f"Heat Map on Inverted Parabolic Shoulder (Al)")
ax_B.set_xlabel("X (mm)")
ax_B.set_ylabel("Y (mm)")
ax_B.set_zlabel("Z (mm)")
ax_B.view_init(elev=30, azim=60)

ax_B.set_xlim(mid_x - max_range, mid_x + max_range)
ax_B.set_ylim(mid_y - max_range, mid_y + max_range)
ax_B.set_zlim(Z_sh_B1.min() - max_range, max_range)

plt.tight_layout()
plt.show()
# === PRINT RESULTS ===
print("\n=== HEAT GENERATION RESULTS ===")
print(f"Material A side:")
print(f"  Total heat        = {Q_total_A/1e3:.3f} kW")
# Need to calculate min/max/avg flux for A and B sides separately
# ... (calculation code would go here)

print(f"Material B side:")
print(f"  Total heat        = {Q_total_B/1e3:.3f} kW")
# Need to calculate min/max/avg flux for B side
# ... (calculation code would go here)

print(f"Combined total heat = {(Q_total_A + Q_total_B)/1e3:.3f} kW")

# === PRINT RESULTS ===
print(f"Total heat generated on Material A side = {Q_total_A/1e3:.2f} kW")
print(f"Total heat generated on Material B side = {Q_total_B/1e3:.2f} kW")
print(f"Total heat generated on both sides     = {(Q_total_A+Q_total_B)/1e3:.2f} kW")
