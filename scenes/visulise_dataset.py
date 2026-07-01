import numpy as np
import matplotlib.pyplot as plt

dataset = np.load("src/applications/plugins/TorsionRigidModel/scenes/dataset.npz")
q = dataset["q"]                       # (N, 6) joint config [theta1, s1, theta2, s2, theta3, s3]
x = dataset["x"]                       # (N, 3) tip position [x, y, z]
pose = dataset["pose"]                 # (N, 7) tip position + quaternion [x, y, z, qx, qy, qz, qw]
manipulability = dataset["manipulability"]  # (N,) sqrt(det(Jp*Jp^T)), tip-point singularity measure

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
sc = ax.scatter(x[:, 0], x[:, 1], x[:, 2], c=manipulability, cmap="viridis", s=2, alpha=0.5)
fig.colorbar(sc, ax=ax, shrink=0.6, label="manipulability sqrt(det(Jp*Jp^T))")
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_zlabel("z (mm)")
ax.set_title("CTR workspace (tip position, coloured by manipulability)")

joint_labels = ["theta1 (rad)", "s1 (mm)", "theta2 (rad)", "s2 (mm)", "theta3 (rad)", "s3 (mm)"]
fig2, axes = plt.subplots(2, 3, figsize=(12, 6))
for i, ax_i in enumerate(axes.flat):
    ax_i.hist(q[:, i], bins=50)
    ax_i.set_xlabel(joint_labels[i])
    ax_i.set_ylabel("count")
fig2.suptitle("Joint configuration space distributions")
fig2.tight_layout()

fig3, ax3 = plt.subplots()
ax3.hist(manipulability, bins=100)
ax3.set_xlabel("manipulability sqrt(det(Jp*Jp^T))")
ax3.set_ylabel("count")
ax3.set_title("Manipulability distribution")

plt.show()