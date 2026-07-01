"""
Distinguish two causes of "multiple q map to nearly the same x" in the CTR
dataset:

  1. Genuine multiplicity: two well-conditioned (non-singular) configurations
     that legitimately reach (almost) the same tip position via different
     joint configs -- a real alternate IK solution.
  2. Near-singularity blow-up: the local Jacobian is near-singular
     (low manipulability), so a small change in x can correspond to a large
     change in q purely due to numerical/kinematic sensitivity -- not a
     distinct solution branch.

Method: for each sample, find its nearest neighbour in workspace position x
(via a KDTree) and compute how far apart the two samples are in joint space.
A large joint-space distance despite a tiny workspace distance is anomalous;
whether it's (1) or (2) is decided by checking manipulability at the pair.

theta1 is excluded from the joint-space distance: U1F1 = (0,0,0) in
RobotParameters.h means theta1 has zero effect on the tip pose, so any two
samples differing only in theta1 would trivially look like "multiplicity"
and swamp the real signal.

Usage:
    python3 singularity_check.py
"""
import numpy as np
from scipy.spatial import cKDTree

DATASET_PATH = "src/applications/plugins/TorsionRigidModel/scenes/dataset.npz"

# top fraction (by q-jump / x-jump ratio) flagged as anomalous nearest-neighbour pairs
RATIO_PERCENTILE = 99.0
# manipulability percentile below which a sample is considered "near-singular"
SINGULAR_MANIP_PERCENTILE = 10.0


def angular_diff(a, b):
    """Shortest distance between two angles on the circle, in [0, pi]."""
    d = np.abs(a - b) % (2 * np.pi)
    return np.minimum(d, 2 * np.pi - d)


def q_distance(q_a, q_b):
    """Distance in the reduced 5-D joint space [s1, theta2, s2, theta3, s3],
    excluding theta1 (dead DOF, see module docstring)."""
    d_s1 = q_a[:, 1] - q_b[:, 1]
    d_s2 = q_a[:, 3] - q_b[:, 3]
    d_s3 = q_a[:, 5] - q_b[:, 5]
    d_t2 = angular_diff(q_a[:, 2], q_b[:, 2])
    d_t3 = angular_diff(q_a[:, 4], q_b[:, 4])
    return np.sqrt(d_s1**2 + d_s2**2 + d_s3**2 + d_t2**2 + d_t3**2)


def main():
    dataset = np.load(DATASET_PATH)
    q = dataset["q"]
    x = dataset["x"]
    manipulability = dataset["manipulability"]
    n = len(x)

    tree = cKDTree(x)
    dist_x, nn_idx = tree.query(x, k=2)  # k=1 is self (dist 0); k=2 is nearest neighbour
    dist_x = dist_x[:, 1]
    nn_idx = nn_idx[:, 1]

    dist_q = q_distance(q, q[nn_idx])
    ratio = dist_q / np.maximum(dist_x, 1e-9)

    ratio_thresh = np.percentile(ratio, RATIO_PERCENTILE)
    manip_thresh = np.percentile(manipulability, SINGULAR_MANIP_PERCENTILE)

    anomalous = ratio >= ratio_thresh
    pair_min_manip = np.minimum(manipulability, manipulability[nn_idx])

    near_singular = anomalous & (pair_min_manip < manip_thresh)
    genuine_multiplicity = anomalous & ~near_singular

    print(f"N samples: {n}")
    print(f"dist_x (nearest neighbour)          min/median/max: "
          f"{dist_x.min():.4f} / {np.median(dist_x):.4f} / {dist_x.max():.4f}")
    print(f"dist_q (reduced, excl. theta1)       min/median/max: "
          f"{dist_q.min():.4f} / {np.median(dist_q):.4f} / {dist_q.max():.4f}")
    print(f"ratio dist_q/dist_x  {RATIO_PERCENTILE:.1f}th percentile threshold: {ratio_thresh:.4f}")
    print(f"manipulability       {SINGULAR_MANIP_PERCENTILE:.1f}th percentile threshold: {manip_thresh:.4f}")
    print()
    print(f"Anomalous nearest-neighbour pairs (top {100 - RATIO_PERCENTILE:.1f}%): {anomalous.sum()}")
    print(f"  -> explained by near-singularity (low manipulability): {near_singular.sum()}")
    print(f"  -> genuine multiplicity (manipulability NOT low):      {genuine_multiplicity.sum()}")

    if genuine_multiplicity.sum() > 0:
        idx = np.where(genuine_multiplicity)[0]
        order = np.argsort(-ratio[idx])[:10]
        print("\nTop genuine-multiplicity candidates (i, neighbour j, dist_x, dist_q, manip_i, manip_j):")
        for k in order:
            i = idx[k]
            j = nn_idx[i]
            print(f"  i={i:6d} j={j:6d}  dist_x={dist_x[i]:.4f}  dist_q={dist_q[i]:.4f}  "
                  f"manip_i={manipulability[i]:.3f}  manip_j={manipulability[j]:.3f}")

    if near_singular.sum() > 0:
        idx = np.where(near_singular)[0]
        order = np.argsort(-ratio[idx])[:10]
        print("\nTop near-singularity candidates (i, neighbour j, dist_x, dist_q, manip_i, manip_j):")
        for k in order:
            i = idx[k]
            j = nn_idx[i]
            print(f"  i={i:6d} j={j:6d}  dist_x={dist_x[i]:.4f}  dist_q={dist_q[i]:.4f}  "
                  f"manip_i={manipulability[i]:.3f}  manip_j={manipulability[j]:.3f}")


if __name__ == "__main__":
    main()
