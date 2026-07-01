"""
Build a (configuration space, workspace) dataset for the 3-tube CTR by driving
the TRMForwardKinematicsEngine DataEngine over a Sobol low-discrepancy sample
of the 6-DOF joint space.

There is no direct Python binding of CTR::ForwardKinematics::FK() — the only
Python-reachable entry point is the SOFA component TRMForwardKinematicsEngine
(see scenes/CTR_visual_fk_test.py). This script builds a single headless node
hosting that engine, then for each sample q writes fk.d_jointConfig and reads
back fk.d_endEffectorPose and fk.d_manipulability (the DataEngine lazily
recomputes on read since d_jointConfig is a registered input).

d_manipulability is the Yoshikawa measure sqrt(det(J*J^T)), computed in C++
from CTR::InverseKinematics::Jacobian(q) (TRMForwardKinematicsEngine.cpp) —
values near 0 indicate a singular / low-dexterity configuration.

Joint convention (RobotParameters.h): q = [theta1, s1, theta2, s2, theta3, s3]
(rad, mm). Per-tube arc-length limits differ (S3_MIN=25 vs 10 for tubes 1/2).

This script builds its own headless Sofa.Core.Node and calls Sofa.Simulation.init
directly — it is a plain Python script, not a SOFA scene, so it must be run with
`python3` (not `runSofa`, which would try to load it as a scene and treat any
trailing arguments as additional scene files).

Usage:
    SOFA_ROOT=/path/to/sofa/build python3 build_dataset.py --n-samples 100000 --output dataset.npz
"""
import argparse
import math

import numpy as np
from scipy.stats import qmc

import Sofa
import Sofa.Simulation

# Joint limits from RobotParameters.h
S1_MIN, S1_MAX = 10.0, 100.0
S2_MIN, S2_MAX = 10.0, 100.0
S3_MIN, S3_MAX = 25.0, 100.0
THETA_MIN, THETA_MAX = -math.pi, math.pi

# q = [theta1, s1, theta2, s2, theta3, s3]
Q_LOW = np.array([THETA_MIN, S1_MIN, THETA_MIN, S2_MIN, THETA_MIN, S3_MIN])
Q_HIGH = np.array([THETA_MAX, S1_MAX, THETA_MAX, S2_MAX, THETA_MAX, S3_MAX])


def build_fk_node():
    root = Sofa.Core.Node("root")
    root.addObject('RequiredPlugin', name='TorsionRigidModel')
    fk = root.addObject('TRMForwardKinematicsEngine', name='fk',
                         d_jointConfig=[0, 50, 0, 50, 0, 50])
    Sofa.Simulation.init(root)
    return root, fk


def sobol_joint_samples(n_samples, seed=None):
    """Sobol sequence over the 6-D joint box, scaled to per-DOF limits."""
    sampler = qmc.Sobol(d=6, scramble=True, seed=seed)
    m = math.ceil(math.log2(max(n_samples, 1)))
    u = sampler.random_base2(m=m)[:n_samples]  # in [0,1]^6
    return qmc.scale(u, Q_LOW, Q_HIGH)


def build_dataset(n_samples, seed=None):
    _root, fk = build_fk_node()

    q = sobol_joint_samples(n_samples, seed=seed)
    x = np.empty((n_samples, 7), dtype=np.float64)      # [x, y, z, qx, qy, qz, qw]
    manipulability = np.empty(n_samples, dtype=np.float64)  # sqrt(det(J*J^T))

    for i in range(n_samples):
        fk.d_jointConfig.value = q[i]
        x[i] = fk.d_endEffectorPose.value
        manipulability[i] = fk.d_manipulability.value

    return q, x, manipulability


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output", type=str, default="ctr_dataset.npz")
    args = parser.parse_args()

    q, x, manipulability = build_dataset(args.n_samples, seed=args.seed)

    np.savez(args.output, q=q, x=x[:, :3], pose=x, manipulability=manipulability)
    print(f"Saved {q.shape[0]} samples to {args.output} "
          f"(q: {q.shape}, position: {x[:, :3].shape}, pose: {x.shape}, "
          f"manipulability: {manipulability.shape})")


if __name__ == "__main__":
    main()
