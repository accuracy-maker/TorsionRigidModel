"""
Verify full-pose IK multiplicity by random seeding.

Method:
  1. Pick ONE target sample from a "high feasible region" of dataset.npz --
     the subset of samples whose manipulability is at or above
     MANIP_PERCENTILE, i.e. away from kinematic singularities where a
     multiplicity claim would be confounded with ill-conditioning -- and
     choose one of those uniformly at random. This is an O(1) filter +
     random pick, not a search over the whole dataset.
  2. Draw N_SEEDS random joint configurations uniformly over the full joint
     box (same bounds build_dataset.py's Sobol sampler uses), independent of
     the dataset.
  3. From each random seed, run a full 6-DOF pose DLS solve (position AND
     orientation, driven by a numerical 6x6 Jacobian -- reimplemented here
     since neither the analytic Jacobian nor a scriptable IK entry point
     with a settable seed is exposed to Python; see CTR_visual_ik_test.py)
     toward the target's pose.
  4. A solve that converges to (near) machine precision on both position and
     orientation, but lands far from q_target in joint space, is a distinct
     branch of the FK map reaching the exact same SE(3) pose -- a second
     seed converging to the same target from a different joint config is
     what "distinct branches" means here, not a numerical artifact of a
     shared local IK solution.

Run headless (like build_dataset.py, NOT runSofa -- there is no GUI here):
    SOFA_ROOT=/path/to/sofa/build python3 verify_multiplicity.py
"""
import math
import os

import numpy as np
from scipy.spatial.transform import Rotation

import Sofa
import Sofa.Simulation

DATASET_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset.npz")

# Joint limits from RobotParameters.h -- same box build_dataset.py samples from.
S1_MIN, S1_MAX = 10.0, 100.0
S2_MIN, S2_MAX = 10.0, 100.0
S3_MIN, S3_MAX = 25.0, 100.0
THETA_MIN, THETA_MAX = -math.pi, math.pi
Q_LOW = np.array([THETA_MIN, S1_MIN, THETA_MIN, S2_MIN, THETA_MIN, S3_MIN])
Q_HIGH = np.array([THETA_MAX, S1_MAX, THETA_MAX, S2_MAX, THETA_MAX, S3_MAX])

# "High feasible region": only target samples whose manipulability is at or
# above this percentile of the dataset are eligible to be picked as the
# fixed target.
MANIP_PERCENTILE = 50

N_SEEDS = 100
RNG_SEED = None  # set to an int for a reproducible run

MAX_ITERS = 1000
POS_TOL_MM = 1e-6
ROT_TOL_RAD = 1e-6
FD_EPS = 1e-6

# A converged solve landing within this joint-space distance of q_target is
# treated as "collapsed back onto the target", not a distinct branch.
DIST_Q_COLLAPSE_THRESH = 1.0


def build_fk_node():
    root = Sofa.Core.Node("root")
    root.addObject('RequiredPlugin', name='TorsionRigidModel')
    fk = root.addObject('TRMForwardKinematicsEngine', name='fk',
                         d_jointConfig=[0, 50, 0, 50, 0, 50])
    Sofa.Simulation.init(root)
    return root, fk


def fk_pose(fk, q):
    """[x, y, z, qx, qy, qz, qw] at joint config q."""
    fk.d_jointConfig.value = q
    return np.array(fk.d_endEffectorPose.value, dtype=np.float64)


def clip_lengths(q):
    q = q.copy()
    q[1] = np.clip(q[1], S1_MIN, S1_MAX)
    q[3] = np.clip(q[3], S2_MIN, S2_MAX)
    q[5] = np.clip(q[5], S3_MIN, S3_MAX)
    return q


def q_distance(q_a, q_b):
    """Wrapped angular distance for the theta_i, plain difference for the s_i."""
    theta_idx = [0, 2, 4]
    length_idx = [1, 3, 5]
    d_theta_raw = np.abs(q_a[theta_idx] - q_b[theta_idx]) % (2 * np.pi)
    d_theta = np.minimum(d_theta_raw, 2 * np.pi - d_theta_raw)
    d_s = q_a[length_idx] - q_b[length_idx]
    return float(np.sqrt(np.sum(d_theta**2) + np.sum(d_s**2)))


def dedupe_qs(qs, thresh):
    """Collapse a list of q's into ones that are pairwise more than thresh
    apart in joint space, so repeat convergence onto the same branch is
    only counted once."""
    unique = []
    for q in qs:
        if all(q_distance(q, u) > thresh for u in unique):
            unique.append(q)
    return unique


def pick_target(dataset, rng):
    """Pick one sample index from the high-feasibility region of the
    dataset -- an O(1) filter + random choice, not a search."""
    manipulability = dataset["manipulability"]
    manip_thresh = np.percentile(manipulability, MANIP_PERCENTILE)
    eligible = np.where(manipulability >= manip_thresh)[0]
    return int(rng.choice(eligible))


def random_seed_q(rng):
    """A joint configuration drawn uniformly over the full joint box,
    independent of the dataset."""
    return rng.uniform(Q_LOW, Q_HIGH)


# --------------------------------------------------------------------------
# Full 6-DOF pose DLS (position + orientation)
# --------------------------------------------------------------------------

def orientation_error(quat_current, quat_target):
    """Rotation vector (rad) that rotates quat_current onto quat_target,
    expressed in quat_current's frame -- the log map of the relative rotation."""
    rel = Rotation.from_quat(quat_current).inv() * Rotation.from_quat(quat_target)
    return rel.as_rotvec()


def pose_jacobian(fk, q):
    """6x6 numerical Jacobian: rows [vx,vy,vz, wx,wy,wz] w.r.t. q (central
    difference). The angular rows are the finite-difference log map of the
    relative rotation between the +/-eps perturbations -- consistent with
    orientation_error's convention since both are evaluated in the same
    local neighbourhood of the current q."""
    J = np.empty((6, 6))
    for k in range(6):
        dq = np.zeros(6)
        dq[k] = FD_EPS
        pose_plus = fk_pose(fk, q + dq)
        pose_minus = fk_pose(fk, q - dq)
        J[:3, k] = (pose_plus[:3] - pose_minus[:3]) / (2 * FD_EPS)
        rel = Rotation.from_quat(pose_minus[3:]).inv() * Rotation.from_quat(pose_plus[3:])
        J[3:, k] = rel.as_rotvec() / (2 * FD_EPS)
    return J


def ik_step_pose(fk, q, target_position, target_quat):
    pose = fk_pose(fk, q)
    e = np.concatenate([target_position - pose[:3], orientation_error(pose[3:], target_quat)])
    J = pose_jacobian(fk, q)
    # Standard Levenberg-Marquardt damping -- theta_i must move freely here
    # since it directly drives the orientation target.
    H = J @ J.T + 0.5 * float(e @ e) * np.eye(6)
    delta_q = J.T @ np.linalg.solve(H, e)
    return clip_lengths(q + delta_q), e


def solve_ik_pose(fk, q0, target_position, target_quat):
    q = q0.copy()
    pos_residual = rot_residual = float("inf")
    for it in range(MAX_ITERS):
        q, e = ik_step_pose(fk, q, target_position, target_quat)
        pos_residual = float(np.linalg.norm(e[:3]))
        rot_residual = float(np.linalg.norm(e[3:]))
        if pos_residual < POS_TOL_MM and rot_residual < ROT_TOL_RAD:
            return q, pos_residual, rot_residual, it + 1
    return q, pos_residual, rot_residual, MAX_ITERS


# --------------------------------------------------------------------------
# Batch driver
# --------------------------------------------------------------------------

def test_seed(fk, q_target, q_seed, target_position, target_quat):
    dist_q_before = q_distance(q_seed, q_target)
    q_final, pos_res, rot_res, iters = solve_ik_pose(fk, q_seed, target_position, target_quat)
    dist_q_after = q_distance(q_final, q_target)

    converged = pos_res < POS_TOL_MM and rot_res < ROT_TOL_RAD
    if converged and dist_q_after > DIST_Q_COLLAPSE_THRESH:
        verdict = "CONFIRMED distinct branch"
    elif converged:
        verdict = "collapsed onto target"
    else:
        verdict = "did not converge"

    return {
        "dist_q_before": dist_q_before,
        "q_final": q_final, "pos_res": pos_res, "rot_res": rot_res,
        "iters": iters, "dist_q_after": dist_q_after,
        "verdict": verdict,
    }


def main():
    dataset = np.load(DATASET_PATH)
    q_all = dataset["q"]
    pose_all = dataset["pose"]

    rng = np.random.default_rng(RNG_SEED)

    target_idx = pick_target(dataset, rng)
    q_target = q_all[target_idx]
    target_position = pose_all[target_idx, :3]
    target_quat = pose_all[target_idx, 3:]

    print(f"[verify_multiplicity] target idx={target_idx}  "
          f"manipulability={dataset['manipulability'][target_idx]:.4e}  "
          f"(>= p{MANIP_PERCENTILE} of dataset)")
    print(f"  x_target = {target_position}")
    print(f"  q_target = {q_target}")
    print(f"  testing {N_SEEDS} random joint-space seeds against this single target\n")

    _root, fk = build_fk_node()

    results = []
    for trial in range(N_SEEDS):
        q_seed = random_seed_q(rng)
        r = test_seed(fk, q_target, q_seed, target_position, target_quat)
        results.append(r)

        print(f"[{trial + 1}/{N_SEEDS}] dist_q(seed,target)={r['dist_q_before']:6.2f}  "
              f"pos_res={r['pos_res']:.2e}mm  rot_res={np.degrees(r['rot_res']):.2e}deg  "
              f"iters={r['iters']:4d}  dist_q(final,target)={r['dist_q_after']:.4f}  "
              f"-> {r['verdict']}")

    n_confirmed = sum(1 for r in results if r["verdict"] == "CONFIRMED distinct branch")
    n_collapsed = sum(1 for r in results if r["verdict"] == "collapsed onto target")
    n_noconverge = sum(1 for r in results if r["verdict"] == "did not converge")

    confirmed_qs = [r["q_final"] for r in results if r["verdict"] == "CONFIRMED distinct branch"]
    distinct_qs = dedupe_qs(confirmed_qs, DIST_Q_COLLAPSE_THRESH)

    print("\n" + "=" * 70)
    print(f"Summary over {N_SEEDS} random seeds against target idx={target_idx}:")
    print(f"  confirmed distinct branch : {n_confirmed}")
    print(f"  collapsed onto target     : {n_collapsed}")
    print(f"  did not converge          : {n_noconverge}")
    print(f"\n  x_target = {target_position}")
    print(f"  q_target = {q_target}")
    if distinct_qs:
        most_distinct_q = max(distinct_qs, key=lambda q: q_distance(q, q_target))
        print(f"  most distinct q ({len(distinct_qs)} branch(es) found, "
              f"dist_q={q_distance(most_distinct_q, q_target):.4f}): {most_distinct_q}")
    else:
        print("  no distinct branch found")


if __name__ == "__main__":
    main()
