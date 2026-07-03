"""
Overlay two joint-space solutions (q1, q2) that map to (nearly) the same
end-effector position, found by singularity_check.py, as two CTR backbones
in the same SOFA scene.

[verify_multiplicity] target idx=70132  manipulability=1.7088e+03  (>= p50 of dataset)
  x_target = [  8.81241693,  26.73094764, 167.88485817]
  q_target = [-1.27955292, 32.35630317,  1.27590315, 94.7628035, -0.3715183,  46.95589992]
    q = [ 4.56071326, 26.85460515,  0.95307159, 86.14985872, -0.51352684, 60.80439022]

Run with runSofa (this is a GUI scene, not a headless script):
    runSofa -l SofaPython3 src/applications/plugins/TorsionRigidModel/scenes/visualise_multiplicity.py
"""
import os

import numpy as np

DATASET_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset.npz")

EXAMPLES = {
    "multiplicity": (22248, 99591),
    "near_singularity": (4736424, 7665160),
    "confirmed_full_pose": (
        [ 1.14387978, 52.14283627, -2.60092807, 71.24652113,  2.01370363, 64.13945619], 
        [ 0.84345381, 41.80899274, -3.12119069, 65.22329285, -4.37676958, 80.30195603],)
}
CASE = "confirmed_full_pose"  # or "multiplicity" or "near_singularity"
ENTRY = EXAMPLES[CASE]


def fmt(arr):
    """Comma-separated values, e.g. for pasting into EXAMPLES above."""
    return "[" + ", ".join(f"{v:.6f}" for v in np.asarray(arr)) + "]"


def createScene(rootNode):
    rootNode.dt = 0.01
    rootNode.gravity = [0, 0, 0]

    rootNode.addObject('RequiredPlugin', name='TorsionRigidModel')
    rootNode.addObject('RequiredPlugin', name='Sofa.Component.Visual')
    rootNode.addObject('RequiredPlugin', name='Sofa.Component.AnimationLoop')

    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels')

    rootNode.addObject('InteractiveCamera', name='camera',
                       position=[0, 0, 300],
                       lookAt=[0, 0, 50],
                       distance=300,
                       fieldOfView=45,
                       zNear=0.1,
                       zFar=2000)

    if isinstance(ENTRY[0], (int, np.integer)):
        # Dataset-index case ("multiplicity" / "near_singularity"): look up
        # both joint configs, plus their pose/manipulability, from dataset.npz.
        i1, i2 = ENTRY
        dataset = np.load(DATASET_PATH)
        q1 = dataset["q"][i1]
        q2 = dataset["q"][i2]
        x1 = dataset["x"][i1]
        x2 = dataset["x"][i2]
        pose1 = dataset["pose"][i1]
        pose2 = dataset["pose"][i2]
        manip1 = dataset["manipulability"][i1]
        manip2 = dataset["manipulability"][i2]

        print(f"[visualise_multiplicity] case = {CASE}")
        print(f"[visualise_multiplicity] q1 (i={i1}) = {fmt(q1)}  manipulability={manip1:.3f}")
        print(f"[visualise_multiplicity] x1 (i={i1}) = {fmt(x1)}")
        print(f"[visualise_multiplicity] pose1 [x,y,z,qx,qy,qz,qw] = {fmt(pose1)}")
        print(f"[visualise_multiplicity] q2 (j={i2}) = {fmt(q2)}  manipulability={manip2:.3f}")
        print(f"[visualise_multiplicity] x2 (j={i2}) = {fmt(x2)}")
        print(f"[visualise_multiplicity] pose2 [x,y,z,qx,qy,qz,qw] = {fmt(pose2)}")
    else:
        # Hardcoded-q case ("confirmed_full_pose"): q1/q2 are literal joint
        # configs (e.g. a verify_multiplicity.py solve result), not dataset
        # samples -- no pose/manipulability lookup available for these.
        q1 = np.array(ENTRY[0])
        q2 = np.array(ENTRY[1])

        print(f"[visualise_multiplicity] case = {CASE}  (hardcoded q vectors)")
        print(f"[visualise_multiplicity] q1 = {fmt(q1)}")
        print(f"[visualise_multiplicity] q2 = {fmt(q2)}")

    # CTR 1: default blue/green/red per-section colouring, opaque
    rootNode.addObject('TRMVisualModel',
                       name='ctr1_solid',
                       d_jointConfig=list(q1),
                       d_nSamples=30)

    # CTR 2: flat orange, semi-transparent ghost overlay of the alternate solution
    rootNode.addObject('TRMVisualModel',
                       name='ctr2_ghost',
                       d_jointConfig=list(q2),
                       d_nSamples=30,
                       d_useFlatColor=True,
                       d_color=[1.0, 0.5, 0.0, 0.4])


# runSofa -l SofaPython3 src/applications/plugins/TorsionRigidModel/scenes/visualise_multiplicity.py