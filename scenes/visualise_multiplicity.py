"""
Overlay two joint-space solutions (q1, q2) that map to (nearly) the same
end-effector position, found by singularity_check.py, as two CTR backbones
in the same SOFA scene.

CTR 1 (q1): default per-section blue/green/red colouring, opaque.
CTR 2 (q2): flat orange, semi-transparent -- overlaid as a "ghost" solution.

Set CASE below to pick which of singularity_check.py's two failure modes to
look at:
  "multiplicity"      i=81389 / j=56213 -- dist_x=0.0244, dist_q=25.73,
                       manip_i=6.81, manip_j=8.38 (both well above the
                       singularity threshold ~5.2: a genuine alternate
                       solution, not a conditioning artifact).
  "near_singularity"  i=17438 / j=84653 -- dist_x=0.0915, dist_q=58.15,
                       manip_i=2.40, manip_j=4.64 (both below the
                       singularity threshold: the large q-jump is explained
                       by the local Jacobian being near-singular here, not
                       a distinct solution branch).

Run with runSofa (this is a GUI scene, not a headless script):
    runSofa -l SofaPython3 src/applications/plugins/TorsionRigidModel/scenes/visualise_multiplicity.py
"""
import os

import numpy as np

DATASET_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset.npz")

EXAMPLES = {
    "multiplicity": (75187, 68921),
    "near_singularity": (17438, 84653),
}
CASE = "multiplicity"  # or "near_singularity"
I1, I2 = EXAMPLES[CASE]


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

    dataset = np.load(DATASET_PATH)
    q1 = dataset["q"][I1]
    q2 = dataset["q"][I2]
    manip1 = dataset["manipulability"][I1]
    manip2 = dataset["manipulability"][I2]

    print(f"[visualise_multiplicity] case = {CASE}")
    print(f"[visualise_multiplicity] q1 (i={I1}) = {q1}  manipulability={manip1:.3f}")
    print(f"[visualise_multiplicity] q2 (j={I2}) = {q2}  manipulability={manip2:.3f}")

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