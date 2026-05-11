import Sofa.Core
import numpy as np

# Scale Geomagic workspace (mm) to CTR tip workspace (mm).
# Tune this if the device motion feels too large/small.
SCALE = 1

PRINT_EVERY = 100  # frames between console prints

# Rotation matrix: maps Geomagic delta -> CTR delta.
# Each column answers: "when I move the device along axis X/Y/Z,
# which CTR axis should move and in which direction?"
#
# Identity (no remapping) — tune after running the calibration procedure:
#   1. Hold device still, note the initial print.
#   2. Push device along one axis, watch which CTR coordinate changes.
#   3. Fill in columns accordingly (swap axes or flip signs as needed).
#
# Example — device Z maps to CTR Z, device X maps to CTR -X:
#   FRAME_ROT = np.array([[-1, 0, 0],
#                          [ 0, 1, 0],
#                          [ 0, 0, 1]], dtype=float)
# Calibrated from observations:
#   device forward (-X) -> CTR +Z  (push = deeper along z-axis)
#   device left   (+Z)  -> CTR +X  (bend left in x-z plane)
#   device up     (+Y)  -> CTR +Y  (unchanged)
FRAME_ROT = np.array([[ 0,  1,  0],   # CTR X = device Z
                       [ 0,  0, 1],   # CTR Y = device Y
                       [1,  0,  0]])  # CTR Z = -device X

# Initial IK target: CTR initially extended along Z (x-z plane, tangent = z-axis)
INITIAL_TARGET = [0.0, 0.0, 120.0]


class GeomagicIKBridge(Sofa.Core.Controller):
    """
    Relative mapping: on the first frame, captures the Geomagic initial position
    and the actual CTR tip position from FK. Subsequent frames move the IK target
    by the delta from that initial device position, scaled by SCALE.
    This ensures the device neutral position always maps to the CTR's initial tip.
    """

    def __init__(self, geomagic_driver, ik_controller, fk_engine, target_mo, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.listening = True
        self.device = geomagic_driver
        self.ik = ik_controller
        self.fk = fk_engine
        self.target_mo = target_mo
        self._initialized = False
        self._initial_tip = None
        self._initial_device_pos = None
        self._frame = 0

    def _device_pos(self):
        pose = self.device.positionDevice.value
        return [pose[0], pose[1], pose[2]]

    def onAnimateEndEvent(self, _):
        if self._initialized:
            return
        self._initial_tip = list(INITIAL_TARGET)
        self._initial_device_pos = self._device_pos()
        self._initialized = True
        print(f"[GeomagicIK] Initial CTR tip     : {self._initial_tip}")
        print(f"[GeomagicIK] Initial device pos  : {self._initial_device_pos}")

    def onAnimateBeginEvent(self, _):
        if not self._initialized:
            return
        dp = self._device_pos()
        delta = np.array(dp) - np.array(self._initial_device_pos)
        mapped = FRAME_ROT @ delta * SCALE
        target = (np.array(self._initial_tip) + mapped).tolist()
        self.ik.d_targetPosition.value = target
        self.target_mo.position.value = [target]

        self._frame += 1
        if self._frame % PRINT_EVERY == 0:
            tip = self.fk.d_endEffectorPose.value[:3]
            err = ((target[0]-tip[0])**2 + (target[1]-tip[1])**2 + (target[2]-tip[2])**2) ** 0.5
            print(f"[CTR]  tip: ({tip[0]:7.2f}, {tip[1]:7.2f}, {tip[2]:7.2f}) mm"
                  f"  |  target: ({target[0]:7.2f}, {target[1]:7.2f}, {target[2]:7.2f}) mm"
                  f"  |  err: {err:.2f} mm")


def createScene(rootNode):
    rootNode.dt = 0.01
    rootNode.gravity = [0, 0, 0]

    rootNode.addObject('RequiredPlugin', pluginName=[
        'TorsionRigidModel',
        'Geomagic',
        'Sofa.Component.Visual',
        'Sofa.Component.AnimationLoop',
        'Sofa.Component.StateContainer',
        'Sofa.Component.Collision.Geometry',
    ])

    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showCollisionModels')

    rootNode.addObject('InteractiveCamera', name='camera',
                       position=[0, 0, 300],
                       lookAt=[0, 0, 50],
                       distance=300,
                       fieldOfView=45,
                       zNear=0.1,
                       zFar=2000)

    geomagic = rootNode.addObject('GeomagicDriver',
                                  name='GeomagicDevice',
                                  deviceName='Default Device',
                                  scale=1,
                                  positionBase=[0, 0, 0],
                                  orientationBase=[0, 0.707, 0, -0.707],
                                  drawDeviceFrame=1,
                                  manualStart=False)

    # IK starts with a reasonable initial config: CTR extended along Z
    ik = rootNode.addObject('TRMInverseKinematicsController',
                             name='ik',
                             d_targetPosition=INITIAL_TARGET)

    fk = rootNode.addObject('TRMForwardKinematicsEngine',
                             name='fk',
                             d_jointConfig=ik.d_jointConfig.getLinkPath())

    rootNode.addObject('TRMVisualModel',
                       name='visual',
                       d_jointConfig=ik.d_jointConfig.getLinkPath(),
                       d_nSamples=30)

    target_node = rootNode.addChild('TargetSphere')
    target_mo = target_node.addObject('MechanicalObject', template='Vec3d',
                                      name='targetMO',
                                      position=[INITIAL_TARGET])
    target_node.addObject('SphereCollisionModel', radius=3.0, color='1 0 0 1', group=-1)

    rootNode.addObject(GeomagicIKBridge(geomagic_driver=geomagic,
                                        ik_controller=ik,
                                        fk_engine=fk,
                                        target_mo=target_mo,
                                        name='bridge'))
