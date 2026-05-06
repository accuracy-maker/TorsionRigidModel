import Sofa.Core
from Sofa.constants import Key

DP = 5.0          # mm per keypress
PRINT_EVERY = 50  # frames between console prints


class IKTargetController(Sofa.Core.Controller):
    """
    Keyboard control for IK target position.

    Key layout:
      1 / 2 / 3    select X / Y / Z axis
      UP  / DOWN   increase / decrease along selected axis
    """

    def __init__(self, ik_controller, fk_engine, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.listening = True
        self.ik = ik_controller
        self.fk = fk_engine
        self.active_axis = 2  # default: Z
        self._frame = 0

    def _tip(self):
        # d_endEffectorPose is Rigid3dTypes::Coord -> flat [x,y,z, qx,qy,qz,qw]
        pose = self.fk.d_endEffectorPose.value
        return pose[:3]  # position only

    def _print_status(self):
        p = list(self.ik.d_targetPosition.value)
        t = self._tip()
        err = ((p[0]-t[0])**2 + (p[1]-t[1])**2 + (p[2]-t[2])**2) ** 0.5
        print(f"[IK]  target: x={p[0]:7.2f}  y={p[1]:7.2f}  z={p[2]:7.2f} mm"
              f"   |   tip: x={t[0]:7.2f}  y={t[1]:7.2f}  z={t[2]:7.2f} mm"
              f"   |   err: {err:.2f} mm")

    def onAnimateEndEvent(self, _):
        self._frame += 1
        if self._frame % PRINT_EVERY == 0:
            self._print_status()

    def onKeypressedEvent(self, event):
        key = event['key']

        # --- axis selection ---
        if key == Key.KP_1 or key == '1':
            self.active_axis = 0
            print("[IK] Active axis: X")
            return
        if key == Key.KP_2 or key == '2':
            self.active_axis = 1
            print("[IK] Active axis: Y")
            return
        if key == Key.KP_3 or key == '3':
            self.active_axis = 2
            print("[IK] Active axis: Z")
            return

        # --- move target ---
        p = list(self.ik.d_targetPosition.value)

        if key == Key.uparrow:
            p[self.active_axis] += DP
        elif key == Key.downarrow:
            p[self.active_axis] -= DP
        else:
            return

        self.ik.d_targetPosition.value = p
        self._print_status()


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

    # IK controller: runs one NR step per animation frame
    ik = rootNode.addObject('TRMInverseKinematicsController',
                             name='ik',
                             d_targetPosition=[0, 0, 100])

    # FK engine: reads IK joint config -> end-effector pose for tip readout
    fk = rootNode.addObject('TRMForwardKinematicsEngine',
                             name='fk',
                             d_jointConfig=ik.d_jointConfig.getLinkPath())

    # Visual model reads joint config directly from the IK output
    rootNode.addObject('TRMVisualModel',
                       name='visual',
                       d_jointConfig=ik.d_jointConfig.getLinkPath(),
                       d_nSamples=30)

    rootNode.addObject(IKTargetController(ik_controller=ik, fk_engine=fk,
                                          name='targetController'))
