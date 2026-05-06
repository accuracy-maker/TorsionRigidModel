import Sofa.Core
from Sofa.constants import Key
import math

# Joint limits matching RobotParameters.h
S_MIN   = [10.0, 10.0, 25.0]
S_MAX   = [100.0, 100.0, 100.0]
THETA_MIN = -math.pi
THETA_MAX =  math.pi

DS     = 5.0          # mm per keypress
DTHETA = 0.1          # rad per keypress (~5.7 deg)


class CTRKeyboardController(Sofa.Core.Controller):
    """
    Keyboard control for 3-tube CTR forward kinematics.

    Key layout:
      1 / 2 / 3       select active tube
      UP  / DOWN       extend / retract arc length s
      LEFT / RIGHT     rotate theta CCW / CW
    """

    def __init__(self, fk_engine, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.listening = True
        self.fk = fk_engine
        self.active_tube = 0           # 0-indexed: 0=tube1, 1=tube2, 2=tube3

    def onKeypressedEvent(self, event):
        key = event['key']

        # --- tube selection ---
        if key == Key.KP_1 or key == '1':
            self.active_tube = 0
            print(f"[CTR] Active tube: 1")
            return
        if key == Key.KP_2 or key == '2':
            self.active_tube = 1
            print(f"[CTR] Active tube: 2")
            return
        if key == Key.KP_3 or key == '3':
            self.active_tube = 2
            print(f"[CTR] Active tube: 3")
            return

        # --- read current config ---
        # q = [theta1, s1, theta2, s2, theta3, s3]
        q = list(self.fk.d_jointConfig.value)
        i = self.active_tube
        theta_idx = i * 2
        s_idx     = i * 2 + 1

        changed = False

        if key == Key.uparrow:
            q[s_idx] = min(q[s_idx] + DS, S_MAX[i])
            changed = True
        elif key == Key.downarrow:
            q[s_idx] = max(q[s_idx] - DS, S_MIN[i])
            changed = True
        elif key == Key.leftarrow:
            q[theta_idx] = max(q[theta_idx] - DTHETA, THETA_MIN)
            changed = True
        elif key == Key.rightarrow:
            q[theta_idx] = min(q[theta_idx] + DTHETA, THETA_MAX)
            changed = True

        if changed:
            self.fk.d_jointConfig.value = q
            print(f"[CTR] tube={i+1}  theta={q[theta_idx]:.3f} rad  s={q[s_idx]:.1f} mm")


def createScene(rootNode):
    rootNode.dt = 0.01
    rootNode.gravity = [0, 0, 0]

    rootNode.addObject('RequiredPlugin', name='TorsionRigidModel')
    rootNode.addObject('RequiredPlugin', name='Sofa.Component.Visual')
    rootNode.addObject('RequiredPlugin', name='Sofa.Component.AnimationLoop')

    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels')

    fk = rootNode.addObject('TRMForwardKinematicsEngine',
                             name='fk',
                             d_jointConfig=[0, 5, 0, 5, 0, 5])

    rootNode.addObject('TRMVisualModel',
                       name='visual',
                       d_jointConfig=fk.d_jointConfig.getLinkPath(),
                       d_nSamples=30)

    rootNode.addObject(CTRKeyboardController(fk_engine=fk, name='kbController'))
