# TorsionRigidModel

A SOFA plugin implementing forward and inverse kinematics for a **3-tube Concentric Tube Robot (CTR)** using Lie group mathematics.

## Overview

Concentric Tube Robots are slender continuum robots made of nested pre-curved elastic tubes. Each tube can rotate and translate independently, giving the robot 6 degrees of freedom. By exploiting the tubes' elastic interactions, the robot can navigate complex curved anatomical paths — making CTRs well-suited for minimally invasive surgery.

This plugin provides:

- **Lie group operations** on SE(3) and se(3): hat/vee operators, exponential map, adjoint transforms, and the A-matrix
- **Forward kinematics (FK)** via the Product of Matrix Exponentials (PoME) method — wrapped as a SOFA `DataEngine`
- **Inverse kinematics (IK)** via a Jacobian-based Newton-Raphson solver — wrapped as a per-frame SOFA `BaseObject` controller
- **Visual model** rendering the CTR backbone as a coloured polyline
- **Unit tests** for all core components using GoogleTest

---

## Code Structure

```
TorsionRigidModel/
├── CMakeLists.txt
├── LICENSE
└── src/TorsionRigidModel/
    ├── initTorsionRigidModel.cpp       plugin entry point
    ├── math/
    │   └── Lie.h                       Lie group / algebra operations (SE(3), se(3))
    ├── core/
    │   ├── RobotParameters.h           physical constants: stiffness, curvatures, joint limits
    │   ├── ForwardKinematics.h/cpp     FK algorithm (PoME)
    │   └── InverseKinematics.h/cpp     IK algorithm (Newton-Raphson)
    ├── Engine/
    │   └── TRMForwardKinematicsEngine.h/cpp    SOFA DataEngine wrapping FK
    ├── Controller/
    │   └── TRMInverseKinematicsEngine.h/cpp    SOFA BaseObject controller wrapping IK
    ├── Visual/
    │   └── TRMVisualModel.h/cpp        SOFA VisualModel rendering the CTR backbone
    └── test/
        ├── test_lie.cpp                Lie group unit tests
        ├── test_fk.cpp                 FK validity tests (SE(3) membership)
        ├── test_ik.cpp                 Jacobian shape + IK convergence tests
        └── test_sofa_engine.cpp        SOFA integration tests

scenes/
├── CTR_visual_fk_test.py              interactive FK demo (keyboard control)
└── CTR_visual_ik_test.py              interactive IK demo (target teleoperation)
```

---

## SOFA Components

### `TRMForwardKinematicsEngine` — FK DataEngine

| | |
|---|---|
| **Base class** | `sofa::core::DataEngine` |
| **Namespace** | `TRMCTR::engine` |
| **Trigger** | Lazy — recomputes when `d_jointConfig` changes |

| Field | Direction | Type | Description |
|---|---|---|---|
| `d_jointConfig` | input | `Vec<6, double>` | `[θ₁, s₁, θ₂, s₂, θ₃, s₃]` (rad, mm) |
| `d_endEffectorPose` | output | `Rigid3Coord` | SE(3) tip pose (position + quaternion) |

---

### `TRMInverseKinematicsController` — IK Controller

| | |
|---|---|
| **Base class** | `sofa::core::objectmodel::BaseObject` |
| **Namespace** | `TRMCTR::controller` |
| **Trigger** | `AnimateBeginEvent` — one Newton-Raphson step per animation frame |

| Field | Direction | Type | Description |
|---|---|---|---|
| `d_targetPosition` | input | `Vec3d` | target end-effector position (mm) |
| `d_jointConfig` | output | `Vec<6, double>` | updated joint config after one IK step |

Internal state `m_jointConfig` (init: `[0, 50, 0, 50, 0, 50]`) carries the solution across frames, enabling iterative convergence during teleoperation.

---

### `TRMVisualModel` — Backbone Visual Model

| | |
|---|---|
| **Base class** | `sofa::core::visual::VisualModel` |
| **Namespace** | `TRMCTR::visual` |

| Field | Direction | Type | Description |
|---|---|---|---|
| `d_jointConfig` | input | `Vec<6, double>` | joint config (link from FK engine or IK controller) |
| `d_nSamples` | input | `int` | samples per section (default 20) |

Renders the backbone as a coloured polyline: **blue** (section 1), **green** (section 2), **red** (section 3).

---

## Scene Files

### FK Test — `CTR_visual_fk_test.py`

Direct joint-space control. Keyboard maps to joint config → FK engine → visual model.

```
CTRKeyboardController ──► d_jointConfig ──► TRMForwardKinematicsEngine ──► TRMVisualModel
```

| Key | Action |
|---|---|
| `1` / `2` / `3` | Select active tube |
| `↑` / `↓` | Extend / retract arc length s |
| `←` / `→` | Rotate θ CCW / CW |

### IK Test — `CTR_visual_ik_test.py`

Target teleoperation. Keyboard sets target position → IK controller converges one step per frame → FK computes tip for readout → visual model renders result.

```
IKTargetController ──► d_targetPosition ──► TRMInverseKinematicsController
                                                       │
                                               d_jointConfig
                                              ╱               ╲
                        TRMForwardKinematicsEngine         TRMVisualModel
                            (tip readout)                  (backbone render)
```

| Key | Action |
|---|---|
| `1` / `2` / `3` | Select X / Y / Z axis |
| `↑` / `↓` | Move target ±5 mm along selected axis |

Console prints target position, current tip position, and Euclidean error every 50 frames and on every keypress.

---

## Robot Parameters

Defined in `core/RobotParameters.h`:

| Parameter | Description |
|---|---|
| `K1, K2, K3` | Stiffness matrices per tube (3×3 diagonal, Nm²) |
| `U1F1, U2F2, U3F3` | Intrinsic curvature vectors per tube (1/mm) |
| `S{1,2,3}_{MIN,MAX}` | Arc length limits per tube (mm) |
| `THETA_{MIN,MAX}` | Rotation limits (±π rad) |

Joint configuration convention: **`q = [θ₁, s₁, θ₂, s₂, θ₃, s₃]`**

---

## Dependencies

| Dependency | Purpose |
|---|---|
| **SOFA Framework** (`Sofa.Core`, `Sofa.GL`) | simulation framework |
| **Sofa.Simulation.Core** | `AnimateBeginEvent` for per-frame IK |
| **Eigen3** | linear algebra |
| **GoogleTest** | unit testing |

---

## Build

This plugin is built as part of SOFA. From the SOFA root:

```bash
# Configure — enable the plugin
cmake -S src -B build -G Ninja -DPLUGIN_TORSIONRIGIDMODEL=ON

# Build
cmake --build build

# Run tests
ctest --test-dir build -V -R TorsionRigidModel
```

To run a single test binary:

```bash
./build/bin/TorsionRigidModel_test --gtest_filter=InverseKinematicsTest.InverseKinematicsConverges
```

---

## Author

Haitao Gao — haitao.gao@unsw.edu.au

## License

GNU General Public License v3 — see [LICENSE](LICENSE).
