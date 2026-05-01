# TorsionRigidModel

A SOFA plugin implementing forward and inverse kinematics for a **3-tube Concentric Tube Robot (CTR)** using Lie group mathematics.

## Overview

Concentric Tube Robots are slender continuum robots made of nested pre-curved elastic tubes. Each tube can rotate and translate independently, giving the robot 6 degrees of freedom. By exploiting the tubes' elastic interactions, the robot can navigate complex curved anatomical paths — making CTRs well-suited for minimally invasive surgery.

This plugin provides the mathematical core for simulating a 3-tube CTR within the SOFA simulation framework:

- **Lie group operations** on SE(3) and se(3): hat/vee operators, exponential map, adjoint transforms, and the specialised A-matrix
- **Forward kinematics (FK)** via the Product of Matrix Exponentials (PoME) method
- **Inverse kinematics (IK)** via a Jacobian-based Newton-Raphson solver with joint-limit enforcement
- **Unit tests** for all three components using GoogleTest

## Code Structure

```
TorsionRigidModel/
├── CMakeLists.txt
├── LICENSE                         (GPL v3)
└── src/
    └── TorsionRigidModel/
        ├── math/
        │   └── Lie.h               Lie group / algebra operations (SE(3), se(3))
        ├── core/
        │   ├── RobotParameters.h   Physical constants: stiffness, curvatures, joint limits
        │   ├── ForwardKinematics.h
        │   ├── ForwardKinematics.cpp
        │   ├── InverseKinematics.h
        │   └── InverseKinematics.cpp
        ├── test/
        │   ├── test_lie.cpp        9 unit tests for Lie operations
        │   ├── test_fk.cpp         FK validity tests (SE(3) membership)
        │   └── test_ik.cpp         Jacobian shape + IK convergence tests
        ├── Constraint/             (placeholder — not yet implemented)
        ├── Engine/                 (placeholder — not yet implemented)
        ├── ForceField/             (placeholder — not yet implemented)
        └── Mapping/                (placeholder — not yet implemented)
```

### Key Components

#### `math/Lie.h` — Lie Group Mathematics

All operations live in the `Lie` namespace.

| Function | Description |
|---|---|
| `hat(Vector3d)` → `Matrix3d` | so(3) hat: vector → skew-symmetric matrix |
| `vee(Matrix3d)` → `Vector3d` | so(3) vee: inverse of hat |
| `hat(Matrix<6,1>)` → `Matrix4d` | se(3) hat: twist → 4×4 matrix |
| `vee(Matrix4d)` → `Matrix<6,1>` | se(3) vee: inverse of hat |
| `exp(Matrix<6,1>)` → `Matrix4d` | Exponential map se(3) → SE(3) |
| `Adg(Matrix4d)` → `Matrix<6,6>` | Adjoint representation of a group element |
| `aMatrix(Matrix<6,1>, double)` | Integration matrix for curved-section Jacobian |

#### `core/ForwardKinematics` — FK via Product of Matrix Exponentials

Given a joint configuration **q = [θ₁, s₁, θ₂, s₂, θ₃, s₃]**, computes the tip pose as a 4×4 SE(3) transformation:

```
g = exp(ξ₁·s₁) · exp(ξ₂·s₂) · exp(ξ₃·s₃)
```

Each twist ξᵢ encodes the resultant curvature of that tube section, calculated from the tubes' stiffness matrices **Kᵢ** and intrinsic curvatures **uᵢ** rotated by θᵢ.

#### `core/InverseKinematics` — Jacobian-based IK

Given a target position **P ∈ ℝ³** and an initial guess **q₀**, returns an updated joint configuration:

```
δq = W · Jᵀ · (J · W · Jᵀ + εI)⁻¹ · e
q  = q₀ + δq
```

where **e** is the position error and **W** is a diagonal weight matrix that emphasises translation DOFs. Joint limits (extension and rotation bounds from `RobotParameters`) are enforced by clamping.

#### `core/RobotParameters.h` — Physical Constants

Defines the static parameters for the 3-tube robot:

| Parameter | Description |
|---|---|
| `K1, K2, K3` | Stiffness matrices (100I, 10I, 0) |
| `U1F1, U2F2, U3F3` | Intrinsic curvature vectors per tube |
| `S{1,2,3}_{MAX,MIN}` | Tube extension limits (mm) |
| `THETA_{MAX,MIN}` | Rotation limits (±π) |

## Dependencies

- **C++17**
- **Eigen3** — linear algebra
- **GoogleTest** — unit testing

On Ubuntu/Debian:

```bash
sudo apt install libeigen3-dev libgtest-dev cmake ninja-build
```

## Build and Test

Build output goes into a `build/` directory inside the plugin folder.

```bash
# 1. Create the build directory
mkdir -p /home/xxx/sofa/src/applications/plugins/TorsionRigidModel/build
cd /home/xxx/sofa/src/applications/plugins/TorsionRigidModel/build

# 2. Configure with CMake
cmake ..

# 3. Compile
make

# 4. Run all tests via CTest
ctest --output-on-failure
```

Or run each test binary directly:

```bash
# Lie group / algebra tests
./test_lie

# Forward kinematics tests
./test_fk

# Inverse kinematics tests
./test_ik
```

To run a single named test:

```bash
./test_ik --gtest_filter=InverseKinematicsTest.InverseKinematicsConverges
```

## Author

Haitao Gao — haitao.gao@unsw.edu.au

## License

GNU General Public License v3 — see [LICENSE](LICENSE).
