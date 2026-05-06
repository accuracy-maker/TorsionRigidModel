# TorsionRigidModel Dev Log

---

## 2026-05-07 09:25

**Issue:** `d_currentJointConfig` was declared as a `Data<Vec6>` input, but for teleoperation the joint config needs to persist as internal state between frames — it should not be an external input.

**Solution:** Removed `d_currentJointConfig` as a `Data` field. Replaced with private member `Vec6 m_jointConfig` (init: `0,50,0,50,0,50`). `doUpdate()` now reads `m_jointConfig` as the initial guess, runs one NR step, and writes the result back — state carries over automatically each frame.

---

## 2026-05-07 09:30

**Issue:** IK was a `DataEngine`, which is lazy/reactive — `doUpdate` only fires when inputs go dirty. For teleoperation the IK must run every frame unconditionally regardless of whether the target moved. `DataEngine` is the wrong abstraction for stateful per-frame computation.

**Solution:** Refactored `TRMInverseKinematicsEngine` → `TRMInverseKinematicsController` (namespace `TRMCTR::controller`), now inheriting `sofa::core::objectmodel::BaseObject`. `init()` sets `f_listening = true`; `handleEvent()` checks for `AnimateBeginEvent` and runs one NR step per frame. Registered in `initTorsionRigidModel.cpp`. FK stays as `DataEngine` (reactive transform, correct fit).

---

## 2026-05-07 09:31

**Issue:** Build error — `sofa/simulation/AnimateBeginEvent.h: No such file or directory`. `Sofa.Simulation.Core` was not listed as a dependency in `CMakeLists.txt`.

**Solution:** Added `find_package(Sofa.Simulation.Core REQUIRED)` and `Sofa.Simulation.Core` to `target_link_libraries` in `CMakeLists.txt`.

---

## 2026-05-07 09:37

**Issue:** `TRMInverseKinematicsController` not found in SOFA factory at runtime. Used deprecated `sofa::core::RegisterObject` API (removed since SOFA v24.12), causing the registration to silently fail and no components to load from the plugin.

**Solution:** Replaced with `factory->registerObjects(sofa::core::ObjectRegistrationData(...).add<T>())`, matching the pattern used by `TRMForwardKinematicsEngine`. Also moved `#include <sofa/core/ObjectFactory.h>` from the header to the `.cpp`.

---

## 2026-05-07 09:38

**Issue:** Build error — `sofa::core::ObjectFactory` not found in header. The registration function declaration `void registerTRMInverseKinematicsController(sofa::core::ObjectFactory*)` was left in the header after `ObjectFactory.h` was moved to the `.cpp`.

**Solution:** Removed the declaration from the header. `initTorsionRigidModel.cpp` already provides its own `extern` declaration, so the header declaration is not needed.

---

## 2026-05-07 09:43

**Issue:** `IndexError: invalid index to scalar variable` when reading tip position. Assumed `d_endEffectorPose.value` returns `[[x,y,z], [qx,qy,qz,qw]]` but SofaPython3 returns a flat array `[x, y, z, qx, qy, qz, qw]`, so `pose[0]` was a scalar.

**Solution:** Changed `pose[0]` to `pose[:3]` to slice the position from the flat array.
