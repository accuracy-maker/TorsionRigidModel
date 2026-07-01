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

---

## 2026-07-01

**Task:** Build a (configuration space, workspace) dataset for the CTR via Sobol sampling, and analyse it for singularities and multiple-solution ambiguity.

**Solution:** `scenes/build_dataset.py` drives a headless `TRMForwardKinematicsEngine` node (no `runSofa`/animation loop needed — reads trigger the `DataEngine`'s lazy recompute) over `scipy.stats.qmc.Sobol` samples of the 6-D joint box, using the true per-tube limits from `RobotParameters.h` (S1/S2 ∈ [10,100], S3 ∈ [25,100], θ ∈ [-π,π]). Saves `q`, tip position `x`, full pose, and `manipulability` to `dataset.npz`. `scenes/visulise_dataset.py` plots the 3-D workspace scatter (coloured by manipulability), per-joint histograms, and a manipulability histogram.

---

**Issue:** Wanted a `sqrt(det(J·Jᵀ))` manipulability/singularity measure per sample, but `CTR::InverseKinematics::Jacobian()` was never exposed to Python — only `TRMForwardKinematicsEngine::d_endEffectorPose` was.

**Solution:** Added `d_manipulability` output to `TRMForwardKinematicsEngine`. Using the full 6×6 twist Jacobian directly produced `NaN` for most samples: `theta1` **and** `theta3` are structurally dead DOFs given the current `RobotParameters.h` (`U1F1 = U3F3 = (0,0,0)`), so two columns of the 6×6 Jacobian are always exactly zero and `det(J·Jᵀ) ≡ 0`, with float round-off randomly flipping its sign before `sqrt()`. Fixed by extracting a new `CTR::InverseKinematics::PositionJacobian()` (the 3×6 translational sub-Jacobian `[-hat(p), I₃]·J(q)`, already computed inline inside `IK()`) and using `sqrt(max(0, det(Jp·Jpᵀ)))` instead — generically full-rank, no more NaNs. `IK()` was refactored to call the same new method instead of duplicating the formula; its solve behaviour (position-only, translation error, no orientation term) was not changed.

---

**Task:** Distinguish "multiple q map to the same x" cases caused by genuine kinematic multiplicity vs. proximity to a singularity.

**Solution:** `scenes/singularity_check.py` builds a KDTree over workspace positions `x`, finds each sample's nearest neighbour, and computes the corresponding joint-space distance (5-D, excluding the dead `theta1`, circular metric for angles). Pairs with a small `dist_x` but disproportionately large `dist_q` are anomalous; whether the pair's `manipulability` is low (near-singularity) or not (genuine multiplicity) disambiguates the cause. On the 100k dataset: ~318 distinct genuine-multiplicity pairs (real alternate solution branches, manipulability well above threshold) vs. ~182 distinct near-singularity pairs (explained by low manipulability). Noted follow-up: condition number (`σ_max/σ_min`) would be a more targeted classifier than raw manipulability here, since it isolates rank-loss specifically rather than conflating it with overall ellipsoid scale.

---

**Task:** Visualise two joint-space solutions side by side in SOFA to sanity-check the multiplicity/near-singularity findings.

**Solution:** `TRMVisualModel` had no way to recolour or fade a CTR instance — colours were hardcoded per tube section. Added `d_useFlatColor` + `d_color` (RGBA) Data fields; transparency needed no extra plumbing since `DrawToolGL::setMaterial()` already auto-enables GL blending when `color[3] < 1`. `scenes/visualise_multiplicity.py` overlays a solid CTR (`q1`) with a transparent orange "ghost" CTR (`q2`), switchable between the `"multiplicity"` and `"near_singularity"` example pairs found by `singularity_check.py`.
