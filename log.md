# TorsionRigidModel Dev Log

---

## 2026-05-07 09:25

**Issue:** `d_currentJointConfig` was declared as a `Data<Vec6>` input, but for teleoperation the joint config needs to persist as internal state between frames — it should not be an external input.

**Solution:** Removed `d_currentJointConfig` as a `Data` field. Replaced with private member `Vec6 m_jointConfig` (init: `0,50,0,50,0,50`). `doUpdate()` now reads `m_jointConfig` as the initial guess, runs one NR step, and writes the result back — state carries over automatically each frame.
