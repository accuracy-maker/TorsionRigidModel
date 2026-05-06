/**
 * @file TRMInverseKinematicsEngine.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Data Engine of Sofa Framework based on the implementation of IK of Torsion Rigid model.
 * @version 0.1
 * @date 2026-05-01
 *
 * @copyright Copyright (c) 2026
 */

#pragma once

#include <sofa/core/DataEngine.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>

#include "core/InverseKinematics.h"

namespace TRMCTR::engine
{

using sofa::core::objectmodel::Data;

/**
 * SOFA DataEngine wrapping CTR::InverseKinematics::IK().
 *
 * Runs one Newton-Raphson step per SOFA frame. The joint config is kept as
 * internal state (m_jointConfig) so each frame's result feeds the next as the
 * initial guess — suitable for real-time teleoperation.
 *
 * Input
 *   d_targetPosition  target end-effector position (x, y, z)  [mm]
 *
 * Output
 *   d_jointConfig     updated joint configuration after one IK step
 */
class TRMInverseKinematicsEngine : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(TRMInverseKinematicsEngine, sofa::core::DataEngine);

    using Vec3 = sofa::type::Vec3d;
    using Vec6 = sofa::type::Vec<6, double>;

    // ------------------------------------------------------------------ Data
    Data<Vec3> d_targetPosition; ///< input:  target end-effector position (mm)
    Data<Vec6> d_jointConfig;    ///< output: updated joint configuration after one IK step

    // -------------------------------------------------------- SOFA life-cycle
    void init()     override;
    void reinit()   override;
    void doUpdate() override;

protected:
    TRMInverseKinematicsEngine();
    ~TRMInverseKinematicsEngine() override = default;

private:
    Vec6 m_jointConfig; ///< internal state: carries the IK solution across frames
};

} // namespace TRMCTR::engine
