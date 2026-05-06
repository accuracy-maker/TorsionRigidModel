/**
 * @file TRMInverseKinematicsEngine.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Per-frame IK controller for the Torsion Rigid Model.
 * @version 0.2
 * @date 2026-05-07
 *
 * @copyright Copyright (c) 2026
 */

#pragma once

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/type/Vec.h>

#include "core/InverseKinematics.h"

namespace TRMCTR::controller
{

using sofa::core::objectmodel::Data;

/**
 * SOFA BaseObject controller wrapping CTR::InverseKinematics::IK().
 *
 * Fires on AnimateBeginEvent — one Newton-Raphson step per animation frame.
 * Joint config is kept as internal state (m_jointConfig) so each frame's
 * result feeds the next as the initial guess, suitable for teleoperation.
 *
 * Input
 *   d_targetPosition  target end-effector position (x, y, z)  [mm]
 *
 * Output
 *   d_jointConfig     updated joint configuration after one IK step
 */
class TRMInverseKinematicsController : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(TRMInverseKinematicsController, sofa::core::objectmodel::BaseObject);

    using Vec3 = sofa::type::Vec3d;
    using Vec6 = sofa::type::Vec<6, double>;

    // ------------------------------------------------------------------ Data
    Data<Vec3> d_targetPosition; ///< input:  target end-effector position (mm)
    Data<Vec6> d_jointConfig;    ///< output: updated joint config after one IK step

    // -------------------------------------------------------- SOFA life-cycle
    void init() override;
    void handleEvent(sofa::core::objectmodel::Event* event) override;

protected:
    TRMInverseKinematicsController();
    ~TRMInverseKinematicsController() override = default;

private:
    Vec6 m_jointConfig; ///< internal state: carries the IK solution across frames
};

} // namespace TRMCTR::controller
