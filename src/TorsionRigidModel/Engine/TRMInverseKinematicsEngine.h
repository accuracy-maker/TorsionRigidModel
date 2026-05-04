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

/**
 * SOFA DataEngine wrapping CTR::InverseKinematics::IK().
 *
 * Runs one Newton-Raphson step per update. At each time step d_currentJointConfig
 * should be wired from the previous output to serve as the initial guess.
 *
 * Inputs
 *   d_targetPosition      target end-effector position (x, y, z)  [mm]
 *   d_currentJointConfig  current joint config used as initial guess [θ1,s1,θ2,s2,θ3,s3]
 *
 * Output
 *   d_jointConfig         updated joint configuration after one IK step
 */
class TRMInverseKinematicsEngine : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(TRMInverseKinematicsEngine, sofa::core::DataEngine);

    using Vec3 = sofa::type::Vec3d;
    using Vec6 = sofa::type::Vec<6, double>;

    // ------------------------------------------------------------------ Data
    Data<Vec3> d_targetPosition;     ///< input:  target end-effector position (mm)
    Data<Vec6> d_currentJointConfig; ///< input:  current joint config as initial guess
    Data<Vec6> d_jointConfig;        ///< output: updated joint config after one IK step

    // -------------------------------------------------------- SOFA life-cycle
    void init()     override;
    void reinit()   override;
    void doUpdate() override;

protected:
    TRMInverseKinematicsEngine();
    ~TRMInverseKinematicsEngine() override = default;
};

} // namespace TRMCTR::engine
