/**
 * @file TRMInverseKinematicsEngine.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Implement inverse kinematics based on Sofa data types
 * @version 0.1
 * @date 2026-05-04
 *
 * @copyright Copyright (c) 2026
 */

#include "TRMInverseKinematicsEngine.h"

namespace TRMCTR::engine
{

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
TRMInverseKinematicsEngine::TRMInverseKinematicsEngine()
    : d_targetPosition(initData(&d_targetPosition,
                                Vec3(0, 0, 0),
                                "d_targetPosition",
                                "Target end-effector position (mm)"))
    , d_jointConfig(initData(&d_jointConfig,
                             "d_jointConfig",
                             "Updated joint configuration after one IK step"))
    , m_jointConfig(0, 50, 0, 50, 0, 50)
{}

// ----------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------
void TRMInverseKinematicsEngine::init()
{
    addInput(&d_targetPosition);
    addOutput(&d_jointConfig);
    setDirtyValue();
}

// ----------------------------------------------------------------------------
// reinit
// ----------------------------------------------------------------------------
void TRMInverseKinematicsEngine::reinit()
{
    update();
}

// ----------------------------------------------------------------------------
// doUpdate
// ----------------------------------------------------------------------------
void TRMInverseKinematicsEngine::doUpdate()
{
    const Vec3& P = d_targetPosition.getValue();

    Eigen::Vector3d P_eigen(P[0], P[1], P[2]);

    Eigen::Matrix<double,6,1> q0_eigen;
    for (int i = 0; i < 6; ++i) q0_eigen(i) = m_jointConfig[i];

    // one Newton-Raphson step; result becomes the next frame's initial guess
    Eigen::Matrix<double,6,1> q_eigen = CTR::InverseKinematics::IK(P_eigen, q0_eigen);

    for (int i = 0; i < 6; ++i) m_jointConfig[i] = q_eigen(i);
    d_jointConfig.setValue(m_jointConfig);
}

} // namespace TRMCTR::engine
