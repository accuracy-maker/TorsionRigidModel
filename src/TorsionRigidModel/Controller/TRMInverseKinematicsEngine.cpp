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
    , d_currentJointConfig(initData(&d_currentJointConfig,
                                    Vec6(0, 50, 0, 50, 0, 50),
                                    "d_currentJointConfig",
                                    "Current joint config as initial guess [theta1,s1,theta2,s2,theta3,s3]"))
    , d_jointConfig(initData(&d_jointConfig,
                             "d_jointConfig",
                             "Updated joint configuration after one IK step"))
{}

// ----------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------
void TRMInverseKinematicsEngine::init()
{
    addInput(&d_targetPosition);
    addInput(&d_currentJointConfig);
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
    const Vec3& P  = d_targetPosition.getValue();
    const Vec6& q0 = d_currentJointConfig.getValue();

    // convert SOFA types to Eigen
    Eigen::Vector3d P_eigen(P[0], P[1], P[2]);

    Eigen::Matrix<double,6,1> q0_eigen;
    for (int i = 0; i < 6; ++i) q0_eigen(i) = q0[i];

    // one Newton-Raphson step
    Eigen::Matrix<double,6,1> q_eigen = CTR::InverseKinematics::IK(P_eigen, q0_eigen);

    // convert result back to SOFA type and write output
    Vec6 q;
    for (int i = 0; i < 6; ++i) q[i] = q_eigen(i);
    d_jointConfig.setValue(q);
}

} // namespace TRMCTR::engine
