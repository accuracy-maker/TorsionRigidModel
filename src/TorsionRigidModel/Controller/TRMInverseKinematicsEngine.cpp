/**
 * @file TRMInverseKinematicsEngine.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Per-frame IK controller for the Torsion Rigid Model.
 * @version 0.2
 * @date 2026-05-07
 *
 * @copyright Copyright (c) 2026
 */

#include "TRMInverseKinematicsEngine.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/AnimateBeginEvent.h>

namespace TRMCTR::controller
{

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
TRMInverseKinematicsController::TRMInverseKinematicsController()
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
void TRMInverseKinematicsController::init()
{
    f_listening.setValue(true);
}

// ----------------------------------------------------------------------------
// handleEvent — fires every animation frame
// ----------------------------------------------------------------------------
void TRMInverseKinematicsController::handleEvent(sofa::core::objectmodel::Event* event)
{
    if (!sofa::simulation::AnimateBeginEvent::checkEventType(event))
        return;

    const Vec3& P = d_targetPosition.getValue();

    Eigen::Vector3d P_eigen(P[0], P[1], P[2]);

    Eigen::Matrix<double,6,1> q0_eigen;
    for (int i = 0; i < 6; ++i) q0_eigen(i) = m_jointConfig[i];

    // one Newton-Raphson step; result becomes next frame's initial guess
    Eigen::Matrix<double,6,1> q_eigen = CTR::InverseKinematics::IK(P_eigen, q0_eigen);

    for (int i = 0; i < 6; ++i) m_jointConfig[i] = q_eigen(i);
    d_jointConfig.setValue(m_jointConfig);
}

// ----------------------------------------------------------------------------
// Registration
// ----------------------------------------------------------------------------
void registerTRMInverseKinematicsController(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
        "CTR inverse kinematics controller — one NR step per frame.")
        .add<TRMInverseKinematicsController>());
}

} // namespace TRMCTR::controller
