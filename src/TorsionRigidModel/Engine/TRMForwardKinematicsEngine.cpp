/**
 * @file TRMForwardKinematicsEngine.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Implementation of FK based on Sofa functions
 * @version 0.1
 * @date 2026-05-01
 *
 * @copyright Copyright (c) 2026
 */

#include "TRMForwardKinematicsEngine.h"


namespace TRMCTR::engine
{

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
TRMForwardKinematicsEngine::TRMForwardKinematicsEngine()
    : d_jointConfig(initData(&d_jointConfig,
                             Vec6(0, 50, 0, 50, 0, 50),
                             "d_jointConfig",
                             "[theta1, s1, theta2, s2, theta3, s3] (rad, mm)"))
    , d_endEffectorPose(initData(&d_endEffectorPose,
                                 "d_endEffectorPose",
                                 "End-effector SE(3) pose as position + quaternion"))
{}

// ----------------------------------------------------------------------------
// init — register inputs / outputs with the DataEngine dependency graph
// ----------------------------------------------------------------------------
void TRMForwardKinematicsEngine::init()
{
    addInput(&d_jointConfig);
    addOutput(&d_endEffectorPose);
    setDirtyValue();  // force doUpdate() on first simulation step
}

// ----------------------------------------------------------------------------
// reinit — called when a Data is edited in the SOFA GUI
// ----------------------------------------------------------------------------
void TRMForwardKinematicsEngine::reinit()
{
    update();
}

// ----------------------------------------------------------------------------
// doUpdate — core computation, called automatically when d_jointConfig changes
// ----------------------------------------------------------------------------
void TRMForwardKinematicsEngine::doUpdate()
{
    // 1. read   const Vec6& q = d_jointConfig.getValue();
    const Vec6& q = d_jointConfig.getValue();

    Eigen::Matrix<double,6,1> q_eigen;
    for (int i = 0; i < 6; ++i) q_eigen(i) = q[i];

    Eigen::Matrix4d g = CTR::ForwardKinematics::FK(q_eigen);
    // 3. write  d_endEffectorPose.setValue(toRigid3(g));
    d_endEffectorPose.setValue(toRigid3(g));
}

// ----------------------------------------------------------------------------
// toRigid3 — convert 4×4 SE(3) Eigen matrix to SOFA Rigid3Coord
// ----------------------------------------------------------------------------
TRMForwardKinematicsEngine::Rigid3Coord
TRMForwardKinematicsEngine::toRigid3(const Eigen::Matrix4d& g)
{
    // 1. Extract translation:  sofa::type::Vec3d pos(T(0,3), T(1,3), T(2,3));
    Eigen::Vector3d p = g.block<3,1>(0,3);
    // 2. Extract rotation 3×3 block and convert to sofa::type::Quatd
    Eigen::Matrix3d R = g.block<3,3>(0,0);
    Eigen::Matrix<double,4,1> Q = Lie::Rot2Quat(R);

    // Lie::Rot2Quat returns (w, x, y, z); sofa::type::Quatd stores (x, y, z, w)
    sofa::type::Quatd quat(Q(1), Q(2), Q(3), Q(0));
    sofa::type::Vec3d pos(p(0), p(1), p(2));

    return Rigid3Coord(pos, quat);
}

} // namespace TRMCTR::engine
