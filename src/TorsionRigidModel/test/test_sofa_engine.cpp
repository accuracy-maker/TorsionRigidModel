/**
 * @file test_sofa_engine.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Test SOFA DataEngine wrappers for FK and IK:
 *        (1) Data field types are correct SOFA types
 *        (2) Computation produces physically valid results
 * @version 0.1
 * @date 2026-05-04
 *
 * @copyright Copyright (c) 2026
 */

#include <gtest/gtest.h>
#include <type_traits>
#include <Eigen/Dense>

#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>

#include "Engine/TRMForwardKinematicsEngine.h"
#include "Engine/TRMInverseKinematicsEngine.h"
#include "core/ForwardKinematics.h"
#include "core/RobotParameters.h"

// ----------------------------------------------------------------------------
// Expose protected constructors for testing without modifying production code
// ----------------------------------------------------------------------------
struct TestFKEngine : TRMCTR::engine::TRMForwardKinematicsEngine {};
struct TestIKEngine : TRMCTR::engine::TRMInverseKinematicsEngine {};

static constexpr double EPS = 1e-8;

// ============================================================================
// FK Engine — type checks
// ============================================================================
class FKEngineTypeTest : public ::testing::Test {
protected:
    using FKEngine = TRMCTR::engine::TRMForwardKinematicsEngine;
    void SetUp() override {}
};

// d_jointConfig must be Data<Vec<6,double>>
TEST_F(FKEngineTypeTest, JointConfigIsDataVec6)
{
    constexpr bool ok = std::is_same_v<
        std::remove_reference_t<decltype(std::declval<FKEngine>().d_jointConfig)>,
        sofa::core::objectmodel::Data<sofa::type::Vec<6, double>>
    >;
    EXPECT_TRUE(ok) << "d_jointConfig should be Data<Vec<6,double>>";
}

// d_endEffectorPose must be Data<Rigid3dTypes::Coord>
TEST_F(FKEngineTypeTest, EndEffectorPoseIsDataRigid3Coord)
{
    constexpr bool ok = std::is_same_v<
        std::remove_reference_t<decltype(std::declval<FKEngine>().d_endEffectorPose)>,
        sofa::core::objectmodel::Data<sofa::defaulttype::Rigid3dTypes::Coord>
    >;
    EXPECT_TRUE(ok) << "d_endEffectorPose should be Data<Rigid3dTypes::Coord>";
}

// ============================================================================
// FK Engine — computation checks
// ============================================================================
class FKEngineComputeTest : public ::testing::Test {
protected:
    TestFKEngine engine;
    void SetUp() override { engine.init(); }
};

// output quaternion must have unit norm
TEST_F(FKEngineComputeTest, OutputQuaternionIsUnit)
{
    engine.d_jointConfig.setValue(sofa::type::Vec<6,double>(0, 50, 0, 50, 0, 50));
    engine.update();

    const auto& quat = engine.d_endEffectorPose.getValue().getOrientation();
    double norm2 = quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3];
    EXPECT_NEAR(norm2, 1.0, EPS) << "Output quaternion must be unit norm";
}

// output position must be finite (no NaN / Inf)
TEST_F(FKEngineComputeTest, OutputPositionIsFinite)
{
    engine.d_jointConfig.setValue(sofa::type::Vec<6,double>(0, 50, 0, 50, 0, 50));
    engine.update();

    const auto& pos = engine.d_endEffectorPose.getValue().getCenter();
    EXPECT_TRUE(std::isfinite(pos[0]));
    EXPECT_TRUE(std::isfinite(pos[1]));
    EXPECT_TRUE(std::isfinite(pos[2]));

    std::cout << "FK output position: ["
              << pos[0] << ", " << pos[1] << ", " << pos[2] << "] mm\n";
}

// engine output must match direct CTR::ForwardKinematics::FK call
TEST_F(FKEngineComputeTest, OutputMatchesDirectFKCall)
{
    sofa::type::Vec<6,double> q_sofa(0, 50, 0, 50, 0, 50);
    engine.d_jointConfig.setValue(q_sofa);
    engine.update();

    // reference result from core FK
    Eigen::Matrix<double,6,1> q_eigen;
    for (int i = 0; i < 6; ++i) q_eigen(i) = q_sofa[i];
    Eigen::Matrix4d g = CTR::ForwardKinematics::FK(q_eigen);

    const auto& pos = engine.d_endEffectorPose.getValue().getCenter();
    EXPECT_NEAR(pos[0], g(0,3), EPS);
    EXPECT_NEAR(pos[1], g(1,3), EPS);
    EXPECT_NEAR(pos[2], g(2,3), EPS);
}

// ============================================================================
// IK Engine — type checks
// ============================================================================
class IKEngineTypeTest : public ::testing::Test {
protected:
    using IKEngine = TRMCTR::engine::TRMInverseKinematicsEngine;
    void SetUp() override {}
};

// d_targetPosition must be Data<Vec3d>
TEST_F(IKEngineTypeTest, TargetPositionIsDataVec3)
{
    constexpr bool ok = std::is_same_v<
        std::remove_reference_t<decltype(std::declval<IKEngine>().d_targetPosition)>,
        sofa::core::objectmodel::Data<sofa::type::Vec3d>
    >;
    EXPECT_TRUE(ok) << "d_targetPosition should be Data<Vec3d>";
}

// d_currentJointConfig must be Data<Vec<6,double>>
TEST_F(IKEngineTypeTest, CurrentJointConfigIsDataVec6)
{
    constexpr bool ok = std::is_same_v<
        std::remove_reference_t<decltype(std::declval<IKEngine>().d_currentJointConfig)>,
        sofa::core::objectmodel::Data<sofa::type::Vec<6, double>>
    >;
    EXPECT_TRUE(ok) << "d_currentJointConfig should be Data<Vec<6,double>>";
}

// d_jointConfig (output) must be Data<Vec<6,double>>
TEST_F(IKEngineTypeTest, OutputJointConfigIsDataVec6)
{
    constexpr bool ok = std::is_same_v<
        std::remove_reference_t<decltype(std::declval<IKEngine>().d_jointConfig)>,
        sofa::core::objectmodel::Data<sofa::type::Vec<6, double>>
    >;
    EXPECT_TRUE(ok) << "d_jointConfig output should be Data<Vec<6,double>>";
}

// ============================================================================
// IK Engine — computation checks
// ============================================================================
class IKEngineComputeTest : public ::testing::Test {
protected:
    TestIKEngine engine;
    void SetUp() override { engine.init(); }
};

// output joint lengths must stay within robot joint limits
TEST_F(IKEngineComputeTest, OutputJointConfigWithinLimits)
{
    engine.d_targetPosition.setValue(sofa::type::Vec3d(0, 0, 120));
    engine.update();

    const auto& q = engine.d_jointConfig.getValue();
    EXPECT_GE(q[1], CTR::RobotParameters::S1_MIN);
    EXPECT_LE(q[1], CTR::RobotParameters::S1_MAX);
    EXPECT_GE(q[3], CTR::RobotParameters::S2_MIN);
    EXPECT_LE(q[3], CTR::RobotParameters::S2_MAX);
    EXPECT_GE(q[5], CTR::RobotParameters::S3_MIN);
    EXPECT_LE(q[5], CTR::RobotParameters::S3_MAX);
}

// one IK step must reduce the position error toward the target
TEST_F(IKEngineComputeTest, OneIKStepReducesPositionError)
{
    sofa::type::Vec3d       target(0, 0, 120);
    sofa::type::Vec<6,double> q0(0, 50, 0, 50, 0, 50);

    engine.d_targetPosition.setValue(target);
    engine.d_currentJointConfig.setValue(q0);
    engine.update();

    const auto& q_new = engine.d_jointConfig.getValue();

    // run FK on initial and updated configs
    auto toEigen = [](const sofa::type::Vec<6,double>& v) {
        Eigen::Matrix<double,6,1> e;
        for (int i = 0; i < 6; ++i) e(i) = v[i];
        return e;
    };

    Eigen::Vector3d P_target(target[0], target[1], target[2]);
    Eigen::Vector3d pos0 = CTR::ForwardKinematics::FK(toEigen(q0)).block<3,1>(0,3);
    Eigen::Vector3d pos1 = CTR::ForwardKinematics::FK(toEigen(q_new)).block<3,1>(0,3);

    double err0 = (P_target - pos0).norm();
    double err1 = (P_target - pos1).norm();

    std::cout << "Position error before IK step: " << err0 << " mm\n";
    std::cout << "Position error after  IK step: " << err1 << " mm\n";

    EXPECT_LT(err1, err0) << "One IK step should reduce position error";
}
