/**
 * @file test_sofa_engine.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Test SOFA DataEngine wrappers for FK:
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
#include "core/ForwardKinematics.h"
#include "core/RobotParameters.h"

// ----------------------------------------------------------------------------
// Expose protected constructors for testing without modifying production code
// ----------------------------------------------------------------------------
struct TestFKEngine : TRMCTR::engine::TRMForwardKinematicsEngine {};

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

