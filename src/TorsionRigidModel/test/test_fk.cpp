/**
 * @file test_fk.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Testing for ForwardKinematics
 * @version 0.1
 * @date 2026-04-30
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include <gtest/gtest.h>
#include "core/ForwardKinematics.h"
#include <Eigen/Dense>

const double EPS = 1e-8;

/**
 * @brief Test for ForwardKinematics
*/

class ForwardKinematicsTest : public ::testing::Test {
protected:
    void SetUp() override {

    }
};

// 1. check if output is a valid SE(3)
TEST_F(ForwardKinematicsTest, OutputIsValidSE3) {
    Eigen::Matrix<double,6,1> q;
    q << 0, 1, 0, 1, 0, 1;
    Eigen::Matrix4d g = CTR::ForwardKinematics::FK(q);

    // check if R^T R = I
    Eigen::Matrix3d R = g.block<3,3>(0,0);
    EXPECT_NEAR((R * R.transpose() - Eigen::Matrix3d::Identity()).norm(), 0.0, EPS);
    EXPECT_NEAR(R.determinant(), 1.0, EPS);

    // check if the last row is [0,0,0,1]
    EXPECT_EQ(g.row(3), Eigen::RowVector4d(0,0,0,1));

    // output g
    std::cout << "Output SE(3) 4x4 matrix g is:\n" << g << std::endl;
}

// 2. test backbone sampling
TEST_F(ForwardKinematicsTest, BackboneSize) {
    Eigen::Matrix<double,6,1> q;
    q << 0, 50, 0, 50, 0, 50;
    int nSamples = 10;

    std::vector<Eigen::Vector3d> pts = CTR::ForwardKinematics::backbone(q, nSamples);

    EXPECT_EQ((int)pts.size(), 3 * nSamples + 1);
}

TEST_F(ForwardKinematicsTest, BackboneFirstPointIsOrigin) {
    Eigen::Matrix<double,6,1> q;
    q << 0, 50, 0, 50, 0, 50;

    std::vector<Eigen::Vector3d> pts = CTR::ForwardKinematics::backbone(q, 10);

    EXPECT_NEAR(pts.front().norm(), 0.0, EPS);
}

TEST_F(ForwardKinematicsTest, BackboneLastPointMatchesFKTip) {
    Eigen::Matrix<double,6,1> q;
    q << 0, 50, 0, 50, 0, 50;

    std::vector<Eigen::Vector3d> pts = CTR::ForwardKinematics::backbone(q, 20);
    Eigen::Vector3d tip_backbone = pts.back();
    Eigen::Vector3d tip_fk = CTR::ForwardKinematics::FK(q).block<3,1>(0,3);

    EXPECT_NEAR((tip_backbone - tip_fk).norm(), 0.0, EPS);
}

TEST_F(ForwardKinematicsTest, BackboneAllPointsFinite) {
    Eigen::Matrix<double,6,1> q;
    q << 0, 50, 0, 50, 0, 50;

    std::vector<Eigen::Vector3d> pts = CTR::ForwardKinematics::backbone(q, 10);

    for (const auto& p : pts)
        EXPECT_TRUE(p.allFinite());
}

TEST_F(ForwardKinematicsTest, BackbonePrint) {
    Eigen::Matrix<double,6,1> q;
    q << 0, 50, 0, 50, 0, 50;

    std::vector<Eigen::Vector3d> pts = CTR::ForwardKinematics::backbone(q, 10);

    for (size_t i = 0; i < pts.size(); ++i)
        std::cout << "pt[" << i << "]: " << pts[i].transpose() << "\n";
}