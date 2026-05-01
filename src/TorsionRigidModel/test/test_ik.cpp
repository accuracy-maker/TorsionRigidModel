/**
 * @file test_ik.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Testing for InverseKinematics
 * @version 0.1
 * @date 2026-05-01
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include <gtest/gtest.h>
#include "core/InverseKinematics.h"
#include "core/ForwardKinematics.h"
#include <Eigen/Dense>

const double EPS = 1e-8;

class InverseKinematicsTest : public ::testing::Test {
protected:
    void SetUp() override {

    }
};

// 1. test Jacobian Matrix
TEST_F(InverseKinematicsTest, JacobianIsValid) {
    Eigen::Matrix<double,6,1> q(0,1,0,1,0,1);

    Eigen::Matrix<double,6,6> J = CTR::InverseKinematics::Jacobian(q);

    // check the shape of the matrix
    EXPECT_EQ(J.rows(), 6);
    EXPECT_EQ(J.cols(), 6);

    // output the matrix
    std::cout << "Input the joint control is:\n" << q << std::endl;
    std::cout << "Output Jacobian Matrix is:\n" << J << std::endl;

}

// 2. test Inverse Kinematics
TEST_F(InverseKinematicsTest, InverseKinematicsConverges) {
    Eigen::Matrix<double,6,1> q(0,50,0,50,0,50);

    // Target end-effector position
    Eigen::Vector3d P(0, 0, 120);

    Eigen::Matrix<double,6,1> q_new = CTR::InverseKinematics::IK(P, q);

    Eigen::Matrix4d gn = CTR::ForwardKinematics::FK(q_new);

    Eigen::Vector3d P_c = gn.block<3,1>(0,3);

    std::cout << "Old Joint configurations are: \n" << q << std::endl;
    std::cout << "New Joint configurations are: \n" << q_new << std::endl;
    std::cout << "Current Location is:\n" << P_c << std::endl;
}