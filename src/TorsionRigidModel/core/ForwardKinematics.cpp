/**
 * @file ForwardKinematics.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Forward Kinematics implementation of Product of Matrix Exponential
 * @version 0.1
 * @date 2026-04-30
 * 
 * @copyright Copyright (c) 2026
 * 
 */
#include "ForwardKinematics.h"
#include <cmath>
#include <iostream>

namespace CTR {

Eigen::Matrix3d ForwardKinematics::rotz(const double theta)
{
    Eigen::Matrix3d R;
    R << cos(theta), -sin(theta), 0,
         sin(theta),  cos(theta), 0,
         0,           0,          1;
    return R;
}

Eigen::Matrix4d ForwardKinematics::FK(const Eigen::Matrix<double,6,1>& q){
    // 1. extract variables
    double theta1 = q(0), s1 = q(1), theta2 = q(2), s2 = q(3), theta3 = q(4), s3 = q(5);


    // 2. compute the resultant curvature in the world frame (eq 6 in Leo's paper)
    Eigen::Matrix3d K123 = RobotParameters::K1 + RobotParameters::K2 + RobotParameters::K3;
    Eigen::Matrix3d K23 = RobotParameters::K2 + RobotParameters::K3;

    Eigen::Matrix3d K123inv = K123.inverse();
    Eigen::Matrix3d K23inv = K23.inverse();

    Eigen::Vector3d ufW1 = K123inv*(RobotParameters::K1*rotz(theta1)*RobotParameters::U1F1 +
                                    RobotParameters::K2*rotz(theta2)*RobotParameters::U2F2 +
                                    RobotParameters::K3*rotz(theta3)*RobotParameters::U3F3);

    Eigen::Vector3d ufW2 = K23inv*(RobotParameters::K2*rotz(theta2)*RobotParameters::U2F2 +
                                   RobotParameters::K3*rotz(theta3)*RobotParameters::U3F3);

    Eigen::Vector3d ufW3 = rotz(theta3)*RobotParameters::U3F3;

    // 3. create the twist xi = [u, v], v = [0,0,1]^T
    Eigen::Matrix<double,6,1> xi1, xi2, xi3;

    xi1 << ufW1,0,0,1;
    xi2 << ufW2,0,0,1;
    xi3 << ufW3,0,0,1;

    // 4. product of matrix expoentials
    Eigen::Matrix4d g = Lie::exp(xi1*s1) * Lie::exp(xi2*s2) * Lie::exp(xi3*s3);

    return g;
}

} // namespace CTR