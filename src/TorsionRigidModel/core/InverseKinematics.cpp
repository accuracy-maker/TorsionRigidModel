/**
 * @file InverseKinematics.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Inverse Kinematics based on Jacobian Matrix
 * @version 0.1
 * @date 2026-04-30
 * 
 * @copyright Copyright (c) 2026
 * 
 */
#include "InverseKinematics.h"


#include <cmath>
#include <iostream>

namespace CTR {

Eigen::Matrix<double,6,6> InverseKinematics::Jacobian(const Eigen::Matrix<double,6,1>& q) {
    
    // extract the variables
    double theta1 = q(0), s1 = q(1), theta2 = q(2), s2 = q(3), theta3 = q(4), s3 = q(5);

    // compute the resultant curvature
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

    // create the twist xi = [u, v], where v = [0,0,1]^T
    Eigen::Matrix<double,6,1> xi1, xi2, xi3;

    xi1 << ufW1,0,0,1;
    xi2 << ufW2,0,0,1;
    xi3 << ufW3,0,0,1;

    // build the g1 and g2 for section 1 and 2
    Eigen::Matrix4d g1, g2;

    g1 = Lie::exp(xi1 * s1);
    g2 = Lie::exp(xi2 * s2);

    // build "A" matrix that is used to compute spaital velocity vee(g^{dot}g^{-1})
    Eigen::Matrix<double,6,6> A1, A2, A3;
    A1 = Lie::aMatrix(xi1, s1);
    A2 = Lie::aMatrix(xi2, s2);
    A3 = Lie::aMatrix(xi3, s3);

    // compute the derivative of xi w.r.t theta
    Eigen::Matrix<double,6,1> b11, b12, b13, b22, b23, b33; // b_{ij} ith section w.r.t jth tube

    Eigen::Matrix3d Jtheta1, Jtheta2, Jtheta3;

    Jtheta1 << -std::sin(theta1), -std::cos(theta1),0,
			std::cos(theta1), -std::sin(theta1), 0,
			0, 0, 0;
	Jtheta2 << -std::sin(theta2), -std::cos(theta2),0,
			std::cos(theta2), -std::sin(theta2), 0,
			0, 0, 0;
	Jtheta3 << -std::sin(theta3), -std::cos(theta3),0,
			std::cos(theta3), -std::sin(theta3), 0,
			0, 0, 0;

    b11 << K123inv*K1*Jtheta1*u1F1_,0,0,0;
	b12 << K123inv*K2*Jtheta2*u2F2_,0,0,0;
	b13 << K123inv*K3*Jtheta3*u3F3_,0,0,0;

	b22 << K23inv*K2*Jtheta2*u2F2_,0,0,0;
	b23 << K23inv*K3*Jtheta3*u3F3_,0,0,0;

	b33 << Jtheta3*u3F3_,0,0,0;

    // build adjoint transform
    Eigen::Matrix<double,6,6> Adg1 = Lie::Adg(g1), Adg1g2 = Lie::Adg(g1*g2);

    // compute the final Jacobian Matrix J.
    Eigen::Matrix<double,6,6> J;

    J << A1*b11, xi1,A1*b12+Adg1*A2*b22, Adg1*xi2, A1*b13+Adg1*A2*b23+Adg1g2*A3*b33, Adg1g2*xi3;

    return J;
}

Eigen::Matrix<double,6,1> InverseKinematics::IK(const Eigen::Vector3d& P, const Eigen::Vector3d& q0) {

}

}