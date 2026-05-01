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

    Eigen::Vector3d ufW1 = K123inv*(RobotParameters::K1*ForwardKinematics::rotz(theta1)*RobotParameters::U1F1 +
                                    RobotParameters::K2*ForwardKinematics::rotz(theta2)*RobotParameters::U2F2 +
                                    RobotParameters::K3*ForwardKinematics::rotz(theta3)*RobotParameters::U3F3);

    Eigen::Vector3d ufW2 = K23inv*(RobotParameters::K2*ForwardKinematics::rotz(theta2)*RobotParameters::U2F2 +
                                   RobotParameters::K3*ForwardKinematics::rotz(theta3)*RobotParameters::U3F3);

    Eigen::Vector3d ufW3 = ForwardKinematics::rotz(theta3)*RobotParameters::U3F3;

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

    Eigen::Matrix3d Jtheta1, Jtheta2, Jtheta3; // derivative of rotz w.r.t theta
    Jtheta1 << -std::sin(theta1), -std::cos(theta1),0,
			std::cos(theta1), -std::sin(theta1), 0,
			0, 0, 0;
	Jtheta2 << -std::sin(theta2), -std::cos(theta2),0,
			std::cos(theta2), -std::sin(theta2), 0,
			0, 0, 0;
	Jtheta3 << -std::sin(theta3), -std::cos(theta3),0,
			std::cos(theta3), -std::sin(theta3), 0,
			0, 0, 0;

    b11 << K123inv*RobotParameters::K1*Jtheta1*RobotParameters::U1F1,0,0,0;
	b12 << K123inv*RobotParameters::K2*Jtheta2*RobotParameters::U2F2,0,0,0;
	b13 << K123inv*RobotParameters::K3*Jtheta3*RobotParameters::U3F3,0,0,0;

	b22 << K23inv*RobotParameters::K2*Jtheta2*RobotParameters::U2F2,0,0,0;
	b23 << K23inv*RobotParameters::K3*Jtheta3*RobotParameters::U3F3,0,0,0;

	b33 << Jtheta3*RobotParameters::U3F3,0,0,0;

    // build adjoint transform
    Eigen::Matrix<double,6,6> Adg1 = Lie::Adg(g1), Adg1g2 = Lie::Adg(g1*g2);

    // compute the final Jacobian Matrix J.
    Eigen::Matrix<double,6,6> J;

    J << A1*b11, xi1,A1*b12+Adg1*A2*b22, Adg1*xi2, A1*b13+Adg1*A2*b23+Adg1g2*A3*b33, Adg1g2*xi3;

    return J;
}

Eigen::Matrix<double,6,1> InverseKinematics::IK(const Eigen::Vector3d& P, const Eigen::Matrix<double,6,1>& q0) {
    /**
     * @brief one-iteration updating
     * q_{k+1} = q_k + \delta-q
     */

    // translation error of the P_{desired} - q_{current}
    Eigen::Vector3d e_Z = Eigen::Vector3d::Ones(3,1); 
    
    // Translation part of the matrix g
    Eigen::Vector3d Pe; 

    // Initialise q
    Eigen::Matrix<double,6,1> q = q0;

    // gn \in SE(3)
    Eigen::Matrix4d gn;

    // q_{new} = q + delta_q
    Eigen::Matrix<double,6,1> delta_q;

    Eigen::Matrix<double,6,6> J0, W;

    Eigen::Matrix<double,3,6> J_Z, J_ZPe;
    Eigen::Matrix3d H;

    // Forward a small step
    gn = ForwardKinematics::FK(q);
    
    // translation error
    e_Z = P - gn.block<3,1>(0,3);

    J0 = Jacobian(q);

    Pe = gn.block<3,1>(0,3);

    // Jacobian of Pe
    J_ZPe << -Lie::hat(Pe), Eigen::Matrix3d::Identity();

    J_Z = J_ZPe * J0;

    W = Eigen::Matrix<double,6,6>::Identity();

    W(0,0) = 1;
    W(1,1) = 10000;
    W(2,2) = 1;
    W(3,3) = 10000;
    W(4,4) = 1;
    W(5,5) = 10000;

    H = J_Z*W*J_Z.transpose()+(0.5*e_Z.norm()*e_Z.norm())*Eigen::Matrix3d::Identity();

    delta_q = W*J_Z.transpose()*H.inverse()*e_Z;
    
    q = q + delta_q;

    // length constraint
    if (q(1) > RobotParameters::S1_MAX) {
        q(1) = RobotParameters::S1_MAX;
    } else if (q(1) < RobotParameters::S1_MIN) {
        q(1) = RobotParameters::S1_MIN;
    }

    if (q(3) > RobotParameters::S2_MAX) {
        q(3) = RobotParameters::S2_MAX;
    } else if (q(3) < RobotParameters::S2_MIN) {
        q(3) = RobotParameters::S2_MIN;
    }

    if (q(5) > RobotParameters::S3_MAX) {
        q(5) = RobotParameters::S3_MAX;
    } else if (q(5) < RobotParameters::S3_MIN) {
        q(5) = RobotParameters::S3_MIN;
    }


    return q;
}

}