/**
 * @file InverseKinematics.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Inverse Kinematics based on Jacobian Matrix head file
 * @version 0.1
 * @date 2026-04-30
 * 
 * @copyright Copyright (c) 2026
 * 
 */
#pragma once

#include "RobotParameters.h"
#include "math/Lie.h"

#include <Eigen/Dense>

namespace CTR {
    class InverseKinematics{
    public:
        /**
         * @brief IK if a function that maps end-effector configuration back to the joint configuration space
         * IK: SE(3) -> Q
         * @param: P <- target point
         * @param: q0 <- current position
         * @return: joint configuration \delta q
        */

        static Eigen::Matrix<double,6,1> IK(const Eigen::Vector3d& P, const Eigen::Vector3d& q0);
    
    private:
        /**
         * @brief Jacobian is a matrix that maps individual joint's velocity to the end-effector's velocity
         * Traditionally, it's just a differential representation of FK mapping
         * In screw motion, we use g^{dot}g^{-1} to represent the instantaneous spatial velocity so that we have a different formulas of Jacobian matrix.
        */
        static Eigen::Matrix<double,6,6> Jacobian(const Eigen::Matrix<double,6,1>& q);
    }
}