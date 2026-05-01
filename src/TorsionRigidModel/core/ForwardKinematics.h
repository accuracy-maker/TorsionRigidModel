/**
 * @file ForwardKinematics.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Forward Kinematics head file.
 * @version 0.1
 * @date 2026-04-30
 * 
 * @copyright Copyright (c) 2026
 * 
 */
#pragma once

#include "math/Lie.h"
#include "RobotParameters.h"

#include <Eigen/Dense>

namespace CTR {
    class ForwardKinematics {
    public:
        /**
         * @brief FK is a function that maps the joint configuration space Q to end-effector configuration SE(3)
         * FK: g_st: Q -> SE(3)
        */
        static Eigen::Matrix4d FK(const Eigen::Matrix<double,6,1>& q);

        /// Rotation matrix about the z-axis by angle theta (radians)
        static Eigen::Matrix3d rotz(double theta);
    };
} // namespace CTR