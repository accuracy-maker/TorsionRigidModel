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

    private:
        /// Rotation matrix about the z-axis by angle theta (radians)
        static Eigen::Matrix3d rotz(double theta);
    };
} // namespace CTR