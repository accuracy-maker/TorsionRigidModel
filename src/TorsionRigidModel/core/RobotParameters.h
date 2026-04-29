#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace CTR {
    /**
     * @brief Physical and Operational parameters for the 3-tube CTR
    */
    struct RobotParameters {
        // stiffness matrix
        inline static const Eigen::Matrix3d K1 = Eigen::Matrix3d::Identity() * 100.0;
        inline static const Eigen::Matrix3d K2 = Eigen::Matrix3d::Identity() * 10.0;
        inline static const Eigen::Matrix3d K3 = Eigen::Matrix3d::Identity() * 0.0;

        // curvature vector
        inline static const Eigen::Vector3d U1F1 = Eigen::Vector3d(0, 0, 0);
        inline static const Eigen::Vector3d U2F2 = Eigen::Vector3d(0, 1.0 / 173.58, 0);
        inline static const Eigen::Vector3d U3F3 = Eigen::Vector3d(0, 0, 0);

        // Joint Limit
        static constexpr double S1_MAX = 100.0;
        static constexpr double S1_MIN = 10.0;

        static constexpr double S2_MAX = 100.0;
        static constexpr double S2_MIN = 10.0;

        static constexpr double S3_MAX = 100.0;
        static constexpr double S3_MIN = 25.0;

        // Rotation Limit
        static constexpr double THETA_MAX = M_PI;
        static constexpr double THETA_MIN = -M_PI;
    };
} // namespace CTR