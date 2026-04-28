#pragma once

#include <Eigen/Dense>

/* Implementation of Lie Group SE(3) and associated Lie algebra se(3) */
namespace Lie {
    /**
     * @brief The "hat" operator for so(3)
     * Maps a 3D vector to a 3x3 skew-symmetric matrix: R^3 -> R^{3x3}
    */
    inline Eigen::Matrix3d hat(const Eigen::Vector3d& omega) {
        Eigen::Matrix3d Omega_hat;

        Omega_hat << 0, -omega(2), omega(1),
                    omega(2), 0, -omega(0),
                    -omega(1), omega(0), 0;
        
        return Omega_hat;
    }

    /**
     * @brief The "vee" operator for so(3)
     * Maps a 3x3 Matrix back to its 3D vector: R^{3x3} -> R^3
    */
    inline Eigen::Vector3d vee(const Eigen::Matrix3d& Omega_hat) {
        Eigen::Vector3d omega;

        omega << Omega_hat(2, 1), Omega_hat(0, 2), Omega_hat(1, 0);

        return omega;
    }

    /**
     * @brief The "hat" operator for se(3)
     * Maps a 6D twist \xi=[omega, v] to 4x4 matrix: R^6 -> R^{4x4}
    */
    inline Eigen::Matrix4d hat(const Eigen::Matrix<double, 6, 1>& xi) {
        Eigen::Vector3d omega = xi.head<3>();
        Eigen::Vector3d v = xi.tail<3>();

        Eigen::Matrix4d Xi_hat = Eigen::Matrix4d::Zero();
        Xi_hat.block<3, 3>(0, 0) = hat(omega);
        Xi_hat.block<3, 1>(0, 3) = v;
        return Xi_hat;
    }

    /**
     * @brief The "vee" operator for se(3)
     * Maps a 4x4 Matrix back to its 6D twist \xi=[omega, v]: R^{4x4} -> R^6
    */
    inline Eigen::Matrix<double, 6, 1> vee(const Eigen::Matrix4d& Xi_hat) {
        Eigen::Matrix<double, 6, 1> xi;
        xi.head<3>() = vee(Eigen::Matrix3d(Xi_hat.block<3, 3>(0, 0)));
        xi.tail<3>() = Xi_hat.block<3, 1>(0, 3);
        return xi;
    }

    /**
     * @brief Exponentail Map: se(3) -> SE(3)
    */
    inline Eigen::Matrix4d exp(const Eigen::Matrix<double, 6, 1>& xi) {
        Eigen::Vector3d omega = xi.head<3>();
        Eigen::Vector3d v = xi.tail<3>();
        double theta = omega.norm();

        Eigen::Matrix4d g = Eigen::Matrix4d::Identity();

        // pure translation
        if (theta < 1e-8) {
            g.block<3,1>(0,3) = v;
        } else {
            // rotation + translation
            Eigen::Vector3d axis = omega / theta;
            Eigen::Matrix3d Omega_hat = hat(axis);
            // R = e^{Omega_hat*\theta} = I + Omega_hat*sin(theta) + Omega_hat^2 * (1 - cos(theta))
            Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + 
                                Omega_hat * std::sin(theta) + 
                                Omega_hat * Omega_hat * (1 - std::cos(theta));

            // V = I*theta + (1-cos(theta))/theta * Omega_hat * theta + (theta - sin(theta)) / (theta^2) * (Omega_hat*theta)^2
            Eigen::Matrix3d V = Eigen::Matrix3d::Identity() * theta + 
                                (1 - std::cos(theta)) * Omega_hat + 
                                (theta - std::sin(theta)) * (Omega_hat * Omega_hat);

            // p
            Eigen::Vector3d p = V * (v / theta);

            g.block<3,3>(0,0) = R;
            g.block<3,1>(0,3) = p;
        }
        return g;
    }

} // namespace Lie
