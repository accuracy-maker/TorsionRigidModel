/* Test the math/lie.h file by google testing */
#include <gtest/gtest.h>
#include "math/Lie.h"
#include <Eigen/Dense>

const double EPS = 1e-8;

/**
 * @brief Text fixture for Lie Group Operation
*/
class LieTest : public ::testing::Test {
protected:
    void SetUp() override {

    }
};

// 1. check the hat operator for so(3)
TEST_F(LieTest, TestSo3HatAndPrint) {
    // 1. Define a 3D vector as input
    Eigen::Vector3d omega(1.0, 2.0, 3.0);

    // 2. Call the so(3) hat operator
    Eigen::Matrix3d Omega_hat = Lie::hat(omega);

    // 3. Check if the output dimensions are 3x3
    // For fixed-size Eigen matrices, these are always 3
    EXPECT_EQ(Omega_hat.rows(), 3);
    EXPECT_EQ(Omega_hat.cols(), 3);

    // 4. Verify skew-symmetry: Omega_hat should be -Omega_hat^T
    Eigen::Matrix3d sum = Omega_hat + Omega_hat.transpose();
    EXPECT_TRUE(sum.isMuchSmallerThan(Eigen::Matrix3d::Identity()));

    // 5. Print the matrix to the console
    std::cout << "\nInput Vector (omega):\n" << omega << std::endl;
    std::cout << "Output so(3) Matrix (Omega_hat):\n" << Omega_hat << std::endl;
}

// 2. check the vee operator for so(3)
TEST_F(LieTest, TestSo3VeeAndPrint) {
    Eigen::Matrix3d Omega_hat;

    Omega_hat << 0, -3, 2,
                3, 0, -1,
                -2, 1, 0;

    Eigen::Vector3d omega = Lie::vee(Omega_hat);

    std::cout << "Input so(3) Matrix:\n" << Omega_hat << std::endl;
    std::cout << "Recovered Vector (vee):\n" << omega << std::endl;

    EXPECT_NEAR(omega(0), 1.0, EPS);
    EXPECT_NEAR(omega(1), 2.0, EPS);
    EXPECT_NEAR(omega(2), 3.0, EPS);
}

// 3. Test the Hat and Vee operators for so(3)
TEST_F(LieTest, So3HatVeeConsistency) {
    Eigen::Vector3d omega(0.1, -0.2, 0.5);
    Eigen::Matrix3d Omega_hat = Lie::hat(omega);
    Eigen::Vector3d omega_recovered = Lie::vee(Omega_hat);

    // check vee(hat(omega)) == omega
    EXPECT_NEAR((omega - omega_recovered).norm(), 0.0, EPS);

    // check skew-symmetric: Omega_hat^T = -Omega_hat
    EXPECT_NEAR((Omega_hat.transpose() + Omega_hat).norm(), 0.0, EPS);
}

// 4. Test the Hat operator for se(3)
TEST_F(LieTest, TestSe3HatAndPrint) {
    Eigen::Vector3d omega(1.0, 2.0, 3.0);
    Eigen::Vector3d v(1.0, 2.0, 3.0);

    Eigen::Matrix<double, 6, 1> xi;
    xi.head<3>() = omega;
    xi.tail<3>() = v;

    Eigen::Matrix4d Xi_hat = Lie::hat(xi);

    EXPECT_EQ(Xi_hat.rows(), 4);
    EXPECT_EQ(Xi_hat.cols(), 4);

    std::cout << "Input 6d twist se(3):\n" << xi << std::endl;
    std::cout << "Output 4x4 matrix SE(3):\n" << Xi_hat << std::endl;
}

// 5. Test the vee operator for SE(3)
TEST_F(LieTest, TestSE3VeeAndPrint) {
    Eigen::Matrix4d Xi_hat;
    Xi_hat << 0, -3, 2, 1,
              3, 0, -1, 2,
              -2, 1, 0, 3,
              0, 0, 0, 0;

    Eigen::Matrix<double, 6, 1> xi = Lie::vee(Xi_hat);

    EXPECT_EQ(xi.rows(), 6);
    EXPECT_EQ(xi.cols(), 1);

    Eigen::Vector3d omega = xi.head<3>();
    Eigen::Vector3d v = xi.tail<3>();

    EXPECT_NEAR(omega(0), 1.0, EPS);
    EXPECT_NEAR(omega(1), 2.0, EPS);
    EXPECT_NEAR(omega(2), 3.0, EPS);
    EXPECT_NEAR(v(0), 1.0, EPS);
    EXPECT_NEAR(v(1), 2.0, EPS);
    EXPECT_NEAR(v(2), 3.0, EPS);

    std::cout << "Input 4x4 SE(3) Matrix:\n" << Xi_hat << std::endl;
    std::cout << "Output 6D twist se(3): \n" << xi << std::endl;
}

// 6. Test the Matrix Exponential
TEST_F(LieTest, MatrixExponential) {
    Eigen::Vector3d omega(0, 2*M_PI, 0);
    Eigen::Vector3d v(0, 0, 1);
    Eigen::Matrix<double, 6, 1> xi;
    xi.head<3>() = omega;
    xi.tail<3>() = v;

    Eigen::Matrix4d g = Lie::exp(xi);

    Eigen::Matrix3d R = g.block<3,3>(0,0);
    Eigen::Vector3d p = g.block<3,1>(0,3);

    EXPECT_NEAR(R.norm(), Eigen::Matrix3d::Identity().norm(), EPS);
    EXPECT_NEAR(p(0), 0.0, EPS);
    EXPECT_NEAR(p(1), 0.0, EPS);
    EXPECT_NEAR(p(2), 0.0, EPS);

    std::cout << "Input 6d se(3) twist:\n" << xi << std::endl;
    std::cout << "Output 4x4 SE(3) Matrix:\n" << g << std::endl;
    std::cout << "Rotation Component:\n" << R << std::endl;
    std::cout << "translation Component:\n" << p << std::endl;
}