/**
 * @file test_lie.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Testing for Lie group and lie algebra
 * @version 0.1
 * @date 2026-04-30
 * 
 * @copyright Copyright (c) 2026
 * 
 */
#include <gtest/gtest.h>
#include "math/Lie.h"
#include <Eigen/Dense>

const double EPS = 1e-8;

/**
 * @brief Test fixture for Lie Group Operation
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

// 7. Test adjoint transform
TEST_F(LieTest, AdjointTransform) {
    Eigen::Vector3d omega(0, 2*M_PI, 0);
    Eigen::Vector3d v(0, 0, 1);
    Eigen::Matrix<double, 6, 1> xi;
    xi.head<3>() = omega;
    xi.tail<3>() = v;

    Eigen::Matrix4d g = Lie::exp(xi);

    Eigen::Matrix<double,6,6> ad = Lie::Adg(g);

    // check if it is a 6 x 6 matrix
    EXPECT_EQ(ad.rows(), 6);
    EXPECT_EQ(ad.cols(), 6);

    // check if different components are correct
    Eigen::Matrix3d R = ad.block<3,3>(0,0);
    Eigen::Matrix3d z = ad.block<3,3>(0,3);
    Eigen::Matrix3d p_hat_R = ad.block<3,3>(3,0);
    Eigen::Matrix3d R_ = ad.block<3,3>(3,3);

    EXPECT_EQ(R, R_);
    EXPECT_NEAR(R.norm(), Eigen::Matrix3d::Identity().norm(), EPS);
    EXPECT_NEAR(R_.norm(), Eigen::Matrix3d::Identity().norm(), EPS);
    EXPECT_EQ(z, Eigen::Matrix3d::Zero());
    EXPECT_NEAR(p_hat_R.norm(), 0.0, EPS);

    // print the matrix
    std::cout << "Matrix Exponential g:\n" << g << std::endl;
    std::cout << "Adjoint Transform Matrix ad:\n" << ad << std::endl;

}

// 8. Test "A" Matrix with norm = 0
TEST_F(LieTest, AMatrixIdentity) {
    // test if the norm is zero
    Eigen::Vector3d w(0,0,0);
    Eigen::Vector3d v(0,0,1);

    Eigen::Matrix<double,6,1> xi;
    xi.head<3>() = w;
    xi.tail<3>() = v;

    double s = 2;

    Eigen::Matrix<double,6,6> aM = Lie::aMatrix(xi, s);

    // check if aM is a 6x6 matrix
    EXPECT_EQ(aM.rows(), 6);
    EXPECT_EQ(aM.cols(), 6);

    // check if it's an Identity matrix
    const double identityNorm = Eigen::Matrix<double,6,6>::Identity().norm();
    EXPECT_NEAR(aM.norm(), s*identityNorm, EPS);

    std::cout << "Input xi:\n" << xi << std::endl;
    std::cout << "Input s:\n" << s << std::endl;
    std::cout << "Output aM:\n" << aM << std::endl;
}

// 9. Test "A" Matrix in general
TEST_F(LieTest, AMatrixGeneral) {
    Eigen::Vector3d w(0,1,0);
    Eigen::Vector3d v(0,0,1);

    Eigen::Matrix<double,6,1> xi;
    xi.head<3>() = w;
    xi.tail<3>() = v;

    double s = 2*M_PI;

    Eigen::Matrix<double,6,6> aM = Lie::aMatrix(xi, s);

    // check if aM is a 6x6 matrix
    EXPECT_EQ(aM.rows(), 6);
    EXPECT_EQ(aM.cols(), 6);

    // check the result is equal to what we expect
    // (4 - 2π·0 - 4·1)/(2·1²) = 0
    // (4·2π - 5·0 + 2π·1)/(2·1³) = 10π/2 = 5π
    // (2 - 2π·0 - 2·1)/(2·1⁴) = 0
    // (2·2π - 3·0 + 2π·1)/(2·1⁵) = 6π/2 = 3π
    // aM = s·I + 5π·Ω² + 3π·Ω⁴
    Eigen::Matrix<double,6,6> Omega;
    Omega << Lie::hat(w), Eigen::Matrix3d::Zero(),
            Lie::hat(v), Lie::hat(w);
    
    const double expected = (Eigen::Matrix<double,6,6>::Identity() * s + 5*M_PI*Omega*Omega + 3*M_PI*Omega*Omega*Omega*Omega).norm();
    EXPECT_NEAR(aM.norm(), expected, EPS);

    std::cout << "Input xi:\n" << xi << std::endl;
    std::cout << "Input s:\n" << s << std::endl;
    std::cout << "Output aM:\n" << aM << std::endl;
}