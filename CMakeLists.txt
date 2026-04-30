cmake_minimum_required(VERSION 3.14)
project(TorsionRigidModel)

# Use C++17 for Eigen and GTest compatibility
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Find Dependencies
find_package(Eigen3 REQUIRED)
# This looks for GTest on your system (e.g., installed via apt, brew, or vcpkg)
find_package(GTest REQUIRED)

# 2. Include Directories
# This allows you to #include "math/Lie.h" from your test files
include_directories(src/TorsionRigidModel)

# 3. Create the Test Executable
add_executable(test_lie src/TorsionRigidModel/test/test_lie.cpp)

# 4. Link Libraries
# We link Eigen3 and the Google Test main entry point
target_link_libraries(test_lie 
    Eigen3::Eigen 
    GTest::gtest_main
)

add_executable(test_fk
    src/TorsionRigidModel/test/test_fk.cpp
    src/TorsionRigidModel/core/ForwardKinematics.cpp
)

target_link_libraries(test_fk
    Eigen3::Eigen
    GTest::gtest_main
)

# 5. Enable CTest
enable_testing()
add_test(NAME LieMathTest COMMAND test_lie)
add_test(NAME ForwardKinematicsTest COMMAND test_fk)