cmake_minimum_required(VERSION 3.22)

project(TorsionRigidModel VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Sofa.Framework REQUIRED)
find_package(Sofa.Simulation.Core REQUIRED)
find_package(Sofa.GL REQUIRED)
find_package(Eigen3 REQUIRED)

set(TRM_SOURCE_DIR src/TorsionRigidModel)

set(HEADER_FILES
    ${TRM_SOURCE_DIR}/math/Lie.h
    ${TRM_SOURCE_DIR}/core/RobotParameters.h
    ${TRM_SOURCE_DIR}/core/ForwardKinematics.h
    ${TRM_SOURCE_DIR}/core/InverseKinematics.h
    ${TRM_SOURCE_DIR}/Engine/TRMForwardKinematicsEngine.h
    ${TRM_SOURCE_DIR}/Controller/TRMInverseKinematicsEngine.h
    ${TRM_SOURCE_DIR}/Visual/TRMVisualModel.h
)

set(SOURCE_FILES
    ${TRM_SOURCE_DIR}/initTorsionRigidModel.cpp
    ${TRM_SOURCE_DIR}/core/ForwardKinematics.cpp
    ${TRM_SOURCE_DIR}/core/InverseKinematics.cpp
    ${TRM_SOURCE_DIR}/Engine/TRMForwardKinematicsEngine.cpp
    ${TRM_SOURCE_DIR}/Controller/TRMInverseKinematicsEngine.cpp
    ${TRM_SOURCE_DIR}/Visual/TRMVisualModel.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME}
    PUBLIC  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
            "$<INSTALL_INTERFACE:include>"
    PRIVATE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/TorsionRigidModel>"
)

target_link_libraries(${PROJECT_NAME}
    Sofa.Core
    Sofa.Simulation.Core
    Sofa.GL
    Eigen3::Eigen
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "-DSOFA_BUILD_TORSIONRIGIDMODEL"
)

cmake_dependent_option(TRM_BUILD_TESTS "Build TorsionRigidModel tests" ON
    "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)
if(TRM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(src/TorsionRigidModel/test)
endif()
