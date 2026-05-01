/**
 * @file test_fk_sofa_engine.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Test the forward kinematics based on sofa engine framework
 * @version 0.1
 * @date 2026-05-01
 *
 * @copyright Copyright (c) 2026
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>

#include "Engine/TRMForwardKinematicsEngine.h"

/**
 * @brief Test for Sofa Engine FK
 */
class SofaEngineTest : public ::testing::Test {
protected:
    using EngineType  = TRMCTR::engine::TRMForwardKinematicsEngine;
    using ExpectedIn  = sofa::core::objectmodel::Data<sofa::type::Vec<6, double>>;
    using ExpectedOut = sofa::core::objectmodel::Data<sofa::defaulttype::Rigid3dTypes::Coord>;

    void SetUp() override {}
};
