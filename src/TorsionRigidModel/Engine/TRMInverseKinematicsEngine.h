/**
 * @file TRMInverseKinematicsEngine.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Data Engine of Sofa Framework based on the implementation of IK of Torsion Rigid model.
 * @version 0.1
 * @date 2026-05-01
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once

#include <sofa/core/DataEngine.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>

#include "core/InverseKinematics.h"

namespace TRMCTR::engine
{

/**
 * SOFA DataEngine wrapping CTR::InverseKinematics::IK().
 * 
 * Reads the target pose (position + quaternion) as a Data In.
 * Output the 6 DoF joint configuration as Data Output.
 */

class TRMInverseKinematicsEngine : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(TRMInverseKinematicsEngine, sofa::core::DataEngine);

    using Rigid3Coord = sofa::defaulttype::Rigid3dTypes::Coord;

    using Vec6        = sofa::type::Vec<6, double>;

    Data<Rigid3Coord> target_pos;
    Data<Vec6> d_jointConfig;

    // SOFA live-cycle
    void init() override;
    void reinit() override;
    void doUpdate() override;

protected:
    TRMInverseKinematicsEngine();
    ~TRMInverseKinematicsEngine() override = default;

private:
    
}


}