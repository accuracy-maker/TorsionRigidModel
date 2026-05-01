/**
 * @file TRMForwardKinematicsEngine.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Data Engine of Sofa Framework based on the implementation of FK of Torsion Rigid Model (TRM)
 * @version 0.1
 * @date 2026-05-01
 *
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <sofa/core/DataEngine.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>

#include "core/ForwardKinematics.h"

namespace TRMCTR::engine
{

/**
 * SOFA DataEngine wrapping CTR::ForwardKinematics::FK().
 *
 * Reads the 6-DOF joint configuration every time it changes and writes
 * the SE(3) end-effector pose (position + quaternion) as a Data output.
 *
 * Scene-graph usage:
 * @code
 * <TRMForwardKinematicsEngine
 *     name="fkEngine"
 *     d_jointConfig="0 50 0 50 0 50" />
 * @endcode
 *
 * Inputs
 *   d_jointConfig     [θ1, s1, θ2, s2, θ3, s3]   (rad, mm)
 *
 * Outputs
 *   d_endEffectorPose  SE(3) pose as (position, quaternion)
 */
class TRMForwardKinematicsEngine : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(TRMForwardKinematicsEngine, sofa::core::DataEngine);

    using Rigid3Coord = sofa::defaulttype::Rigid3dTypes::Coord;
    using Vec6        = sofa::type::Vec<6, double>;

    // ------------------------------------------------------------------ Data
    Data<Vec6>        d_jointConfig;     ///< input:  [θ1, s1, θ2, s2, θ3, s3]
    Data<Rigid3Coord> d_endEffectorPose; ///< output: position + quaternion

    // -------------------------------------------------------- SOFA life-cycle
    void init()     override;
    void reinit()   override;
    void doUpdate() override;

protected:
    TRMForwardKinematicsEngine();
    ~TRMForwardKinematicsEngine() override = default;

private:
    /// Convert the 4×4 SE(3) Eigen matrix returned by FK() to a SOFA Rigid3Coord.
    static Rigid3Coord toRigid3(const Eigen::Matrix4d& T);
};

} // namespace TRMCTR::engine
