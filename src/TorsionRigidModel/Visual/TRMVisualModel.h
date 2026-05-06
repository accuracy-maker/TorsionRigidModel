/**
 * @file TRMVisualModel.h
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Draw the CTR backbone as a coloured polyline using SOFA DrawTool
 * @version 0.1
 * @date 2026-05-06
 *
 * @copyright Copyright (c) 2026
 */

#pragma once

#include <sofa/core/visual/VisualModel.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/type/Vec.h>
#include <sofa/type/RGBAColor.h>

#include "core/ForwardKinematics.h"

namespace TRMCTR::visual
{

using sofa::core::objectmodel::Data;

/**
 * Draws the 3-section CTR backbone every render frame.
 *
 * Wire d_jointConfig from a TRMForwardKinematicsEngine (or IK controller):
 * @code
 * <TRMVisualModel name="visual"
 *     d_jointConfig="@fk.d_jointConfig"
 *     d_nSamples="20" />
 * @endcode
 *
 * Each section is drawn in a distinct colour:
 *   Section 1 (all tubes)   — blue
 *   Section 2 (tubes 2+3)   — green
 *   Section 3 (tube 3 only) — red
 */
class TRMVisualModel : public sofa::core::visual::VisualModel
{
public:
    SOFA_CLASS(TRMVisualModel, sofa::core::visual::VisualModel);

    using Vec6 = sofa::type::Vec<6, double>;

    // ------------------------------------------------------------------ Data
    Data<Vec6> d_jointConfig; ///< input: [θ1,s1,θ2,s2,θ3,s3] wired from FK engine
    Data<int>  d_nSamples;   ///< samples per section (default 20)

    // -------------------------------------------------------- SOFA life-cycle
    void doInitVisual(const sofa::core::visual::VisualParams* vparams) override;
    void doDrawVisual(const sofa::core::visual::VisualParams* vparams) override;

protected:
    TRMVisualModel();
    ~TRMVisualModel() override = default;
};

} // namespace TRMCTR::visual
