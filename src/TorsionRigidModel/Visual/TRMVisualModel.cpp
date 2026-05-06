/**
 * @file TRMVisualModel.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Compute the visual model based on the joint configuration at each time step
 * @version 0.1
 * @date 2026-05-06
 *
 * @copyright Copyright (c) 2026
 */

#include "TRMVisualModel.h"
#include <sofa/core/ObjectFactory.h>

namespace TRMCTR::visual
{

// ----------------------------------------------------------------------------
// Constructor — register Data fields
// ----------------------------------------------------------------------------
TRMVisualModel::TRMVisualModel()
    : d_jointConfig(initData(&d_jointConfig,
                             Vec6(0, 50, 0, 50, 0, 50),
                             "d_jointConfig",
                             "[theta1,s1,theta2,s2,theta3,s3] (rad, mm)"))
    , d_nSamples(initData(&d_nSamples, 20, "d_nSamples",
                          "Number of sample points per section"))
{}

// ----------------------------------------------------------------------------
// doInitVisual — one-time setup (textures, shaders, etc.)
// Nothing needed for procedural line drawing.
// ----------------------------------------------------------------------------
void TRMVisualModel::doInitVisual(const sofa::core::visual::VisualParams* /*vparams*/)
{
    // TODO (optional): load a texture or set up a shader if you want
    //                  a fancier look than plain coloured lines.

}

// ----------------------------------------------------------------------------
// doDrawVisual — called every render frame
// ----------------------------------------------------------------------------
void TRMVisualModel::doDrawVisual(const sofa::core::visual::VisualParams* vparams)
{
    // -------------------------------------------------------------------------
    // 1. Read inputs
    // -------------------------------------------------------------------------
    const Vec6& q_sofa  = d_jointConfig.getValue();
    const int   N       = d_nSamples.getValue();

    Eigen::Matrix<double,6,1> q;
    for (int i = 0; i < 6; ++i) q(i) = q_sofa[i];

    // -------------------------------------------------------------------------
    // 2. Sample backbone → 3*N+1 points: [0..N] sec1, [N..2N] sec2, [2N..3N] sec3
    // -------------------------------------------------------------------------
    std::vector<Eigen::Vector3d> points = CTR::ForwardKinematics::backbone(q, N);

    // -------------------------------------------------------------------------
    // 3. Slice into per-section vectors and convert Eigen → sofa::type::Vec3d
    // -------------------------------------------------------------------------
    std::vector<Eigen::Vector3d> p1(points.begin(),       points.begin() + N + 1);
    std::vector<Eigen::Vector3d> p2(points.begin() + N,   points.begin() + 2*N + 1);
    std::vector<Eigen::Vector3d> p3(points.begin() + 2*N, points.end());

    std::vector<sofa::type::Vec3d> sec1, sec2, sec3;

    for (int i = 0; i <= N; ++i)
    {
        sec1.push_back({p1[i].x(), p1[i].y(), p1[i].z()});
        sec2.push_back({p2[i].x(), p2[i].y(), p2[i].z()});
        sec3.push_back({p3[i].x(), p3[i].y(), p3[i].z()});
    }

    // -------------------------------------------------------------------------
    // 4. Draw each section as cylinders and a tip sphere
    // -------------------------------------------------------------------------
    auto* dt = vparams->drawTool();

    auto drawTube = [&](const std::vector<sofa::type::Vec3d>& pts,
                        float radius,
                        const sofa::type::RGBAColor& color)
    {
        for (std::size_t i = 0; i + 1 < pts.size(); ++i)
            dt->drawCylinder(pts[i], pts[i + 1], radius, color);
    };

    drawTube(sec1, 1.5f, sofa::type::RGBAColor::blue());
    drawTube(sec2, 1.0f, sofa::type::RGBAColor::green());
    drawTube(sec3, 0.75f, sofa::type::RGBAColor::red());

    sofa::type::Vec3d tip = sec3.back();
    dt->drawSpheres({tip}, 1.5f, sofa::type::RGBAColor::white());
}


// ----------------------------------------------------------------------------
// computeBBox — gives qglviewer a correct scene radius so near/far clip planes
// are computed properly; without this the viewer clips geometry at some angles.
// ----------------------------------------------------------------------------
void TRMVisualModel::computeBBox(const sofa::core::ExecParams* /*params*/, bool /*onlyVisible*/)
{
    const Vec6& q_sofa = d_jointConfig.getValue();
    const int   N      = d_nSamples.getValue();

    Eigen::Matrix<double,6,1> q;
    for (int i = 0; i < 6; ++i) q(i) = q_sofa[i];

    std::vector<Eigen::Vector3d> points = CTR::ForwardKinematics::backbone(q, N);

    constexpr double margin = 2.0; // accounts for the widest tube radius (1.5 mm)
    sofa::type::Vec3d minBB( 1e10,  1e10,  1e10);
    sofa::type::Vec3d maxBB(-1e10, -1e10, -1e10);

    for (const auto& p : points)
    {
        minBB[0] = std::min(minBB[0], p.x() - margin);
        minBB[1] = std::min(minBB[1], p.y() - margin);
        minBB[2] = std::min(minBB[2], p.z() - margin);
        maxBB[0] = std::max(maxBB[0], p.x() + margin);
        maxBB[1] = std::max(maxBB[1], p.y() + margin);
        maxBB[2] = std::max(maxBB[2], p.z() + margin);
    }

    this->f_bbox.setValue(sofa::type::BoundingBox(minBB, maxBB));
}

void registerTRMVisualModel(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("CTR backbone visual model: draws the 3-section tube backbone as a coloured polyline.")
        .add<TRMVisualModel>());
}

} // namespace TRMCTR::visual
