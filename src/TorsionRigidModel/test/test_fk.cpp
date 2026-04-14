#include "types/TubeParameters.h"
#include "types/Section.h"
#include "types/Frame.h"
#include "core/ForwardKinematics.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

using namespace ctr;

constexpr double NITINOL_E  = 58.0e9;
constexpr double NITINOL_NU = 0.33;
constexpr double mm = 1.0e-3;
constexpr double cm = 1.0e-2;

/* Test 1: Straight tube -> tip should be at (0,0,L) */
void test_straight_tube() {
    double L = 0.10;  // 10cm
    std::vector<TubeParameters> tubes = {
        makeTube(2.0*mm, 1.5*mm, {{0, 0, L}}, 1.0)
    };
 
    ForwardKinematics fk(tubes);
    double alpha[] = {0.0};
    double beta[]  = {0.0};
 
    Frame tip = fk.computeTipFrame(alpha, beta, 1);
 
    std::cout << "  Tip: " << tip << std::endl;
    std::cout << "  Expected: p=(0, 0, " << L << ")" << std::endl;
 
    assert(std::abs(tip.p.x) < 1e-12);
    assert(std::abs(tip.p.y) < 1e-12);
    assert(std::abs(tip.p.z - L) < 1e-12);
 
    // Orientation should be identity (no bending)
    assert(std::abs(tip.R.m[0][0] - 1.0) < 1e-12);
    assert(std::abs(tip.R.m[1][1] - 1.0) < 1e-12);
    assert(std::abs(tip.R.m[2][2] - 1.0) < 1e-12);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 2: Single curved tube -> circular arc */
void test_single_curved_tube() {
    std::cout << "=== Test: Single curved tube (circular arc) ===" << std::endl;

    double kappa_y = 20.0; // ky = 20.0, rotate in xz plane
    double radius = 1.0 / kappa_y; 
    double L = 0.10; // 10cm arc length
    double theta = kappa_y * L; // bending angle

    std::vector<TubeParameters> tubes = {
        makeTube(2.0*mm, 1.5*mm, {{0, kappa_y, L}}, 1.0)
    };

    ForwardKinematics fk(tubes);

    double alpha[] = {0.0};
    double beta[] = {0.0};

    Frame tip = fk.computeTipFrame(alpha, beta, 1);

    // Expected tip position: circular arc in xz plane
    double expected_x = radius * (1 - std::cos(theta));
    double expected_y = 0.0;
    double expected_z = radius * std::sin(theta);

    std::cout << "  κ = " << kappa_y << " (1/m), radius = " << radius*100 << " cm" << std::endl;
    std::cout << "  Arc length = " << L*100 << " cm, bend angle = "
              << theta * 180 / M_PI << " deg" << std::endl;
    std::cout << "  Tip: " << tip << std::endl;
    std::cout << "  Expected: p=(" << expected_x << ", " << expected_y
              << ", " << expected_z << ")" << std::endl;

    assert(std::abs(tip.p.x - expected_x) < 1e-10);
    assert(std::abs(tip.p.y - expected_y) < 1e-10);
    assert(std::abs(tip.p.z - expected_z) < 1e-10);

    // check tangent direction
    // R[0][2] = ky*sinθ/n = sinθ, R[2][2] = cosθ
    Vec3 tangent = tip.tangent();
    double expected_tan_x = std::sin(theta);
    double expected_tan_z = std::cos(theta);
    std::cout << "  Tangent: (" << tangent.x << ", " << tangent.y << ", " << tangent.z << ")" << std::endl;
 
    assert(std::abs(tangent.x - expected_tan_x) < 1e-10);
    assert(std::abs(tangent.y) < 1e-10);
    assert(std::abs(tangent.z - expected_tan_z) < 1e-10);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 3: 90 degree bending -> tip tangent should be perpendicular to z */
void test_90_degree_bend() {
    std::cout << "=== Test: 90-degree bend ===" << std::endl;

    double kappa = 20.0;
    double radius = 1.0/kappa;
    // Arc length for 90 degrees: s = (π/2) / κ
    double L = (M_PI / 2.0) / kappa;

    std::vector<TubeParameters> tubes = {
        makeTube(2.0*mm, 1.5*mm, {{0, kappa, L}}, 1.0)
    };

    ForwardKinematics fk(tubes);
    double alpha[] = {0.0};
    double beta[]  = {0.0};
 
    Frame tip = fk.computeTipFrame(alpha, beta, 1);
 
    // After 90° bend with κy: tip is at (radius, 0, radius)
    std::cout << "  Tip: " << tip << std::endl;
    std::cout << "  Expected: p=(" << radius << ", 0, " << radius << ")" << std::endl;
 
    assert(std::abs(tip.p.x - radius) < 1e-10);
    assert(std::abs(tip.p.y) < 1e-10);
    assert(std::abs(tip.p.z - radius) < 1e-10);
 
    // Tangent should point in +x direction (z-column of R after 90° bend about x)
    Vec3 tangent = tip.tangent();
    std::cout << "  Tangent: (" << tangent.x << ", " << tangent.y << ", " << tangent.z << ")" << std::endl;
 
    assert(std::abs(tangent.x - 1.0) < 1e-10);
    assert(std::abs(tangent.y) < 1e-10);
    assert(std::abs(tangent.z) < 1e-10);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 4: 180 degree bending -> semicircular */
void test_180_degree_bend() {
    std::cout << "=== Test: 180-degree bend (semicircle) ===" << std::endl;
 
    double kappa = 20.0;
    double radius = 1.0/kappa;
    double L = M_PI / kappa;
 
    std::vector<TubeParameters> tubes = {
        makeTube(2.0*mm, 1.5*mm, {{0, kappa, L}}, 1.0)
    };
 
    ForwardKinematics fk(tubes);
    double alpha[] = {0.0};
    double beta[]  = {0.0};
 
    Frame tip = fk.computeTipFrame(alpha, beta, 1);
 
    // Semicircle: tip at (2*radius, 0, 0) — diameter across
    double expected_x = 2.0 * radius;
    std::cout << "  Tip: " << tip << std::endl;
    std::cout << "  Expected: p=(" << expected_x << ", 0, 0)" << std::endl;
 
    assert(std::abs(tip.p.x - expected_x) < 1e-10);
    assert(std::abs(tip.p.y) < 1e-10);
    assert(std::abs(tip.p.z) < 1e-10);
 
    // Tangent should point in -z (reversed)
    Vec3 tangent = tip.tangent();
    assert(std::abs(tangent.x) < 1e-10);
    assert(std::abs(tangent.z - (-1.0)) < 1e-10);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 5: Rotate alpha changes bending plane -> 3D curve */
void test_rotation_changes_bending_plane() {
    std::cout << "=== Test: Rotation changes bending plane ===" << std::endl;

    double kappa_y = 20.0;
    double radius = 1.0/kappa_y;
    double L = (M_PI / 2.0) / kappa_y; // 90 degree bend

    std::vector<TubeParameters> tubes = {
        makeTube(2.0*mm, 1.5*mm, {{0,kappa_y,L}}, 1.0)
    };

    ForwardKinematics fk(tubes);

    // alpha = 0, bending in xz plane
    {
        double alpha[] = {0.0};
        double beta[] = {0.0};

        Frame tip = fk.computeTipFrame(alpha, beta, 1);
        std::cout << "  alpha=0:    " << tip << std::endl;
        assert(std::abs(tip.p.x - radius) < 1e-10);
        assert(std::abs(tip.p.y) < 1e-10);
    }

    // alpha = pi/2, bending in the yz plane
    {
        double alpha[] = {M_PI / 2};
        double beta[] = {0.0};
        Frame tip = fk.computeTipFrame(alpha, beta, 1);
        std::cout << "  alpha=π/2:  " << tip << std::endl;
        assert(std::abs(tip.p.x) < 1e-10);
        assert(std::abs(tip.p.y - radius) < 1e-10);
    }

    // alpha = pi: curvature flipped
    {
        double alpha[] = {M_PI};
        double beta[]  = {0.0};
        Frame tip = fk.computeTipFrame(alpha, beta, 1);
        std::cout << "  alpha=π:    " << tip << std::endl;
        assert(std::abs(tip.p.x - (-radius)) < 1e-10);
        assert(std::abs(tip.p.y) < 1e-10);
    }
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 6: Two-section backbone -> straight followed by a curved section*/
void test_two_section_backbone() {
    std::cout << "=== Test: Two-section backbone ===" << std::endl;

    double L_straight = 0.05; // 5cm
    double kappa_y = 20.0;
    double L_curved = (M_PI / 2.0) / kappa_y; // 90 degree bending
    double radius = 1.0 / kappa_y;

    std::vector<TubeParameters> tubes = {
        makeTube(
            2.0*mm, 1.5*mm, {{0, 0, L_straight}, {0, kappa_y, L_curved}}, 1.0
        )
    };

    ForwardKinematics fk(tubes);

    double alpha[] = {0.0};
    double beta[]  = {0.0};

    Frame tip = fk.computeTipFrame(alpha, beta, 1);

    // Straight section moves to (0, 0, L_straight)
    // Then 90° bend adds (radius, 0, radius) in local frame
    double expected_x = radius;
    double expected_z = L_straight + radius;

    std::cout << "  Tip: " << tip << std::endl;
    std::cout << "  Expected: p=(" << expected_x << ", 0, " << expected_z << ")" << std::endl;
 
    assert(std::abs(tip.p.x - expected_x) < 1e-10);
    assert(std::abs(tip.p.y) < 1e-10);
    assert(std::abs(tip.p.z - expected_z) < 1e-10);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 7: 4 tubes configuration in the paper */
void test_paper_4tubes_fk() {
    std::cout << "=== Test: Paper 4-tube FK (Nitinol) ===" << std::endl;

    std::vector<TubeParameters> tubes = {
        makeTubeFromMaterial(1.684*mm, 1.346*mm,
            {{0.0, 1.0/(11.89*cm), 7.58*cm}}, NITINOL_E, NITINOL_NU),
        makeTubeFromMaterial(1.295*mm, 1.036*mm,
            {{0.0, 1.0/(4.16*cm), 7.58*cm}}, NITINOL_E, NITINOL_NU),
        makeTubeFromMaterial(1.003*mm, 0.813*mm,
            {{0.0, 0.0, 7.58*cm}, {0.0, 1.0/(1.5*cm), 4.71*cm}}, NITINOL_E, NITINOL_NU),
        makeTubeFromMaterial(0.8*mm, 0.615*mm,
            {{0.0, 0.0, 13.3*cm}}, NITINOL_E, NITINOL_NU)
    };
 
    ForwardKinematics fk(tubes);
 
    // All aligned, fully inserted
    double alpha[] = {0.0, 0.0, 0.0, 0.0};
    double beta[]  = {0.0, 0.0, 0.0, 0.0};
 
    // Section frames
    auto section_frames = fk.computeSectionFrames(alpha, beta, 4);
    std::cout << "  Section frames (" << section_frames.size() << "):" << std::endl;
    for (size_t i = 0; i < section_frames.size(); i++) {
        std::cout << "    [" << i << "] " << section_frames[i] << std::endl;
    }
 
    // Base should be at origin
    assert(std::abs(section_frames[0].p.x) < 1e-12);
    assert(std::abs(section_frames[0].p.y) < 1e-12);
    assert(std::abs(section_frames[0].p.z) < 1e-12);
 
    // Tip should be the last frame
    Frame tip = fk.computeTipFrame(alpha, beta, 4);
    assert(std::abs(tip.p.x - section_frames.back().p.x) < 1e-12);
    assert(std::abs(tip.p.y - section_frames.back().p.y) < 1e-12);
    assert(std::abs(tip.p.z - section_frames.back().p.z) < 1e-12);
 
    std::cout << "\n  Tip frame: " << tip << std::endl;
 
    // Sampled frames for visualization
    auto sampled = fk.computeSampledFrames(alpha, beta, 4, 5);
    std::cout << "\n  Sampled frames (" << sampled.size() << "):" << std::endl;
    for (size_t i = 0; i < sampled.size(); i++) {
        std::cout << "    [" << i << "] " << sampled[i] << std::endl;
    }
 
    // Last sampled frame should match tip
    assert(std::abs(sampled.back().p.x - tip.p.x) < 1e-10);
    assert(std::abs(sampled.back().p.y - tip.p.y) < 1e-10);
    assert(std::abs(sampled.back().p.z - tip.p.z) < 1e-10);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

/* Test 8: 3D curve -> two tubes at different rotations */
void test_3d_curve() {
    std::cout << "=== Test: 3D curve (out-of-plane) ===" << std::endl;

    double kappa_y = 20.0;
    double L = 0.05;

    // Tube 1: curved 5cm
    // Tube 2: curved, extended 5cm beyond tube 1
    std::vector<TubeParameters> tubes = {
        makeTube(2.0*mm, 1.5*mm, {{0, kappa_y, L}}, 1.0),
        makeTube(1.5*mm, 1.0*mm, {{0, kappa_y, 2*L}}, 1.0)
    };

    ForwardKinematics fk(tubes);

    // Tube 2 rotates 90 degree from tube 1
    // section 1:resultant 45 degree between x and y
    // section 2: bends ub yz plane

    double alpha[] = {0.0, M_PI/2};
    double beta[] = {0.0, 0.0};

    Frame tip = fk.computeTipFrame(alpha, beta, 2);

    std::cout << "  Tip: " << tip << std::endl;

    // Should have nonzero x, y, and z → truly 3D
    std::cout << "  |x| = " << std::abs(tip.p.x)
              << ", |y| = " << std::abs(tip.p.y)
              << ", |z| = " << std::abs(tip.p.z) << std::endl;
    
    // The second section bends only in tube 2's plane (α=π/2)
    // so y should be nonzero
    assert(std::abs(tip.p.y) > 1e-6);
    // x should also be nonzero from the first section's combined curvature
    assert(std::abs(tip.p.x) > 1e-6);
 
    std::cout << "PASSED" << std::endl << std::endl;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    test_straight_tube();
    test_single_curved_tube();
    test_90_degree_bend();
    test_180_degree_bend();
    test_rotation_changes_bending_plane();
    test_two_section_backbone();
    test_paper_4tubes_fk();
    test_3d_curve();

    std::cout << "=== ALL FK TESTS PASSED ===" << std::endl;
    return 0;
}