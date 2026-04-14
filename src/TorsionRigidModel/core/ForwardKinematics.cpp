#include "ForwardKinematics.h"

#include <cmath>
#include <cassert>

namespace ctr {
    // ============================================================
    // Constructor
    // ============================================================

    ForwardKinematics::ForwardKinematics(const std::vector<TubeParameters>& tubes)
    : sectionComputer(tubes)
    {
    }

    // ============================================================
    // Closed-form section transform: eq (14) from the paper
    // ============================================================
    Frame ForwardKinematics::sectionTransform(double kx, double ky, double s) {
        return sectionTransformAt(kx, ky, s);
    }

    Frame ForwardKinematics::sectionTransformAt(double kx, double ky, double t) {
        double n = std::sqrt(kx*kx + ky*ky); // ||u|| = curvature magnitude

        Frame T;

        if (n < 1e-12) {
            // -------------------------------------------------------
            // Straight section: ||u|| ≈ 0
            // T = [I  [0, 0, t]ᵀ]
            //     [0       1     ]
            // -------------------------------------------------------
            T.R = Mat3::identity();
            T.p = Vec3(0, 0, t);
        } else {
            // -------------------------------------------------------
            // Circular arc section: ||u|| > 0
            //
            // From eq (14) in paper, with θ = t * ||u||:
            //
            // Rotation part R (3x3):
            //   R[0][0] = (κx² + κy²cosθ) / n²
            //   R[0][1] = κxκy(1 - cosθ)  / n²
            //   R[0][2] = κy sinθ / n
            //   R[1][0] = κxκy(1 - cosθ)  / n²
            //   R[1][1] = (κy² + κx²cosθ) / n²
            //   R[1][2] = -κx sinθ / n
            //   R[2][0] = -κy sinθ / n
            //   R[2][1] = κx sinθ / n
            //   R[2][2] = cosθ
            //
            // Translation part p (3x1):
            //   p[0] =  κy(1 - cosθ) / n²
            //   p[1] = -κx(1 - cosθ) / n²
            //   p[2] = sinθ / n
            // -------------------------------------------------------
    
            double theta = t * n;
            double c = std::cos(theta);
            double s = std::sin(theta);
            double n2 = n * n;
    
            // Rotation matrix
            T.R.m[0][0] = (kx*kx + ky*ky*c) / n2;
            T.R.m[0][1] = kx*ky*(1 - c) / n2;
            T.R.m[0][2] = ky*s / n;
    
            T.R.m[1][0] = kx*ky*(1 - c) / n2;
            T.R.m[1][1] = (ky*ky + kx*kx*c) / n2;
            T.R.m[1][2] = -kx*s / n;
    
            T.R.m[2][0] = -ky*s / n;
            T.R.m[2][1] = kx*s / n;
            T.R.m[2][2] = c;
    
            // Translation
            T.p.x =  ky*(1 - c) / n2;
            T.p.y = -kx*(1 - c) / n2;
            T.p.z =  s / n;
        }

        return T;
    }

    // ============================================================
    // Compute tip frame
    // ============================================================
    Frame ForwardKinematics::computeTipFrame(const double alpha[], const double beta[], int num_tubes) const {
        // compute sections with curvature
        auto sections = sectionComputer.computeSections(alpha, beta, num_tubes);

        // chain transform
        Frame tip = Frame::identity();

        for (const auto& sec : sections) {
            double kx = sec.resultant_curvature[0];
            double ky = sec.resultant_curvature[1];
            double len = sec.length();

            Frame T = sectionTransform(kx, ky, len);

            tip = tip * T;
        }

        return tip;
    }

    // ============================================================
    // Compute frames at section boundaries
    // ============================================================
    std::vector<Frame> ForwardKinematics::computeSectionFrames(const double alpha[], const double beta[], int num_tubes) const {
        auto sections = sectionComputer.computeSections(alpha, beta, num_tubes);

        std::vector<Frame> frames;

        frames.reserve(sections.size() + 1);

        // Base Frame
        Frame current = Frame::identity();
        frames.push_back(current);

        // chain section frame
        for (const auto& sec : sections) {
            double kx = sec.resultant_curvature[0];
            double ky = sec.resultant_curvature[1];
            double len = sec.length();

            Frame T = sectionTransform(kx, ky, len);
            current = current * T;
            frames.push_back(current);
        }

        return frames;
    }

    // ============================================================
    // Compute sampled frames along backbone (for visualization)
    // ============================================================
    
    std::vector<Frame> ForwardKinematics::computeSampledFrames(
        const double alpha[],
        const double beta[],
        int num_tubes,
        int points_per_section) const
    {
        assert(points_per_section >= 1);
    
        auto sections = sectionComputer.computeSections(alpha, beta, num_tubes);
    
        std::vector<Frame> frames;
    
        // Base frame
        Frame section_base = Frame::identity();
        frames.push_back(section_base);
    
        for (const auto& sec : sections) {
            double kx = sec.resultant_curvature[0];
            double ky = sec.resultant_curvature[1];
            double len = sec.length();
            double dt = len / points_per_section;
    
            // Sample intermediate points within this section
            for (int i = 1; i <= points_per_section; i++) {
                double t = i * dt;
                Frame T_local = sectionTransformAt(kx, ky, t);
                frames.push_back(section_base * T_local);
            }
    
            // Advance section base for next section
            Frame T_full = sectionTransform(kx, ky, len);
            section_base = section_base * T_full;
        }
    
        return frames;
    }
} // namespace ctr