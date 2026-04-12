#include "SectionComputer.h"

#include <algorithm>
#include <cmath>
#include <cassert>

namespace ctr {
    // ============================================================
    // Constructor
    // ============================================================

    SectionComputer::SectionComputer(const std::vector<TubeParameters>& tubes)
        : tubes(tubes)
    {

    }

    // ============================================================
    // Main entry point
    // ============================================================
    std::vector<Section> SectionComputer::computeSections(
        const double alpha[],
        const double beta[],
        int num_tubes
    ) const {
        // step 1: compute section boundaries from beta values
        auto section = computeBoundaries(beta, num_tubes);

        // step 2: compute resultant curvatures from alpha values
        computeCurvature(section, alpha);

        return sections;
    }

    // ============================================================
    // Step 1: Geometric decomposition
    // ============================================================

    std::vector<SectionComputer::TransitionPoint>
    SectionComputer::collectTransitionPoints(const double beta[], int num_tubes) const {
        std::vector<TransitionPoint> points;

        for (int i = 0; i < num_tubes; i++) {
            double cursor = beta[i];

            // Tube base
            points.push_back({cursor, i, 0});

            // Internal curvature segment boundaries
            for (size_t j = 0; j < tubes[i].curved_sections.size(); j++) {
                cursor += tubes[i].curved_sections[j].length;

                if (j < tubes[i].curved_sections.size() - 1) {
                    points.push_back({cursor, i, 1});
                } else {
                    // tube tip
                    points.push_back({cursor, i, 2});
                }
            }
        }

        return points;
    }

    std::vector<Section> SectionComputer::computeBoundaries(const double beta[], int num_tubes) const {
        assert(num_tubes > 0);
        assert(num_tubes <= static_cast<int>(tube.size()));

        // --------------------------------------------------------
        // 1. Collect all transition points
        // --------------------------------------------------------
        auto transitions = collectTransitionPoints(beta, num_tubes);

        // --------------------------------------------------------
        // 2. Extract unique arc length values and sort
        // --------------------------------------------------------
        std::vector<double> boundaries;

        for (const auto& tp : transitions) {
            boundaries.push_back(tp.s);
        }

        std::sort(boundaries.begin(), boundaries.end());

        // Remove duplicate within merge tolerance
        std::vector<double> unique_boundaries;

        if (!boundaries.empty()) {
            unique_boundaries.push_back(boundaries[0]);
            for (size_t i = 1; i < boundaries.size(); i++) {
                if (boundaries[i] - unique_boundaries.back() > merge_tolerance) {
                    unique_boundaries.push_back(boundaries[i]);
                }
            }
        }
    }

    // --------------------------------------------------------
    // 3. Build sections between consecutive boundaries
    // --------------------------------------------------------
    std::vector<Section> sections;

    for (size_t k = 0; k + 1 < unique_boundaries.size(); k++) {
        double s_start = unique_boundaries[k];
        double s_end = unique_boundaries[k+1];

        // skip the zero-length sections
        if (s_end - s_start < merge_tolerance) {
            continue;
        }

        Section sec;
        sec.s_start = s_start;
        sec.s_end = s_end;
        
        // calculate midpoint
        double s_mid = 0.5 * (s_start + s_end);

        for (int i = 0; i < num_tubes; i++) {
            int seg_idx = getSegmentIndexAt(i, beta[i], s_mid);
            if (seg_idx >= 0) {
                sec.active_tubes.push_back(i);
                sec.segment_indices.push_back(seg_idx);
            }
        }

        // only add sections that have at least one active tube
        if (!sec.active_tubes.empty()) {
            sections.push_back(sec);
        }
    }

    // ============================================================
    // Step 2: Curvature computation (eq 9)
    // ============================================================
    void SectionComputer::computeCurvature(
        std::vector<Section>& sections,
        const double alpah[]
    ) const {
        // stiffness-weighted sum of curvature in world frame
        //  u_f^W = (Σ Kᵢ · ūᵢ^W) / (Σ Kᵢ)
        // ūᵢ^W = R_z(αᵢ) · ūᵢ^F

        double kappa_wx_sum = 0.0;
        double kappa_wy_sum = 0.0;
        double EI_sum = 0.0;

        for (size_t j = 0; j < sec.active_tubes.size(); j++) {
            int tube_idx = sec.active_tube[j];
            int seg_idx = sec.segment_indices[j];

            // Get tube's precurvature in its local frame F_i
            double kx_F = tubes[tube_idx].curved_sections[seg_idx].curvature_x;
            double ky_F = tubes[tube_idx].curved_sections[seg_idx].curvature_y;

            // Rotate into world frame W by angle alpha
            auto [kx_W, ky_W] = rotateCurvature(kx_F, ky_F, alpha[tube_idx]);

            // Accumulate stiffness-weighted curvature
            double EI_i = tubes[tube_idx].EI;
            kappa_wx_sum += EI * kx_W;
            kappa_wy_sum += EI * ky_W;
            EI_sum += EI_i;

            // Resultant curvature (torsion = 0 by assumption)
            if (EI_sum > 0.0) {
                sec.resultant_curvature[0] = kappa_wx_sum / EI_sum;
                sec.resultant_curvature[1] = kappa_wy_sum / EI_sum;
                sec.resultant_curvature[2] = 0.0;
            } else {
                sec.resultant_curvature = {0.0, 0.0, 0.0};
            }
        }
    }

    // ============================================================
    // Helper: find which segment of tube i applies at arc length s
    // ============================================================
    int SectionComputer::getSegmentIndexAt(
        int tube_idx, double beta_i, double s
    ) const {
        // convert backbone arc length s to tube-local arc length local_s = s - beta_i
        double local_s = s - beta_i;

        // check if s is outside this tube's range
        if (local_s < -merge_tolerance) {
            return -1; // before base
        }

        // walk through tube's curved sections
        double cursor = 0.0;
        const auto& segments = tubes[tube_idx].curved_sections;

        for (size_t j = 0; j < segments.size(); j++) {
            double seg_end = cursor + segments[j].length;

            if (local_s < seg_end + merge_tolerance) {
                return static_cast<int>(j);
            }

            cursor = seg_end;
        }

        // past the tube tip
        return -1;
    }

    // ============================================================
    // Helper: get tube curvature at a backbone position
    // ============================================================
    std::array<double, 2> SectionComputer::getTubeCurvatureAt(
        int tube_idx, double beta_i, double s
    ) const {
        int seg_idx = getSegmentIndexAt(tube_idx, beta_i, s);
        if (seg_idx < 0) {
            return {0.0, 0.0};
        }

        const auto& seg = tubes[tube_idx].curved_sections[seg_idx];
        return {seg.curvature_x, seg.curvature_y};
    }

    // ============================================================
    // Helper: rotate curvature from tube frame to world frame
    // ============================================================
    std::array<double, 2> SectionComputer::rotateCurvature(
        double kappa_x, double kappa_y, double alpha
    ) {
        // R_z(a) applied to [k_x, k_y]:
        // k_x^W = cos(a) * k_x - sin(a) * k_y;
        // k_y^W = sin(a) * k_x + cos(a) * k_y;

        double c = std::cos(alpha);
        double s = std::sin(alpha);

        return {
            c * kappa_x - s * kappa_y,
            s * kappa_x + c * kappa_y
        };
    }

} // namespace ctr