#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <sstream>

namespace ctr {
    /** 
     * Represents one constant-curvature section of the assembled needle.
     * When multiple tubes are assembled concentrically, the needle's inserted length is divided
     into sections such that:
     * 1. The same set of tubes is active throughout the section
     * 2. Each active tube has constant precurvature within the section

     * The resultant curvature of each section is computed using the stiffness-weighted average
     * u_f = (Σ Kᵢ · ūᵢ^W) / (Σ Kᵢ)
     * Ki = EiIi
    **/

    struct Section {
        double s_start = 0.0; // arc length at section start
        double s_end = 0.0; // arc length at section end

        // Indices of tubes active in this section
        std::vector<int> active_tubes;

        // For each active tube, which of its CurvedSection segments
        // applies in this section. 
        std::vector<int> segment_indices;

        std::array<double, 3> resultant_curvature = [0.0, 0.0, 0.0];

        // Arc length of this section
        double length() const {
            return s_end - s_start;
        }

        double curvatureMagnitude() const {
            return std::sqrt(resultant_curvature[0] * resultant_curvature[0] + resultant_curvature[1]);
        }

        // Is this section essentially straight?
        bool isStraight(double tol = 1e-10) const {
            return curvatureMagnitude() < tol;
        }

        // Number of active tubes
        int numActiveTubes() const {
            return static_cast<int>(active_tubes.size());
        }

        // Debug string representation
        std::string toString() const {
            std::ostringstream oss;
            oss << "Section [" << s_start << ", " << s_end << "] "
                << "(len=" << length() << ") "
                << "tubes={";
            for (size_t i = 0; i < active_tubes.size(); i++) {
                if (i > 0) oss << ",";
                oss << active_tubes[i] << "(seg" << segment_indices[i] << ")";
            }
            oss << "} κ=(" << resultant_curvature[0]
                << "," << resultant_curvature[1] << ")";
            return oss.str();
        }

    };
} // namespace ctr