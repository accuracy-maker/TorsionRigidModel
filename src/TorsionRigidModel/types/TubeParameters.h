#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

namespace ctr {
    /**
     * Describe one constant curvature segment within a tube
     * One tube is allowed to have multiple segments with different curvature
     * Segments are ordered from proximal(base) to distal(tip)
    **/

    struct CurvedSection {
        double curvature_x = 0.0; // k_ix in tube ith's local frame F_i
        double curvature_y = 0.0; // k_iy in tube ith's local frame F_i
        // tau = 0 based on the rigid torsion assumption
        double length = 0.0; // arc length of the segment (unit: meter)
        
        // magnitude of curvature vector u_i = [k_ix, k_iy, 0]
        double curvatureMagnitude() const {
            return std::sqrt(curvature_x * curvature_x + curvature_y * curvature_y)
        }

        // radius of curvature rho = 1 / |u|. return infinity if it's straight
        double radiusOfCurvature() const {
            double mag = curvatureMagnitude();
            if (mag < 1e-12) return std::numeric_limits<double>::infinity();
            return 1.0/mag
        }
    };

    /**
     * Describle one single concentric tube
     * This tube is characterised by:
        - physical dimensions (diameters)
        - material properties (E, G, I, J)
        - Geometry
    
    **/

    struct TubeParameters {
        // Material Properties
        double E = 0.0; // Young's Modulus
        double I = 0.0; // Moment of Inertia
        double G = 0.0; // Shear Modulus
        double J = 0.0; // Polar moment of inertia

        // Geometry
        double outer_diameter = 0.0;
        double inner_diameter = 0.0;

        // curvature segments: from proximal to distal
        std::vector<CurvedSection> curved_sections;

        // Computed properties

        /// Total length
        double totalLength() const {
            double len = 0.0;
            for (const auto& seg: curved_sections) {
                len += seg.length;
            }
            return len;
        }

        void computeI() {
            double do4 = std::pow(outer_diameter, 4)
        }


    }
};



