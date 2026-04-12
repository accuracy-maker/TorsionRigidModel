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
        double EI = 0.0;
        double G = 0.0; // Shear Modulus
        double J = 0.0; // Polar moment of inertia
        double GJ = 0.0;

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

        // Compute area moment of inertia from diameter
        // for concentric tube: I = (pi/64)*(do4 - di4)
        void computeI() {
            double do4 = std::pow(outer_diameter, 4);
            double di4 = std::pow(inner_diameter, 4);
            I = (M_PI / 64.0) * (do4 - di4);
        }

        // Compute polar moment of inertia from diameter
        // for concentric tube: J = (pi/32)*(do4 - di4)
        void computeJ() {
            double do4 = std::pow(outer_diameter, 4);
            double di4 = std::pow(inner_diameter, 4);
            J = (M_PI / 32.0) * (do4 - di4);
        }


        // Compute EI from E and I
        void computeEI() {
            EI = E * I;
        }

        // Compute GJ from G and J
        void computeGJ() {
            GJ = G * J;
        }

        // compute stiffness
        void computeStiffness(double youngs_modulus, double shear_modulus) {
            E = youngs_modulus;
            G = shear_modulus;
            computeI();
            computeJ();
            computeEI();
            computeGJ();
        }

        // compute stiffness from E and Poisson's ratio v
        // G = E / (2(1+v)) 
        void computeStiffnessFromPoisson(double youngs_modulus, double poisson_ratio) {
            double shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
            computeStiffness(youngs_modulus, shear_modulus);
        }

        std::array<double, 3> stiffnessTensor() const {
            return {EI, EI, GJ};
        }

        // Validate tube parameters
        void validate() const {
            if (outer_diameter <= 0.0)
                throw std::invalid_argument("Outer diameter must be positive");
            
            if (inner_diameter <= 0.0)
                throw std::invalid_argument("Inner diameter must be positive");

            if (inner_diameter >= outer_diameter)
                throw std::invalid_argument("Inner diameter must be less than outer diameter");

            if (EI <= 0.0)
                throw std::invalid_argument("Stiffness EI must be positive");

            if (curved_sections.empty())
                throw std::invalid_argument("Tube must have at least one curved section");

            for(const auto& seg : curved_sections) {
                if (seg.length <= 0.0)
                    throw std::invalid_argument("Section length must be positive");
            }
        }
    };

    /**
     * General factory functions for common tube configuration
     * Accepts any number of curved sections (proximal -> distal order)
     * Stiffness can be set via EI directly or computed from material properties
    **/
    inline TubeParameters makeTube(
        double outer_d, double inner_d,
        const std::vector<CurvedSection>& sections,
        double EI
    ) {
        TubeParameters tube;
        tube.outer_diameter = outer_d;
        tube.inner_diameter = inner_d;
        tube.EI = EI;
        tube.curved_sections = sections;
        return tube;
    }

    inline TubeParameters makeTubeFromMaterial(
        double outer_d, double inner_d,
        const std::vector<CurvedSection>& sections,
        double youngs_modulus, double poisson_ratio
    ) {
        TubeParameters tube;
        tube.outer_diameter = outer_d;
        tube.inner_diameter = inner_d;
        tube.curved_sections = sections;
        tube.computeStiffnessFromPoisson(youngs_modulus, poisson_ratio);

        return tube;
    }



} // namespace ctr



