#pragma once

#include "types/TubeParameters.h"
#include "types/Section.h"

#include <vector>
#include <array>

namespace ctr {
    /** 
     * Computes the constant-curvature sections of an assembled CTR needle

     * 1. Geometric decomposition
          determines section boundaries along the backbone arc 
    
     * 2. Curvature computation
          compute the value of resultant curvature
    **/

    class SectionComputer {
        public:
            explicit SectionComputer(const std::vector<TubeParameters>& tubes);

            /**
            * Compute all sections with their resultant curvatures.
            * This is the main entry point combining both tasks.
            *
            * @param alpha     Rotation angles [α₀, α₁, ..., αₙ₋₁] (radians)
            * @param beta      Insertion depths [β₀, β₁, ..., βₙ₋₁] (meters)
            * @param num_tubes Number of tubes (must match tubes vector size)
            * @return Ordered list of sections from proximal to distal
            */

            std::vector<Section> computeSections(
                const double alpha[],
                const double beta[],
                int num_tubes
            ) const;

            /**
            * Compute section boundaries only (without curvatures).
            * Useful for debugging or when you only need geometry.
            */
            std::vector<Section> computeBoundaries(
                const double beta[],
                int num_tubes) const;

            /**
            * Compute resultant curvatures for pre-computed sections.
            * Applies moment equalibrium to each section using the given α values.
            */
            void computeCurvature(
                std::vector<Section>& sections,
                const double alpha[]
            ) const;

            // Set tolerance for merging very short section
            void setMergeTolerance(double tol) {merge_tolerance_ = tol;}

            // Get the tube parameters
            const std::vector<TubeParameters>& getTubes() const {return tubes;}

        private:
            /**
            * Represents a transition point along the backbone arc length.
            * Stores what caused the transition for debugging.
            */
            struct TransitionPoint {
                double s;
                int tube_idx; 
                int type; // 0 = tube base, 1 = internal curvature change, 2 = tube tip
            };

            /**
            * Collect all transition points from tube endpoints
            * and internal curvature changes.
            */
            std::vector<TransitionPoint> collectTransitionPoints(
                const double beta[], int num_tubes
            ) const;

            /**
            * For a given backbone arc length position, determine which
            * of tube i's CurvedSection segments applies.
            *
            * @param tube_idx   Index of the tube
            * @param beta_i     Insertion depth of tube i
            * @param s          Backbone arc length position
            * @return Index into tubes[tube_idx].curved_sections,
            *         or -1 if tube is not active at position s
            */
            int getSegmentIndexAt(int tube_idx, double beta_i, double s) const;

            /**
            * Get the precurvature of tube i at backbone position s.
            * Returns the curvature from the appropriate CurvedSection segment.
            *
            * @param tube_idx  Index of the tube
            * @param beta_i    Insertion depth of tube i
            * @param s         Backbone arc length position
            * @return {κ_x, κ_y} in tube's local frame F_i
            */
            std::array<double, 2> getTubeCurvatureAt(
                int tube_idx, double beta_i, double s
            ) const;

            /**
             * Rotate a 2D curvature vector by angle a about z-axis
             * Transform from tube frame F_i to world frame W :
             * ū^W = R_z(α) · ū^F
             *  @param kappa_x  Curvature x-component in tube frame
             * @param kappa_y  Curvature y-component in tube frame
             * @param alpha    Rotation angle (radians)
             * @return {κ_x^W, κ_y^W} in world frame
             */
            
             static std::array<double, 2> rotateCurvature(
                double kappa_x, double kappa_y, double alpha
             );

             const std::vector<TubeParameters>& tubes;
             double merge_tolerance_ = 1e-8; // meters
    };
} // namespace ctr