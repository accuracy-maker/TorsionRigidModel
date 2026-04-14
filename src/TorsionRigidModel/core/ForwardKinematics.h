#pragma once

#include "types/TubeParameters.h"
#include "types/Section.h"
#include "types/Frame.h"
#include "SectionComputer.h"

#include <vector>

namespace ctr {
    /*
    Forward Kinematics for a concentric tube robot using torison-rigid model with piecewise constant curvature

    The FK maps the joint variables q = [a1, a2, ... , an, b1, b2, ..., bn] to the world frame W(S)

    Algorithm:
        1. SectionComputer determines sections the curvature from q
        2. For each section, compute the closed form of T(u,s) based on the paper
        3. Chain Transform: W(s) = T1 * T2 * T3 ...
    */

    class ForwardKinematics {
        public:
            /**
            * Construct with tube design parameters.
            * @param tubes  Vector of tube parameters, indexed 0 = outermost
            */
            explicit ForwardKinematics(const std::vector<TubeParameters>& tubes);

            /**
            * Compute the tip frame given joint variables.
            *
            * @param alpha     Rotation angles [α₀, ..., αₙ₋₁] (radians)
            * @param beta      Insertion depths [β₀, ..., βₙ₋₁] (meters)
            * @param num_tubes Number of tubes
            * @return Tip frame (position + orientation at the end of the backbone)
            */

            Frame computeTipFrame(
            const double alpha[],
            const double beta[],
            int num_tubes) const;

            /**
            * Compute frames at all section boundaries.
            * Returns m+1 frames for m sections (including base frame).
            *
            * @param alpha     Rotation angles (radians)
            * @param beta      Insertion depths (meters)
            * @param num_tubes Number of tubes
            * @return Vector of frames at [base, end of sec 1, ..., end of sec m]
            */
            std::vector<Frame> computeSectionFrames(
            const double alpha[],
            const double beta[],
            int num_tubes) const;

            /**
            * Compute sampled frames along the entire backbone.
            * Points are uniformly distributed by arc length within each section.
            * Used for SOFA visualization and collision meshes.
            *
            * @param alpha          Rotation angles (radians)
            * @param beta           Insertion depths (meters)
            * @param num_tubes      Number of tubes
            * @param points_per_section  Number of sample points per section
            * @return Vector of frames sampled along the backbone
            */
            std::vector<Frame> computeSampledFrames(
                const double alpha[],
                const double beta[],
                int num_tubes,
                int points_per_section = 10) const;


            /**
            * Compute the closed-form transform for a single constant-curvature section.
            * This is eq (14) from the paper: T(u, s) = matrix exponential of twist.
            *
            * For curvature u = [κx, κy, 0] and arc length s:
            *   - If |κ| ≈ 0: straight section (pure translation along z)
            *   - If |κ| > 0: circular arc with radius 1/|κ|
            *
            * @param kx  Resultant curvature x-component in world frame
            * @param ky  Resultant curvature y-component in world frame
            * @param s   Arc length of the section
            * @return Frame transform from section start to section end
            */
            static Frame sectionTransform(double kx, double ky, double s);

            /**
            * Compute a frame at arc length t within a constant-curvature section.
            * Same as sectionTransform but for an intermediate point (0 ≤ t ≤ s).
            *
            * @param kx  Resultant curvature x-component
            * @param ky  Resultant curvature y-component
            * @param t   Arc length from section start (0 ≤ t ≤ section length)
            * @return Frame at arc length t within the section
            */
            static Frame sectionTransformAt(double kx, double ky, double t);
        
            /// Access the internal SectionComputer
            const SectionComputer& getSectionComputer() const { return sectionComputer; }

        private:
            SectionComputer sectionComputer;
    };

} // namespace ctr