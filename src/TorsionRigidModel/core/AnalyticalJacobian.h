#pragma once

#include "types/TubeParameters.h"
#include "types/Section.h"
#include "types/Frame.h"
#include "SectionComputer.h"
#include "ForwardKinematics.h"

#include <vector>
#include <array>
#include <cmath>

namespace ctr {
    /*
    Analytical Jacobian for CTR using the product of matrix exponentials

    FK is W^{W(0)} (l_n) = T(u_{f1}, s1) * T(u_{f2}, s2) * ... * T(u_{fm}, sm) [implmented in ForwardKinematics.cpp]

    There is another notation about this from Dr. Leo Wu's paper

    Fk is: g = exp(ξ₁s₁) · exp(ξ₂s₂) · ... · exp(ξₘsₘ); which is exactly same with the above one.

    Therefore, they are two types of spatial Jacobian columns are:

        ∂g/∂αᵢ:  Σ_k  Ad(g₁...gₖ₋₁) · Aₖ · bₖᵢ, where sum is over sections k where tube i is active
        ∂g/∂βᵢ:  computed numerically (section boundaries are non-smooth in β)
    
    where:
        ξₖ = [κx, κy, 0, 0, 0, 1]  is the twist for section k
        gₖ = se3Exp(ξₖ · sₖ)       is the section transform
        Aₖ = aMatrix(ξₖ, sₖ)       is the integral of the Adjoint (6×6)
        bₖᵢ = [∂κₖ/∂αᵢ; 0;0;0]    is the curvature derivative (6×1)

    */

    class AnalyticalJacobian {
        public:
            explicit AnalyticalJacobian(const std::vector<TubeParameters>& tubes);
    }
}