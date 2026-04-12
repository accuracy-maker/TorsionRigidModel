#include "types/TubeParameters.h"
#include "types/Section.h"
#include "core/SectionComputer.h"

#include <iostream>
#include <cassert>
#include <cmath>

using namespace ctr {
    /**
    * Test against the 4-tube design from Sears & Dupont 2006, Table 1 and eq (11).
    *
    * Tube 1: ρ=11.89cm, single curvature, l₁=7.58cm
    * Tube 2: ρ=4.16cm,  single curvature, l₂=7.58cm
    * Tube 3: straight proximal (7.58cm) + curved distal ρ=1.5cm (4.71cm)
    * Tube 4: straight, l₄=13.3cm
    *
    * With l₁=l₂=7.58, l₃=12.29, l₄=13.3 fully inserted (β=[0,0,0,0]):
    *
    *   s₁ = max{0, l₃ - l₃₁} = 7.58   tubes 1-4, tube 3 straight
    *   s₂ = l₂ - s₁ = 0                 (zero length, should be skipped)
    *   s₃ = l₃ - l₂ = 4.71              tubes 3,4, tube 3 curved
    *   s₄ = l₄ - l₃ = 1.01              tube 4 only
    */

    void test_paper_4tube_example() {
        std::cout << "=== Test: Paper 4-tube example ===" << std::endl;

        // Convert cm to meters for consistency

        const double cm = 0.01;
        
    }
}