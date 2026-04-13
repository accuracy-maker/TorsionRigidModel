#include "types/TubeParameters.h"
#include "types/Section.h"
#include "core/SectionComputer.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

using namespace ctr;

//=====================================
// In the paper, they mentioned Nitinol material; So we can use some conventional E and G
//=====================================

constexpr double NITINOL_E = 58.0e9; // Young's Modulus G (Pa)
constexpr double NITINOL_NU = 0.33; // Poisson's ratio

// unit conversions
constexpr double mm = 1.0e-3;
constexpr double cm = 1.0e-2;

// Print Tube information
void printTubeInfo(const std::vector<TubeParameters>& tubes) {
    for (size_t i = 0; i < tubes.size(); i++) {
        std::cout << "  Tube " << i
                  << ": d_o=" << tubes[i].outer_diameter/mm << "mm"
                  << ", d_i=" << tubes[i].inner_diameter/mm << "mm"
                  << ", EI=" << std::scientific << std::setprecision(4) << tubes[i].EI << " N·m²"
                  << ", GJ=" << tubes[i].GJ << " N·m²"
                  << ", I=" << tubes[i].I << " m⁴"
                  << ", J=" << tubes[i].J << " m⁴"
                  << std::fixed << std::setprecision(6)
                  << std::endl;
    }
}

/**
 * Test against the 4-tube design from Sears & Dupont 2006, Table 1.
 *
 * using real Nitinol material properties:
 *   E  = 58 GPa (austenite phase)
 *   ν  = 0.33
 *   G  = E / (2(1+ν)) ≈ 21.8 GPa
 *   I  = (π/64)(d_o⁴ - d_i⁴)  from actual tube dimensions
 *   J  = 2I
 *   EI = E * I  (computed from real geometry)
 *   GJ = G * J  (stored for future torsion-compliant extension)
 */

 void test_paper_4tube_example() {
    std::cout << "=== Test: Paper 4-tube example (Nitinol) ===" << std::endl;

    std::vector<TubeParameters> tubes {
        // Tube 1: d_o = 1.684mm, d_i = 1.346mm, p = 11.89cm, l = 7.58cm
        makeTubeFromMaterial(1.684*mm, 1.346*mm, {{0.0, 1.0/(11.89*cm), 7.58*cm}}, NITINOL_E, NITINOL_NU),

        // Tube 2: d_o = 1.295mm, d_i = 1.036mm, p = 4.16cm, l = 7.58
        makeTubeFromMaterial(1.295*mm, 1.036*mm, {{0.0, 1.0/(4.16*cm), 7.58*cm}}, NITINOL_E, NITINOL_NU),

        // Tube 3: d_o = 1.003mm, d_i = 0.813mm, 
        // proximal: straight, l = 7.58cm
        // distal: p = 1.5cm, l = 4.71cm
        makeTubeFromMaterial(1.003*mm, 0.813*mm, {{0.0, 0.0, 7.58*cm}, {0.0, 1.0/(1.5*cm), 4.71*cm}}, NITINOL_E, NITINOL_NU),

        // Tube 4: d_o = 0.8mm, d_i = 0.615mm, straight, l = 13.3cm
        makeTubeFromMaterial(0.8*mm, 0.615*mm, {{0.0, 0.0, 13.3*cm}}, NITINOL_E, NITINOL_NU)
    };

    // print computed stiffness value
    printTubeInfo(tubes);

    // Verify EI decreases from outer to inner
    assert(tubes[0].EI > tubes[1].EI);
    assert(tubes[1].EI > tubes[2].EI);
    assert(tubes[2].EI > tubes[3].EI);

    // verify J = 2I for each tube
    for (const auto& tube : tubes) {
        assert(std::abs(tube.J - 2.0 * tube.I) < 1e-30);
    }

    // verify GJ/EI ratio = 1/(1+v) = 0.7519
    double expected_ratio = 1.0 / (1.0 + NITINOL_NU);
    for (const auto& tube : tubes) {
        double ratio = tube.GJ / tube.EI;
        std::cout << "  GJ/EI = " << ratio
                  << " (expected " << expected_ratio << ")" << std::endl;
        assert(std::abs(ratio - expected_ratio) < 1e-10);
    }

    // Test section computation
    SectionComputer sc(tubes);

    double alpha[] = {0.0, 0.0, 0.0, 0.0};
    double beta[] = {0.0, 0.0, 0.0, 0.0};

    auto sections = sc.computeSections(alpha, beta, 4);

    std::cout << "\n  Number of sections: " << sections.size() << std::endl;
    for (const auto& sec : sections) {
        std::cout << "  " << sec.toString() << std::endl;
    }

    // verify: there are 3 sections (s2 = 0)
    assert(sections.size() == 3);

    // section 1: s1 = 7.58cm, all 4 tubes are activated
    assert(sections[0].active_tubes.size() == 4);
    assert(std::abs(sections[0].length() - 7.58*cm) < 1e-6);

    // section 1 curvature: siffness weighted average
    // tube 3 and 4 is staright here
    double EI_1 = tubes[0].EI, EI_2 = tubes[1].EI;
    double EI_3 = tubes[2].EI, EI_4 = tubes[3].EI;
    double k1 = 1.0/(11.89*cm), k2 = 1.0/(4.16*cm);
    double expected_k = (EI_1 * k1 + EI_2 * k2) / (EI_1 + EI_2 + EI_3 + EI_4);
    std::cout << "\n  Section 1 curvature: " << sections[0].curvatureMagnitude()
              << " (expected ~" << expected_k << ")" << std::endl;
    assert(std::abs(sections[0].curvatureMagnitude() - expected_k) < 1e-6);

    // Print stiffness relationships
    std::cout << "  EI₁/EI₂ = " << EI_1/EI_2
              << " (balanced if ~1)" << std::endl;
    std::cout << "  (EI₁+EI₂)/EI₃ = " << (EI_1+EI_2)/EI_3
              << " (dominating if >>1)" << std::endl;

    // Section 2: s₃ = 4.71cm, tubes {2,3}
    assert(sections[1].active_tubes.size() == 2);
    assert(sections[1].active_tubes[0] == 2);
    assert(sections[1].active_tubes[1] == 3);
    assert(std::abs(sections[1].length() - 4.71*cm) < 1e-6);
    assert(sections[1].segment_indices[0] == 1);
 
    // Section 3: s₄ = 1.01cm, tube {3} only
    assert(sections[2].active_tubes.size() == 1);
    assert(sections[2].active_tubes[0] == 3);
    assert(std::abs(sections[2].length() - 1.01*cm) < 1e-6);
 
    std::cout << "PASSED" << std::endl << std::endl;
 }

 void test_variable_insertion() {
    // Test that changing insertion depths creates different sections
    std::cout << "=== Test: Variable insertion (s₂ ≠ 0) ===" << std::endl;

    std::vector<TubeParameters> tubes = {
        makeTubeFromMaterial(1.684*mm, 1.346*mm,
            {{0, 1.0/(11.89*cm), 7.58*cm}}, NITINOL_E, NITINOL_NU),
        makeTubeFromMaterial(1.295*mm, 1.036*mm,
            {{0, 1.0/(4.16*cm), 7.58*cm}}, NITINOL_E, NITINOL_NU),
        makeTubeFromMaterial(1.003*mm, 0.813*mm,
            {{0, 0, 7.58*cm}, {0, 1.0/(1.5*cm), 4.71*cm}}, NITINOL_E, NITINOL_NU),
        makeTubeFromMaterial(0.8*mm, 0.615*mm,
            {{0, 0, 13.3*cm}}, NITINOL_E, NITINOL_NU)
    };

    SectionComputer sc(tubes);

    // Push tube 3 forward (β₃ = -2cm), its curved section now
    // overlaps with tubes 1,2 → creates extra sections
    double alpha[] = {0.0, M_PI/4, 0.0, 0.0};
    double beta[]  = {0.0, 0.0, -2.0*cm, 0.0};

    auto sections = sc.computeSections(alpha, beta, 4);

    std::cout << "  Number of sections: " << sections.size() << std::endl;
    for (const auto& sec : sections) {
        std::cout << "  " << sec.toString() << std::endl;
    }

    assert(sections.size() >= 4);

    std::cout << "PASSED" << std::endl << std::endl;
 }


 int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Nitinol: E = " << NITINOL_E/1e9 << " GPa"
              << ", v = " << NITINOL_NU
              << ", G = " << NITINOL_E/(2*(1+NITINOL_NU))/1e9 << " GPa"
              << std::endl << std::endl;

    test_paper_4tube_example();
    test_variable_insertion();

    std::cout << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
}