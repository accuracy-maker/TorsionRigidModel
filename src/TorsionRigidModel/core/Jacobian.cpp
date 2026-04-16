#include "Jacobian.h"

#include <cmath>
#include <cassert>
#include <cstring>

namespace ctr {
    // constructor
    Jacobian::Jacobian(const std::vector<TubeParameters>& tubes)
        : fk(tubes)
    {
    }

    // full 6 x 2n Jacobian via central finite differences
    std::vector<double> Jacobian::compute(
        const double alpha[],
        const double beta[],
        int num_tubes
    ) const {
        int n = num_tubes;
        int num_cols = 2 * n;
        int num_rows = 6;

        std::vector<double> J(num_rows * num_cols, 0.0);

        
    }
}