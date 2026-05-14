#include "bcsSolver.hpp"


int main() {
    bcsSolver solver(14, 0.1);
    solver.solve(1e-6, 1e-6);

    return 0;
}

