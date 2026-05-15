#include <iostream>

#include "bcsSolver.hpp"
#include "exactSolver.hpp"

void whitespace() {
    std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n" << std::endl;
}

int main() {
    bool running = true;

    while (running) {
        whitespace();
        unsigned int choice;
        std::cout << "1. Calculate BCS energy" << std::endl;
        std::cout << "2. Calculate exact energy" << std::endl;
        std::cout << "3. Exit" << std::endl;

        std::cin >> choice;

        if (choice == 3) {
            running = false;
            continue;
        } else if (choice != 1 && choice != 2) {
            continue;
        }
        whitespace();

        unsigned int N;
        double g;
        std::cout << "N: ";
        std::cin >> N;
        std::cout << "g: ";
        std::cin >> g;
        whitespace();

        bcsSolver solver(N, g);
        ExactSolver exactSolver(N, g);
        double E0;
        switch (choice) {
            case 1:
                solver.solve(1e-6, 1e-6);
                E0 = solver.showE0();
                break;
            
            case 2:
                exactSolver.solve();
                E0 = exactSolver.showE0();
                break;
            default:
                break;
        }

        if (solver.isConverged() || exactSolver.isSuccess()) {
            std::cout << "E0 = " << E0 << std::endl;
        } else {
            std::cout << "Calculation not converged" << std::endl;
        }

        char temp;
        std::cout << "Press Enter to continue";
        std::cin.get(temp);
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    }

    whitespace();
    return 0;
}
