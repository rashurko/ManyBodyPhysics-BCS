#include <iostream>
#include <chrono>
#include <fstream>

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
        std::vector<double> gVals;
        std::cout << "N: ";
        std::cin >> N;
        
        char multiple;
        bool save = false;
        std::cout << "Compute for multiple g values? (y/n)";
        std::cin >> multiple;
        if (multiple == 'y') {
            // g interval
            double gMin;
            double gMax;
            double gStep;
            std::cout << "min(g): ";
            std::cin >> gMin;
            std::cout << "max(g): ";
            std::cin >> gMax;
            std::cout << "step: ";
            std::cin >> gStep;
            for (unsigned int i = 0; i <= static_cast<unsigned int>((gMax - gMin)/gStep + 0.5); i++) {
                gVals.push_back(gMin + static_cast<double>(i)*gStep);
            }

            // Save data
            char saveData;
            std::cout << "Save data? (y/n)";
            std::cin >> saveData;
            if (saveData == 'y') {
                save = true;
            }
        } else {
            double g;
            std::cout << "g: ";
            std::cin >> g;
            gVals.push_back(g);
        }
        whitespace();

        std::ofstream outFile("data.txt");
        std::ofstream outFileOccupation("occupation.txt");
        outFile << "g E0 Ep1 Epair duration" << std::endl;
        outFileOccupation << "g occupationNumbers" << std::endl;
        unsigned int count = 0;
        for (auto g : gVals) {
            count += 1;

            auto timeStart = std::chrono::high_resolution_clock::now();
            bcsSolver solver(N, g);
            ExactSolver exactSolver(N, g);
            double E0;
            double Ep1;
            double Epair;
            std::vector<double> occupationNumbers;
            switch (choice) {
                case 1:
                    solver.solve(1e-6, 1e-6);
                    E0 = solver.showE0();
                    Ep1 = solver.showE1p();
                    Epair = solver.showEpairing();
                    solver.getOccupationNumbers();
                    occupationNumbers = solver.showOccupationNumbers();
                    break;
                
                case 2:
                    exactSolver.solve();
                    E0 = exactSolver.showE0();
                    Ep1 = exactSolver.showEp1();
                    Epair = exactSolver.showEpairing();
                    exactSolver.getOccupationNumbers();
                    occupationNumbers = exactSolver.showOccupationNumbers();
                    break;
                default:
                    break;
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart);

            if (solver.isConverged() || exactSolver.isSuccess()) {
                whitespace();
                std::cout << count <<"/" << gVals.size() << std::endl;
                    if (save) {
                        outFile << g << " " << E0 << " " << Ep1 << " " << Epair << " " << duration.count() << std::endl;
                        outFileOccupation << g << " ";
                        for (auto occ : occupationNumbers) {
                            outFileOccupation << occ << " ";
                        }
                        outFileOccupation << std::endl;
                    } else {
                        std::cout << "g=" << g << " " << "E0=" << E0 << " " << "Ep1=" << Ep1 << " " << "Epair=" << Epair << " " << "duration=" << duration.count() << std::endl;
                        for (auto occ : occupationNumbers) {
                            std::cout << occ << " ";
                        }
                        std::cout << std::endl;
                    }
            } else {
                std::cout << "Calculation not converged" << std::endl;
            }
        }

        char temp;
        std::cout << "Press Enter to continue";
        std::cin.get(temp);
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    }

    whitespace();
    return 0;
}
