#ifndef EXACTSOLVER_HPP
#define EXACTSOLVER_HPP

#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

unsigned int factorial(unsigned int val) {
    unsigned int result = 1;
    for (unsigned int i = 0; i < val - 1; i++) {
        result *= val - i;
    }
    return result;
}

struct ExactSystem {
    // Constants
    unsigned int N; // number of particles = number of one-particle energy levels (half-filling)
    double g;       // pairing term

    // One-particle energies
    std::vector<unsigned int> epsAlphas;

    // Basis vectors
    std::vector<std::vector<unsigned int>> basis;
    // Hamiltonian matrix
    unsigned int dimH;
    Eigen::MatrixXd H;

    double E0;   // ground-state energy
};

class ExactSolver {
    private:
        ExactSystem system;
        unsigned int N;
        double g;

        bool success = false;

        void generateBasis(unsigned int start, std::vector<unsigned int> &current) {
            // If vector has already max number of pairs
            if (current.size() == N/2) {
                system.basis.push_back(current);
                return;
            }

            for (unsigned int i = start; i < N; i++) {
                current.push_back(i);

                generateBasis(i+1, current);

                current.pop_back();
            }
        }

        bool connected(const std::vector<unsigned int>& a, const std::vector<unsigned int>& b) const {
            unsigned int common = 0;

            for (auto x : a) {
                for (auto y : b) {
                    if (x == y) {
                        common++;
                    }
                }
            }

            return common == (N/2 - 1);
        }

    public:
        ExactSolver(unsigned int N, double g) {
            system.N = N;
            system.g = g;
            this->N = N;
            this->g = g;

            // Declare the one-particle energies
            // Define one-particle energies and set initial amplitudes
            for (unsigned int i = 1; i < N + 1; i++) {
                system.epsAlphas.push_back(i - 1);
            }

            // Declare the basis vectors
            std::vector<unsigned int> current;
            generateBasis(0, current);

            // Declare the hamiltonian matrix
            system.dimH = factorial(N) / (factorial(N/2)*factorial(N/2));
            system.H = Eigen::MatrixXd(system.dimH, system.dimH);
            system.H.setZero();
            // Fill in the matrix
            for (unsigned int i = 0; i < system.dimH; i++) {
                for (unsigned int j = 0; j < system.dimH; j++) {
                    std::vector<unsigned int> vec1 = system.basis[i];
                    std::vector<unsigned int> vec2 = system.basis[j];

                    if (i == j) {
                        unsigned int energiesSum = 0;
                        for (unsigned int k = 0; k < N/2; k++) {
                            energiesSum += system.epsAlphas[vec1[k]];
                        }
                        system.H(i,j) = 2*static_cast<double>(energiesSum) - g*(static_cast<double>(N/2));
                    } else if (connected(vec1, vec2)) {
                        system.H(i,j) = -g;
                    }

                }
            }
        }

        void solve() {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(system.H);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Diagonalization failed" << std::endl;
                return;
            } else {
                success = true;
            }

            system.E0 = solver.eigenvalues()(0);
        }

        bool isSuccess() const {
            return success;
        }

        double showE0() const {
            return system.E0;
        }
};

#endif