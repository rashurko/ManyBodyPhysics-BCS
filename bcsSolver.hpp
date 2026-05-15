#ifndef BCSSOLVER_HPP
#define BCSSOLVER_HPP

#include <vector>
#include <cmath>
#include <iostream>

struct bcsSystem {
    // Constants
    unsigned int N;                      // number of particles = number of one-particle energy levels (half-filling)
    double g;                            // pairing term
    std::vector<unsigned int> epsAlphas; // N vector of one-particle energy levels

    // Amplitudes
    std::vector<double> vAmplitudes; // N vector of amplitudes v with alpha = 1, ..., N
    std::vector<double> uAmplitudes; // N vector of amplitudes u with alpha = 1, ..., N

    double E0;       // converged ground-state energy
    double E1p;      // single-particle hamiltonian contribution
    double Epairing; // pairing contribution
};


class bcsSolver {
    private:
        bcsSystem system;
        unsigned int N;
        double g;

        bool converged = false;

        double getE0() {
            double E0 = 0.0;
            double E1p = 0.0;
            double Epairing = 0.0;
            double uvSum = 0.0;
            for (unsigned int j = 0; j < N; j++) {
                    double uj = system.uAmplitudes[j];
                    double vj = system.vAmplitudes[j];
                    uvSum += uj*vj;
                }

            for (unsigned int i = 0; i < N; i++) {
                double epsi = system.epsAlphas[i];
                double ui = system.uAmplitudes[i];
                double vi = system.vAmplitudes[i];
                E0 += (2*epsi - g*vi*vi)*vi*vi - g*ui*vi*uvSum;
                E1p += 2*epsi*vi*vi;
            }
            Epairing = E0 - E1p;
            system.E1p = E1p;
            system.Epairing = Epairing;
            return E0;
        }

        double getF(double mu, double gap) const {
            double f = 0.0;
            for (unsigned int i = 0; i < N; i++) {
                double epsi = system.epsAlphas[i];
                double vi = system.vAmplitudes[i];
                f += (sqrt(pow((epsi - g*vi*vi - mu), 2.0) + pow(gap, 2.0)) - (epsi - g*pow(vi, 2.0) - mu)) / (2*sqrt(pow((epsi - g*vi*vi - mu), 2.0) + pow(gap, 2.0)));
            }
            return f;
        }

        double bisection(double gap, double threshold) const {
            double goal = static_cast<double>(N)/2.0;

            // Determine muPos > N/2 and muNeg < N/2
            double muP = 0;
            double muN = 0;
            double fP = getF(muP, gap);
            while (fP <= goal) {
                muN = muP;
                muP += 1;
                fP = getF(muP, gap);
                if (std::abs(fP - goal) < threshold) {
                    return muP;
                }
            }

            double fN = getF(muN, gap);
            while (fN >= goal) {
                muN -= 1;
                fN = getF(muN, gap);
            }
            
            // Determine mu so that |f(mu)-N/2| < threshold, or when max number iteration is exceeded
            double mu = (muP + muN) / 2;
            double f = getF(mu, gap);
            unsigned int nIterations = 0;
            while (std::abs(f - goal) >= threshold) {
                if (f < goal) {
                    muN = mu;
                } else {
                    muP = mu;
                }
                mu = (muP + muN) / 2;
                f = getF(mu, gap);
                nIterations += 1;

                if (nIterations > 1e6) {
                    std::cerr << "Bisection method: max number of iterations reached (1e6)." << std::endl;
                    return mu;
                }
            }
            return mu;
        }

    public:
        bcsSolver(unsigned int N, double g) {
            system.N = N;
            system.g = g;
            system.E0 = 0;

            this->N = N;
            this->g = g;

            // Define one-particle energies and set initial amplitudes
            for (unsigned int i = 1; i < N + 1; i++) {
                system.epsAlphas.push_back(i - 1);
                system.vAmplitudes.push_back(sqrt(0.5));
                system.uAmplitudes.push_back(sqrt(0.5));
            }

            // Calculate the initial ground-state energy
            system.E0 = getE0();
        }

        double showE0() const {
            return system.E0;
        }

        double showE1p() const {
            return system.E1p;
        }

        double showEpairing() const {
            return system.Epairing;
        }

        bool isConverged() const {
            return converged;
        }

        void solve(double bisectionThreshold, double energyThreshold) {
            unsigned int nIterations = 0;
            while (true) {
                // 1) Calculate the gap
                double gap = 0.0;
                for (unsigned int i = 0; i < N; i++) {
                    double ui = system.uAmplitudes[i];
                    double vi = system.vAmplitudes[i];
                    gap += -g*ui*vi;
                }

                // 2) Determine mu
                double mu = bisection(gap, bisectionThreshold);

                // 3 + 4) Update one-particle energies and quasi-particle energy spectrum
                for (unsigned int i = 0; i < N; i++) {
                    double epsi = system.epsAlphas[i];
                    double vi = system.vAmplitudes[i];

                    double chii = epsi - g*pow(vi, 2.0);
                    double Ei = sqrt(pow(chii - mu, 2.0) + pow(gap, 2.0));
                    
                    double newVi = sqrt(std::abs((Ei - chii + mu) / (2*Ei)));
                    double newUi = sqrt(1 - pow(newVi, 2.0));

                    system.vAmplitudes[i] = newVi;
                    system.uAmplitudes[i] = newUi;
                }

                // Calculate new ground-state energy E0
                double newE0 = getE0();

                // Check for convergence
                if (std::abs(newE0 - system.E0) < energyThreshold) {
                    system.E0 = newE0;
                    converged = true;
                    break;
                } else if (nIterations > 1e6) {
                    std::cerr << "Convergence not reached after 1e6 steps. " << newE0 << " " << system.E0 << std::endl;
                    system.E0 = newE0;
                    break;
                }

                system.E0 = newE0;

                nIterations += 1;
            }
        }
};

#endif

