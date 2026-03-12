#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <Eigen/Dense>
#include "lp_oracles/solve_lp.h"   

// Wrapper to match the Point constructor expected by ComputeChebychevBall
struct VolestiPoint : public Eigen::VectorXd { 
    VolestiPoint() : Eigen::VectorXd() {}
    VolestiPoint(int d) : Eigen::VectorXd(d) {}
    
    template<typename Iter>
    VolestiPoint(int d, Iter begin, Iter end) : Eigen::VectorXd(d) {
        int i = 0;
         auto it = begin;
         while(it!=end && i<d){
             (*this)(i) = *it;
             ++it;
             ++i;
         }
    }
};

std::pair<Eigen::MatrixXd, Eigen::VectorXd> generate_random_polytope(int n, int m) {
    Eigen::MatrixXd A(m, n);
    Eigen::VectorXd b(m);

    std::random_device rd;
    std::mt19937 gen(rd());
std::uniform_real_distribution<> a_dist(-10.0, 10.0);
std::uniform_real_distribution<> factor_dist(0.5, 2.0);

    for (int i = 0; i < m; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            A(i, j) = a_dist(gen);
            row_sum += std::abs(A(i, j));
        }
        b(i) = row_sum *  factor_dist(gen); // to normalize the b(i) so it will not be very far or very near
    }
    return {A, b};
}

typedef double NT;
typedef VolestiPoint Point;

int main() {
    std::cout << "================================================================\n";
    std::cout << "   lp_solve Stress Test – Random H‑polytopes\n";
    std::cout << "================================================================\n\n";

    const int min_dim = 10;      
    const int max_dim = 200;     
    const int step = 10;         
    const int instances_per_dim = 10;   
    const double tol = 1e-6;      

    std::cout << std::left << std::setw(6) << "dim"
              << std::setw(12) << "instances"
              << std::setw(12) << "success"
              << std::setw(15) << "avg_time(ms)"
              << std::setw(15) << "avg_radius"
              << "notes\n";
    std::cout << "----------------------------------------------------------------\n";

    for (int n = min_dim; n <= max_dim; n += step) {
        int success_count = 0;
        double total_time = 0.0;
        double total_radius = 0.0;
        std::string notes = "";

        for (int inst = 0; inst < instances_per_dim; ++inst) {

            auto [A, b] = generate_random_polytope(n, 5 * n);

            try {
                auto start = std::chrono::high_resolution_clock::now();

                std::pair<Point, NT> result = ComputeChebychevBall<NT, Point, Eigen::MatrixXd, Eigen::VectorXd>(A, b);

                auto end = std::chrono::high_resolution_clock::now();
                double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

                double radius = result.second;

                bool feasible = true;
                if (radius < 0) feasible = false;
                else {
                    Eigen::VectorXd c = result.first;
                    Eigen::VectorXd residual = A * c - b; // A*c<=b -> A*c -b <=0
                    for (int i = 0; i < residual.size(); ++i) {
                        if (residual(i) > tol) { 
                            feasible = false;
                            break;
                        }
                    }
                }

                if (feasible) {
                    success_count++;
                    total_time += elapsed_ms;
                    total_radius += radius;
                } else {
                    notes = "numerical issue";
                }
            } catch (const std::exception& e) {
                notes = "exception";
            } catch (...) {
                notes = "crash";
            }
        }

        double avg_time = (success_count > 0) ? total_time / success_count : 0.0;
        double avg_radius = (success_count > 0) ? total_radius / success_count : 0.0;

        std::cout << std::left << std::setw(6) << n
                  << std::setw(12) << instances_per_dim
                  << std::setw(12) << success_count
                  << std::setw(15) << std::fixed << std::setprecision(2) << avg_time
                  << std::setw(15) << std::setprecision(6) << avg_radius
                  << notes << "\n";
    }

    return 0;
}