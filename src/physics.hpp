#pragma once

#include <vector>
#include <cmath>
#include <complex>
#include "params.hpp"

// Tridiagonal matrix represented as three vectors (lower, main, upper diagonal)
struct TridiagMatrix {
    std::vector<std::complex<double>> lower;
    std::vector<std::complex<double>> main;
    std::vector<std::complex<double>> upper;
};

// Build the spatial grid from x_min to x_max with N points
std::vector<double> build_grid(const Params& p);

// Build the potential vector V(x) over the spatial grid
std::vector<double> set_potential(const std::vector<double>& x,
                                  const Params& p);

// Initialize the Gaussian wave packet psi(x, t=0)
std::vector<std::complex<double>> init_gauss_psi(const std::vector<double>& x,
                                                 const Params& p);

// Build the left-hand side tridiagonal matrix A of the Crank-Nicolson scheme
TridiagMatrix build_diagA(const std::vector<double>& V,
                          const Params& p);

// Build the right-hand side tridiagonal matrix B of the Crank-Nicolson scheme
TridiagMatrix build_diagB(const std::vector<double>& V,
                          const Params& p);