// params.hpp
// Parameter structure and function declarations for TDSE Crank-Nicolson solver
// Atomic units: hbar=1, m=1

#pragma once

#include <string>

struct Params {

    // Domain and grid
    double L;       // [bohr]      total domain length
    double dx;      // [bohr]      spatial step
    double dt;      // [hbar/Eh]   time step
    double T_tot;   // [hbar/Eh]   total simulation time

    // Derived grid quantities (computed, not read from input)
    int N;          // number of spatial points
    int M;          // number of time steps

    // Gaussian wave packet
    double x0;      // [bohr]      initial position of packet center
    double sigma;   // [bohr]      packet width
    double k0;      // [bohr^-1]   initial momentum

    // Potential barrier
    double x_bar;   // [bohr]      left edge of barrier
    double w_bar;   // [bohr]      barrier width
    double V0;      // [Eh]        barrier height
};

// Read parameters from input file and compute derived quantities
Params read_params(const std::string& filename);

// Print all parameters to stdout
void print_params(const Params& p);
